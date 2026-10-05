# ── Set working directory  ──────────────────────────────────────────────────
root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning"
setwd(root.path)

# ── Required packages  ────────────────────────────────────────────────────────
source("inst/libraries.R")

# ── Load functions from R folder  ────────────────────────────────────────────
source("R/synthetic_data.R")
source("R/utils.R")
source("R/evaluation.R")

# ── General parameters  ───────────────────────────────────────────────────────
seed <- 2026
set.seed(seed)
VFolds <- 2 # folds to split data

n <- 10000
type <- "tree"
alpha <- 0.1
z <- qnorm(1 - alpha/2)

# ── Synthetic data generation  ──────────────────────────────────────────────
## Training observations
exp <- generate_data(n, is_RCT = FALSE, seed = seed, type = type)
df_obs <- exp[[1]] # extract observational data

### Test observations
exp_new_sample <- generate_data(5000, is_RCT = FALSE, seed = seed+1, type = type)
# extract observational data
df_test <- exp_new_sample[[1]]
potential_outcomes <- exp_new_sample[[2]] %>%
  select(starts_with("Potential_outcomes."))
prop_score_new <- exp_new_sample[[4]] 

# ── Define data parameters  ─────────────────────────────────────────────────
covariates_name <- c("X1","X2", "X3", "X4", "X5")
treatment_name <- "A" # name of treatment indicator in dataset
A_new <- df_test[,treatment_name]
levels_A <- levels(A_new) # treatment levels
m <- length(levels_A) # number of treatment levels

outcome_name <- "Y" # name of outcome in dataset
Y <- df_obs[,outcome_name]
Y_new <- df_test[,outcome_name]
ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new)))

# ── Divide data in two folds ──────────────────────────────────────────────────
folds <- SuperLearner::CVFolds(n, id = NULL, Y = Y,
                               cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds,
                                                                                 shuffle = TRUE))
train1 <- df_obs[folds[[1]],] # train set-valued policy
train2 <-  df_obs[folds[[2]],] # estimating SPV

# ── 1. Train a set-valued policy ──────────────────────────────────────────────
set.seed(seed)
lowers <- uppers <- matrix(0, nrow=nrow(df_test), ncol=m)
if(type=="tree"){
  glb.model.grf <- grf::regression_forest(
    X = cbind(train1[,covariates_name], as.numeric(train1[,treatment_name])), 
    Y = train1[,outcome_name], seed = seed)
  
  for (l in as.numeric(levels_A)){
    data_l <- data.frame(df_test[,covariates_name], A=l)
    pred <- stats::predict(glb.model.grf, newdata = data_l, estimate.variance = TRUE)
    se <- sqrt(pred$variance.estimates)
    lowers[,l] <- pred$predictions - z * se
    uppers[,l] <- pred$predictions + z * se
  }
  uppest_lrw_bound <- apply(lowers, 1, max)
}else{
  formula_lm <- stats::as.formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name))
  glb.model.lm <- stats::lm(formula = formula_lm, data = train1)
  
  for (l in as.numeric(levels_A)){
    data_l <- data.frame(df_test[,covariates_name], A=factor(l, levels = levels_A))
    pred <- stats::predict(glb.model.lm, newdata = data_l, se.fit = TRUE)
    se <- pred$se.fit
    lowers[,l] <- (pred$fit - z * se) %>% as.numeric()
    uppers[,l] <- (pred$fit + z * se) %>% as.numeric()
  }
  uppest_lrw_bound <- apply(lowers, 1, max)
}
conf_set_lm <- binary_to_confidence_set(uppers >= uppest_lrw_bound)

# ── 2. SPV estimation ─────────────────────────────────────────────────────────
### Oracular
oracular_SPV <- set_policy_value_plug_in(test_set = conf_set_lm, 
                                         test = df_test, 
                                         Q.all.actions = potential_outcomes, 
                                         gAX.pred = prop_score_new, 
                                         levels = levels_A)


# base arguments for SPV functions
base_args <- list(
  test_set = conf_set_lm,
  test     = df_test, levels = levels_A,
  Y        = Y_new,A = A_new, ab = ab) 

# SL library 
SL.library_cond <- c("SL.randomForest", "SL.ksvm", "SL.mean", "SL.glm", "SL.xgboost")

# ── Train outcome models (Q) ────────────────────────────────────────────────
# Correctly specified 
potential_outcomes_new <- potential_outcomes
  
# Mis-specified 
QAW.reg.train_misspecified <- SuperLearner::SuperLearner(
    Y = train2[, outcome_name], 
    X = train2[, c("X2", "X3", "X4", "X5", treatment_name)],
    SL.library = SL.library_cond, family = "gaussian") 

potential_outcomes_new_misspecified <- do.call(cbind,map(seq_len(m), function(val) {
  d <- df_test[, c("X2", "X3", "X4", "X5", treatment_name)]
  d[[treatment_name]] <- factor(val, levels = levels_A)
  SuperLearner::predict.SuperLearner(QAW.reg.train_misspecified, newdata = d)$pred
}))
  
  
# ── Train propensity score models (g) ─────────────────────────────────────────
# Correctly specified 
gAX.pred <- prop_score_new
  
# Mis-specified
if(type=="tree"){
  approx <- function(t,p){
    stopifnot(p>=0 & p<= 1) 
    sign(t)*(abs(t)^p)
  }
  p<- 0
  sx1_train2 <- approx(t=train2[,"X1"], p=p)
  gAX.train_misspecified <- grf::probability_forest(
    X = cbind(sx1_train2, train2[, c("X2", "X3", "X4", "X5")]), 
    Y = as.factor(train2[, treatment_name]))
  
  sx1_test <- approx(t=df_test[,"X1"], p=p)
  gAX.pred_misspecified <- stats::predict(gAX.train_misspecified, 
                                          newdata = cbind(sx1_test, 
                                                          df_test[, c("X2","X3", "X4", "X5")]))$pred
  
}else{
  gAX.train_misspecified <- grf::probability_forest(
    X =train2[, c("X2", "X3", "X4", "X5")], 
    Y = as.factor(train2[, treatment_name]))
  
  gAX.pred_misspecified <- stats::predict(gAX.train_misspecified, 
                                          newdata = df_test[, c("X2","X3", "X4", "X5")])$pred
  
}
# ── Define Model Combinations ───────────────────────────────────────
models <- list(
  correct   = list(Q = potential_outcomes_new, 
                     g = gAX.pred),
  Q_miss    = list(Q = potential_outcomes_new_misspecified, 
                     g = gAX.pred),
  g_miss    = list(Q = potential_outcomes_new, 
                     g = gAX.pred_misspecified),
  both_miss = list(Q = potential_outcomes_new_misspecified, 
                     g = gAX.pred_misspecified))
  
specs <- list(
  #plug_in   = list(fn = set_policy_value_plug_in, 
  #                   Q = models$correct$Q, g = models$correct$g),
  AIPW      = list(fn = set_policy_value_aipw,    
                     Q = models$correct$Q, g = models$correct$g),
  TMLE      = list(fn = set_policy_value_tmle,    
                     Q = models$correct$Q, g = models$correct$g),
  AIPW_QAX  = list(fn = set_policy_value_aipw,    
                     Q = models$Q_miss$Q,  g = models$Q_miss$g),
  TMLE_QAX  = list(fn = set_policy_value_tmle,    
                     Q = models$Q_miss$Q,  g = models$Q_miss$g),
  AIPW_gAX  = list(fn = set_policy_value_aipw,    
                     Q = models$g_miss$Q,  g = models$g_miss$g),
  TMLE_gAX  = list(fn = set_policy_value_tmle,    
                     Q = models$g_miss$Q,  g = models$g_miss$g),
  AIPW_both = list(fn = set_policy_value_aipw,    
                     Q = models$both_miss$Q, g = models$both_miss$g),
  TMLE_both = list(fn = set_policy_value_tmle,    
                     Q = models$both_miss$Q, g = models$both_miss$g))
  
# ── Evaluation ────────────────────────────────────────────────────────────────
plot_data <- imap_dfr(specs, function(s, name) {
  args <- c(base_args, list(Q.all.actions = s$Q, gAX.pred = s$g))
  res  <- do.call(s$fn, args[intersect(names(args), formalArgs(s$fn))])
  

  extract_df <- function(item, df_name) {
    val   <- as.numeric(item)
    lower <- attr(item, "low") %||% attr(item, "conf.low") %||% NA
    upper <- attr(item, "high") %||% attr(item, "conf.high") %||% NA
    
    data.frame(
      estimator = name,
      SPV = df_name,
      value = val,
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      stringsAsFactors = FALSE
    )}
  
  bind_rows(
    extract_df(res[[1]], "Uniform SPV"),
    extract_df(res[[2]], "Propensity SPV")
  )
})

saveRDS(plot_data, 
        file =  paste0("inst/toy_examples/images_", type,"/dr_ci_oracle.rds"))  


# ── Plot SPVs ─────────────────────────────────────────────────────────────────
hline_data <- data.frame(
  SPV          = c("Uniform SPV", "Propensity SPV"),
  target_value = c(oracular_SPV[[1]], oracular_SPV[[2]])
)


p <- plot_data |>
  mutate(
    estimator = factor(
      estimator,
      levels = c(
        "plug_in",
        "AIPW", "TMLE",
        "AIPW_QAX", "TMLE_QAX",
        "AIPW_gAX", "TMLE_gAX",
        "AIPW_both", "TMLE_both"
      ),
      labels = c(
        "Plug-in",
        "AIPW", "TMLE",
        "AIPW (Q-model mis-specified)", "TMLE (Q-model mis-specified)",
        "AIPW (g-model mis-specified)", "TMLE (g-model mis-specified)",
        "AIPW (Both mis-specified)", "TMLE (Both mis-specified)"
      )
    )
  ) |>
  ggplot(aes(x = estimator, y = value, ymin = lower, ymax = upper, color = estimator)) +
  geom_point(
    data = \(.df) filter(.df, estimator == "Plug-in"), shape=16
  ) +
  geom_errorbar() +
  geom_hline(
    data = hline_data,
    aes(yintercept = target_value), linetype = "dashed", color = "black") +
  facet_wrap(~ SPV) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "Estimator", y = "Value")

ggsave(p, 
       filename =  paste0("inst/toy_examples/images_", type,"/DR_estimators_CIs_oracle.pdf"), 
       width = 12, height = 5)
