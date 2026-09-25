# ── Set working directory  ──────────────────────────────────────────────────
root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning"
setwd(root.path)

# ── Required packages  ────────────────────────────────────────────────────────
source("inst/libraries.R")

# ── Load functions from R folder  ────────────────────────────────────────────
source("inst/toy_example_simple/synthetic_data.R")
source("R/utils.R")
source("R/evaluation.R")

# ── General parameters  ──────────────────────────────────────────────────────────────
seed <- 1801
set.seed(seed)
VFolds <- 4 # folds to split data

n <- 10000
random_rates <- seq(0,1,0.1)
alpha <- 0.1
# ── Simulations for varying sample sizes  ─────────────────────────────────────
SL.out<- list() # list where results will be saved

# ── Synthetic data generation  ──────────────────────────────────────────────
## Training observations
exp <- generate_data(n, is_RCT = FALSE, seed = seed)
# extract observational data
SL.out$df_obs <- exp[[1]]
# extract complete data
df_complete <- exp[[2]]
summary(df_complete)

# extract optimal policy
SL.out$optimal_policy <- exp[[3]]
# extract potential outcomes
SL.out$potential_outcomes <- df_complete %>%
  select(starts_with("Potential_outcomes."))

### Test observations
exp_new_sample <- generate_data(n/2, is_RCT = FALSE, seed = seed+1)
# extract observational data
SL.out$df_new_sample <- exp_new_sample[[1]]
# extract optimal policy
SL.out$optimal_policy_new <- exp_new_sample[[3]]
# extract potential outcomes
SL.out$potential_outcomes <- exp_new_sample[[2]] %>%
  select(starts_with("Potential_outcomes."))
SL.out$prop_score_new <- exp_new_sample[[4]] 


# ── Define data parameters  ─────────────────────────────────────────────────
covariates_name <- c("X1","X2", "X3", "X4", "X5")
X <- SL.out$df_obs[,covariates_name] %>% as.matrix()
X_new <- SL.out$df_new_sample[,covariates_name] %>% as.matrix()

treatment_name <- "A" # name of treatment indicator in dataset
A <- SL.out$df_obs[,treatment_name]
A_new <- SL.out$df_new_sample[,treatment_name]

levels_A <- levels(A) # treatment levels
m <- length(levels(A)) # number of treatment levels

outcome_name <- "Y" # name of outcome in dataset
Y <- SL.out$df_obs[,outcome_name]
Y_new <- SL.out$df_new_sample[,outcome_name]

ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new)))

# ── Train conformal method  ───────────────────────────────────────────────────
# ── 0) Divide data into three even sets ─────────────────────────────────────
# ── Noisy label generation, scoring model & calibration ─────────────────────
SL.out$folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y,
                                      cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds,
                                                                                        shuffle = TRUE))
train1 <- SL.out$df_obs[c(SL.out$folds[[1]],SL.out$folds[[4]]),] # generate noisy labels
train2 <-  SL.out$df_obs[SL.out$folds[[2]],] # score model and nuisances
calibration <-  SL.out$df_obs[SL.out$folds[[3]],] # calibration

optimal_policy_cal <- SL.out$optimal_policy[SL.out$folds[[3]]]
true_potential_outcomes_cal <- df_complete[SL.out$folds[[3]],] %>% 
  select(starts_with("Potential_outcomes."))

# ── 1) Black-box label generation (i.e. estimates of (X,A*)) ───────────────────────
## 1.1) Generate random labels (i.e. A_rd)
A_rd <- apply(data.frame(1:nrow(calibration)),1,function(i)sample(as.numeric(levels_A),size=1))

## 1.2) Sample an optimal treatment per observation to generate oracular (X,A*)
SL.out$true_cal<- apply(data.frame(1:nrow(calibration)),1,function(x){
  el <- optimal_policy_cal[[x]]
  length_el <- length(el)
  if(length_el==1){el}else{el[sample(length_el,1)]}})


## 1.3) Estimate A* (OTR) using experts
# Training performed on train1
# Two predictions:
# (i) on calibration
# (ii) on SL.out$df_new_sample

X_train <- X[c(SL.out$folds[[1]],SL.out$folds[[4]]),]  #X_train <- X[SL.out$folds[[1]],] 
A_train <- A[c(SL.out$folds[[1]],SL.out$folds[[4]])]   #A_train <- A[SL.out$folds[[1]]]
Y_train <- Y[c(SL.out$folds[[1]],SL.out$folds[[4]])]   #Y_train <- Y[SL.out$folds[[1]]]

res_policies <- train_policies(SL.out$df_obs, train1, calibration, seed) # trains GLB inside as well
SL.out$doptFactorPredict_cal <- trained_policies_results$doptFactorPredict_cal
glb.model.grf <- trained_policies_results$glb.model.grf
glb.model.lm <- trained_policies_results$glb.model.lm
results.single.policy <- trained_policies_results$results.policy
results.policy.agg <- trained_policies_results$results.policy.agg
numalgs <- ncol(trained_policies_results$doptFactorPredict_cal)

unweighted_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_cal,
                                           weights =rep(1/numalgs, numalgs),
                                           df_pred = calibration,
                                           levels = as.numeric(levels_A))

unweighted_aggregation <- apply(apply(unweighted_probs, 1, function(x){
  rmultinom(1,1,prob=x)}),
  2, which.max)

SL.out$unweighted_cal <- unweighted_aggregation
SL.out$policy_cal <- SL.out$doptFactorPredict_cal[,"ql.reg.forest"]

# ── 2) Nonconformity score model (i.e. s(X,A)) ───────────────────────
# Training performed on train2
# Two predictions:
# (i) on calibration
# (ii) on SL.out$df_new_sample
SL.out$QAW.reg.train = stats::lm(formula = formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name)),
                                 data = train2)

# ── 3) Build set-valued policies ────────────────────────────────────────────
# ── 3.1) Conformal set-valued policy learning :calibration step ─────────────
# Compute margin score on calibration data 
potential_outcomes_cal <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- calibration[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  stats::predict(SL.out$QAW.reg.train, newdata = new_data)}))

margin_po <-  margin_score(potential_outcomes_cal)
# single policy noisy labels 
r0_scores_policy <- margin_po[cbind(1:nrow(calibration), SL.out$policy_cal)]
# aggregated noisy labels 
r0_scores_aggregation <- margin_po[cbind(1:nrow(calibration), SL.out$unweighted_cal)]
r1_score <- margin_po[cbind(seq_len(nrow(calibration)), A_rd)] # random labels 
r_true <- margin_po[cbind(seq_len(nrow(calibration)), SL.out$true_cal)] # true labels

# Compute margin score on new data 
potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- SL.out$df_new_sample[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  stats::predict(SL.out$QAW.reg.train, newdata = new_data)}))

margin_po_new <-  margin_score(potential_outcomes_new)
SL.out$new_scores <- margin_po_new

#### Obtain set-valued policies  ─────────────────────────────────────────────
# single policy-based noisy labels
quant <- stats::quantile(r0_scores_policy, (1-alpha))
conf_set_policy <- binary_to_confidence_set(margin_po_new < quant)
results_policy <- table.evaluation(conf_set_policy)

# aggregation-based noisy labels
quant_agg <- quantile(r0_scores_aggregation, 1 - alpha)
conf_set_agg <- binary_to_confidence_set(margin_po_new < quant_agg)
results_agg <- table.evaluation(conf_set_agg)

# ── 3.2) Greatest lower bound (GLB) ───────────────────────────────────────────
# ── grf ───────────────────────────────────────────────────────────────────────
lowers <- uppers <- matrix(0, nrow=nrow(SL.out$df_new), ncol=m)
for (l in as.numeric(levels_A)){
  data_l <- data.frame(SL.out$df_new[,covariates_name], Treatment=l)
  pred <- stats::predict(glb.model.grf, newdata = data_l, estimate.variance = TRUE)
  se <- sqrt(pred$variance.estimates)
  lowers[,l] <- pred$predictions - z * se
  uppers[,l] <- pred$predictions + z * se
}
uppest_lrw_bound <- apply(lowers, 1, max)
conf_set_grf <- binary_to_confidence_set(uppers >= uppest_lrw_bound)
results_glb_grf <- table.evaluation(conf_set_grf)

# ── glm ─────────────────────────────────────────────────────────────────────
lowers <- uppers <- matrix(0, nrow=nrow(SL.out$df_new), ncol=m)
for (l in as.numeric(levels_A)){
  data_l <- data.frame(SL.out$df_new[,covariates_name], A=factor(l, levels = levels_A))
  pred <- stats::predict(glb.model.lm, newdata = data_l, se.fit = TRUE)
  se <- sqrt(pred$se.fit)
  lowers[,l] <- (pred$fit - z * se) %>% as.numeric()
  uppers[,l] <- (pred$fit + z * se) %>% as.numeric()
}
uppest_lrw_bound <- apply(lowers, 1, max)
conf_set_lm <- binary_to_confidence_set(uppers>=uppest_lrw_bound)
results_glb_lm <- table.evaluation(conf_set_lm)

# ── 4) Collect results for this bootstrap iteration ────────────────────────
results_list <- list(
  # Exact match (coverage for correctly recommended treatment)
  exact_match = cbind(
    results_glb_lm[[1]], results_glb_grf[[1]],
    results_policy[[1]], results_agg[[1]],
    results.single.policy[[1]], results.policy.agg[[1]]
  ),
  
  # Coverage (full inclusion)
  coverage = cbind(
    results_glb_lm[[3]], results_glb_grf[[3]],
    results_policy[[3]], results_agg[[3]],
    results.single.policy[[3]], results.policy.agg[[3]]
  ),
  
  # Mean cardinality
  cardinality = c(
    results_glb_lm[[2]], results_glb_grf[[2]],
    results_policy[[2]], results_agg[[2]],
    results.single.policy[[2]], results.policy.agg[[2]]
  ),
  
  # SPV uniform
  spv_uniform = c(
    results_glb_lm[[4]], results_glb_grf[[4]],
    results_policy[[4]], results_agg[[4]],
    results.single.policy[[4]], results.policy.agg[[4]]
  ),
  
  # SPV propensity
  spv_propensity = c(
    results_glb_lm[[5]], results_glb_grf[[5]],
    results_policy[[5]], results_agg[[5]],
    results.single.policy[[5]], results.policy.agg[[5]]
))

source("inst/toy_example_simple/table.R")

# plot ECDF of nonconformity scores
data_toghether <- cbind(r0_scores_policy,r0_scores_aggregation) %>%  
  as.data.frame() %>% 
  pivot_longer(cols = everything(), 
               names_to = "Method",
               values_to = "Value")

ggplot(data_toghether, aes(x = Value, colour = Method)) +
  stat_ecdf(geom = "step", linewidth = 1, alpha=0.75) +
  geom_hline(yintercept = 1-alpha, colour = "red") +
  stat_ecdf(
    data = as.data.frame(true_score),
    aes(x = true_score,colour = "Oracular labels"),
    linetype = "dashed", colour="black",
    linewidth = 1.2
  ) +
  stat_ecdf(
    data = as.data.frame(r1_score),
    aes(x = r1_score, colour = "Random labels"),
    linetype = "dashed", colour="gray",
    linewidth = 1.2
  ) +
  labs(y = "ECDF", x = "Value", colour = "Method")

ggplot2::ggsave(filename = paste0("inst/toy_example_simple/images/ecdf_", n,".pdf"), 
                width = 10, height = 8)
