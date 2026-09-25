# ── Set working directory  ──────────────────────────────────────────────────
root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning"
setwd(root.path)

# ── Required packages  ────────────────────────────────────────────────────────
library(SL.ODTR)
library(hitandrun)
library(tidyr)
library(dplyr)
library(lava)
library(purrr)
library(grf)
library(randomForest)
library(gridExtra)
library(SuperLearner)
library(policytree)
library(glmnet)
library(tmle)
library(parallel)
library(caret)
library(polle)
library(viridisLite)

# ── Load functions from R folder  ────────────────────────────────────────────
source("R/synthetic_data.R")
source("R/label-estimation.R")
source("R/utils.R")
source("R/evaluation.R")

# ── General parameters  ──────────────────────────────────────────────────────────────
seed <- 1801
set.seed(seed)
VFolds <- 4 # folds to split data
synthetic_scenario <- TRUE
type <- "complex" # additional name for images (here: type of synthetic scenario)

is_RCT <- ifelse(type=="normal", FALSE, TRUE)
RCT_file<- ifelse(is_RCT==TRUE,"RCT/", "non_RCT/")
n <- 6000
ncov <- 4 #ifelse(type=="normal",4, 2)
n_test <- 100
alpha <- 0.1

# ── Simulations for varying sample sizes  ─────────────────────────────────────
SL.out<- list() # list where results will be saved
  
# ── Synthetic data generation  ──────────────────────────────────────────────
## Training observations
exp <- generate_data(n, ncov = ncov, type=type, is_RCT=is_RCT, seed = seed)
# extract observational data
SL.out$df_obs <- exp[[1]]
# extract complete data
df_complete <- exp[[2]]
# extract optimal policy
SL.out$optimal_policy <- exp[[3]]
# extract potential outcomes
SL.out$potential_outcomes <- df_complete %>%
    select(starts_with("Potential_outcomes."))
  
### Test observations
exp_new_sample <- generate_data(n/2, ncov=ncov, type=type)
# extract observational data
SL.out$df_new_sample <- exp_new_sample[[1]]
# extract optimal policy
SL.out$optimal_policy_new <- exp_new_sample[[3]]
# extract potential outcomes
SL.out$potential_outcomes <- exp_new_sample[[2]] %>%
    select(starts_with("Potential_outcomes."))
  
# ── Define data parameters  ─────────────────────────────────────────────────
covariates_name <- c("X1","X2", "X3", "X4")
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
  
# ── 0) Divide data into three even sets ─────────────────────────────────────
# ── Noisy label generation, scoring model & calibration ─────────────────────
SL.out$folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y,
                                        cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds,
                                                                                          shuffle = TRUE))
  
train1 <- SL.out$df_obs[c(SL.out$folds[[1]],SL.out$folds[[4]]),] # generate noisy labels
train2 <-  SL.out$df_obs[SL.out$folds[[2]],] # score model and nuisances
test <-  SL.out$df_obs[SL.out$folds[[3]],] # calibration
optimal_policy_test <- SL.out$optimal_policy[SL.out$folds[[3]]]
true_potential_outcomes_test <- df_complete[SL.out$folds[[3]],] %>% 
    select(starts_with("Potential_outcomes."))
  
# ── 1) Black-box label generation (i.e. estimates of (X,A*)) ───────────────────────
## 1.1) Generate random labels (i.e. A_rd)
A_rd <- apply(data.frame(1:nrow(test)),1,function(i)sample(as.numeric(levels_A),size=1))
  
## 1.2) Sample an optimal treatment per observation to generate oracular (X,A*)
SL.out$true_cal<- apply(data.frame(1:nrow(test)),1,function(x){
  el <- optimal_policy_test[[x]]
  length_el <- length(el)
  if(length_el==1){el}else{el[sample(length_el,1)]}})
  
## Learn the treatment assignment mechanism by doctor's. 
#g.reg.train_spv <- grf::probability_forest(X = X[c(SL.out$folds[[1]], SL.out$folds[[4]]),], 
#                                                  Y =  A[c(SL.out$folds[[1]],SL.out$folds[[4]])] %>% 
#                                                           as.factor())

#rf.reg.train_spv <- randomForest::randomForest(x = train1[,covariates_name],y = train1[,treatment_name])
#glm.reg.train_spv <- stats::glm(paste0(treatment_name, "~."),data = train1[,c(covariates_name, treatment_name)], family = "binomial")

## 1.3) Estimate A* (OTR) using experts
# Training performed on train1
# Two predictions:
# (i) on test
# (ii) on SL.out$df_new_sample
 
X_train <- X[c(SL.out$folds[[1]],SL.out$folds[[4]]),]  #X_train <- X[SL.out$folds[[1]],] 
A_train <- A[c(SL.out$folds[[1]],SL.out$folds[[4]])]   #A_train <- A[SL.out$folds[[1]]]
Y_train <- Y[c(SL.out$folds[[1]],SL.out$folds[[4]])]   #Y_train <- Y[SL.out$folds[[1]]]

source("inst/toy_example/train_policies.R")

unweighted_probs_top <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_test[,c(2,6)],
                                           weights =rep(1/2, 2),
                                           df_pred = test,
                                           levels = as.numeric(levels_A))

unweighted_aggregation_top <- apply(apply(unweighted_probs_top, 1, function(x){
  rmultinom(1,1,prob=x)}),
  2, which.max)

unweighted_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_test,
                                           weights =rep(1/numalgs, numalgs),
                                           df_pred = test,
                                           levels = as.numeric(levels_A))
unweighted_aggregation <- apply(apply(unweighted_probs, 1, function(x){
  rmultinom(1,1,prob=x)}),
  2, which.max)
  
unweighted_probs_new <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_new,
                                           weights =rep(1/numalgs, numalgs),
                                           df_pred = SL.out$df_new_sample,
                                           levels = as.numeric(levels_A))
unweighted_aggregation_new <- apply(apply(unweighted_probs_new, 1, function(x){
  rmultinom(1,1,prob=x)}),
  2, which.max)

source("inst/toy_example/train_policies_weighted.R")


# ── 2) Nonconformity score model (i.e. s(X,A)) ───────────────────────
# Training performed on train2
# Two predictions:
# (i) on test
# (ii) on SL.out$df_new_sample
SL.out$glm = stats::lm(formula = formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name)),
                        data = train2)

potential_outcomes_test_glm <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- test[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  stats::predict(SL.out$glm, newdata = new_data)}))

potential_outcomes_new_glm <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- SL.out$df_new_sample[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  stats::predict(SL.out$glm, newdata = new_data)}))

SL.out$grf <- grf::regression_forest(X = cbind(
  X[SL.out$folds[[2]],], A[SL.out$folds[[2]]]), 
  Y = Y[SL.out$folds[[2]]], seed = seed)

potential_outcomes_test_grf <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- cbind(X[SL.out$folds[[3]],], factor(val, levels=levels_A) %>% as.numeric())
  stats::predict(SL.out$grf, newdata = new_data)$predictions}))

potential_outcomes_new_grf <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- cbind(X_new, factor(val, levels=levels_A) %>% as.numeric())
  stats::predict(SL.out$grf, newdata = new_data)$predictions}))


###### Accuracy of policy learning and ncs methods 
dopt_test <- cbind(SL.out$doptFactorPredict_test, #glm_ncs, grf_ncs, 
                   unweighted_aggregation_top, 
                   unweighted_aggregation, 
                   weighted_aggregation)

# Oracular policy value 
oracular_pv <- apply(dopt_test %>% as.data.frame(), 2, 
                     function(x)oracular_set_policy_value(test_set= x, 
                                                          test= test, 
                                                          test_potential_outcome = true_potential_outcomes_test,
                                                          covariates = covariates_name))

# Inclusion accuracy wrt oracular 
proportions_test <- apply(dopt_test, 2, function(col) {
  mean(mapply(function(x, ref) x %in% ref, col, optimal_policy_test))
})

g.reg.train <- grf::probability_forest(X = train2[,covariates_name],
                                              Y = train2[,treatment_name])

rf.reg.train <- randomForest::randomForest(x = train2[,covariates_name], 
                                           y = train2[,treatment_name])
  
glm.reg.train <- glmnet::glmnet(x = train2[,covariates_name], 
                                y = train2[,treatment_name],
                                lambda= 0, family = "multinomial")
# 2.2) Predict nonconformity scores (margin score)
evaluate_technique <- function(scores, true_scores, random_scores, alpha = 0.1) {
  # Threshold at 1 - alpha for this technique
  q_val <- stats::quantile(scores, probs = 1 - alpha, na.rm = TRUE)
  
  # Evaluate ECDFs at q_val
  f_hat  <- stats::ecdf(scores)(q_val)
  f_true <- stats::ecdf(true_scores)(q_val)
  f_rd   <- stats::ecdf(random_scores)(q_val)
  
  # Compute r_bar
  r_bar <- (f_hat - f_true) / (f_hat - f_rd)
  
  # Flag violations
  violates <- r_bar > 1
  
  reason <- case_when(
    !violates     ~ "Valid (r_bar <= 1)",
    f_rd > f_true ~ "Violation: Random outperforms Oracular (F_rd > F_true)",
    f_hat < f_rd  ~ "Violation: Technique worse than Random (F_hat < F_rd)",
    TRUE          ~ "Violation: Unknown reason"
  )
  
  tibble(
    Quantile_Threshold = q_val,
    F_hat  = f_hat,
    F_true = f_true,
    F_rd   = f_rd,
    r_bar  = r_bar,
    Violates_Constraint = violates,
    Status = reason
  )
}
# Nonconformity scores on calibration data
################################################################################
###################################  GLM ####################################### 
################################################################################
margin_po_glm <-  margin_score(potential_outcomes_test_glm)
margin_po_glm_new <-  margin_score(potential_outcomes_new_glm)

# Noisy nonconformity scores
r0_scores_glm <- apply(dopt_test, 2, function(x){margin_po_glm[cbind(1:nrow(test), x)]})
r0_scores_glm_new <- apply(SL.out$doptFactorPredict_new, 2, function(x){margin_po_glm_new[cbind(1:nrow(SL.out$df_new_sample), x)]})
# Nonconformity scores on oracular data
true_score_glm <- margin_po_glm[cbind(1:nrow(test), SL.out$true_cal)]
# Nonconformity scores on random labels
r1_glm <- margin_po_glm[cbind(1:nrow(test), A_rd)]


data_toghether_glm <- r0_scores_glm %>%  
  as.data.frame() %>% 
  pivot_longer(cols = everything(), 
               names_to = "Method",
               values_to = "Value")

ggplot(data_toghether_glm, aes(x = Value, colour = Method)) +
  stat_ecdf(geom = "step", linewidth = 1, alpha=0.75) +
  geom_hline(yintercept = 1-alpha, colour = "red") +
  stat_ecdf(
    data = as.data.frame(true_score_glm),
    aes(x = true_score_glm,colour = "Oracular labels"),
    linetype = "dashed", colour="black",
    linewidth = 1.2
  ) +
  stat_ecdf(
    data = as.data.frame(r1_glm),
    aes(x = r1_glm, colour = "Random labels"),
    linetype = "dashed", colour="gray",
    linewidth = 1.2
  ) +
  labs(y = "ECDF", x = "Value", colour = "Method")

ggplot2::ggsave(filename = paste0("inst/images/tryzone/",type,"/glm_ecdf_",type, "_", n,".pdf"), 
                width = 10, height = 8)


# 3. Apply across all techniques in r0_scores_glm
report_summary_glm <- r0_scores_glm %>%
  as.data.frame() %>%
  map_dfr(evaluate_technique, alpha = alpha, 
          true_scores = true_score_glm, 
          random_scores = r1_glm, .id = "Method")


saveRDS(report_summary_glm, file = paste0("inst/images/tryzone/",type,"/report_glm_",type, "_", n,".rds"))


results_glm <- t(rbind(#glm.ncs_accuracy = proportions_test_glm, 
                       oracle_accuracy = proportions_test, 
                       oracular_spv = oracular_pv, 
                       r_bar = report_summary_glm$r_bar 
                       #,problematic.cases= prop_glm
                       )) %>% as.data.frame()



df_long_glm <- results_glm %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Method")%>% 
  pivot_longer(cols = c(
    #"glm.ncs_accuracy", 
    "oracle_accuracy", "oracular_spv", "r_bar", #"problematic.cases"
    ), 
               names_to = "Metric", 
               values_to = "Value")

df_long_glm$Method <- factor(df_long_glm$Method)
method_levels <- levels(df_long_glm$Method)

method_signs_glm <- df_long_glm %>%
  filter(Metric == "r_bar") %>%
  select(Method, Value) %>%
  mutate(Method = factor(Method, levels = method_levels)) %>%
  arrange(Method) %>%
  mutate(
    Color = case_when(
      Value > 1 ~ "red",
      Value < 0 ~ "green",
      TRUE      ~ "black"
    )
  )

y_label_colors_glm <- setNames(method_signs_glm$Color, method_signs_glm$Method)

# 2. Plot with automatic legend mapping
ggplot(df_long_glm) +
  geom_vline(xintercept = c(0, 1)) +
  geom_point(aes(x = Value,y = Method), size = 3, shape = 3, color = "steelblue") +
  facet_grid(cols = vars(Metric), scales = "free_x") +
  scale_x_continuous(limits = function(x) {
    # If the metric's data goes negative, force c(-1, 1), otherwise c(0, 1)
    #if (min(x, na.rm = TRUE) < -1 ) {
    #  return(c(-1, 1))
    #} else 
    if((max(x, na.rm = TRUE) <=1)&(min(x, na.rm = TRUE)>=0)) {
      return(c(0, 1))
    }else{
      return(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
    }}) +
  #scale_x_continuous(limits = function(x) c(quantile(x, 0.001, na.rm=TRUE), quantile(x, 1, na.rm=TRUE))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(color = y_label_colors_glm)) +
  labs(
    x = "",
    y = "Noisy label generation technique",
    title = "Benchmarking Noisy Label Generation Using GLM-Derived Nonconformity Scores"
  )

ggplot2::ggsave(filename = paste0("inst/images/tryzone/",type,"/glm_benchmark_",type, "_", n,".pdf"), 
                width = 12, height = 8)


################################################################################
###################################  GRF ####################################### 
################################################################################
margin_po_grf <-  margin_score(potential_outcomes_test_grf)
margin_po_grf_new <-  margin_score(potential_outcomes_new_grf)

# Noisy nonconformity scores
r0_scores_grf <- apply(dopt_test, 2, function(x){margin_po_grf[cbind(1:nrow(test), x)]})
r0_scores_grf_new <- apply(SL.out$doptFactorPredict_new, 2, function(x){margin_po_grf_new[cbind(1:nrow(SL.out$df_new_sample), x)]})

# Nonconformity scores on oracular data
true_score_grf <- margin_po_grf[cbind(1:nrow(test), SL.out$true_cal)]
# Nonconformity scores on random labels
r1_grf <- margin_po_grf[cbind(1:nrow(test), A_rd)]

#prop_grf <- colSums(incorrect&(r0_scores_grf < true_score_grf))/colSums(incorrect)

data_toghether_grf <- r0_scores_grf %>%  
  as.data.frame() %>% 
  pivot_longer(cols = everything(), 
               names_to = "Method",
               values_to = "Value")

ggplot(data_toghether_grf, aes(x = Value, colour = Method)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  geom_hline(yintercept = 1-alpha, colour = "red") +
  stat_ecdf(
    data = as.data.frame(true_score_grf),
    aes(x = true_score_grf,colour = "Oracular labels"),
    linetype = "dashed", colour="black",
    linewidth = 1.2
  ) +
  stat_ecdf(
    data = as.data.frame(r1_grf),
    aes(x = r1_grf, colour = "Random labels"),
    linetype = "dashed", colour="gray",
    linewidth = 1.2
  ) +
  labs(y = "ECDF", x = "Value", colour = "Method")
ggplot2::ggsave(filename = paste0("inst/images/tryzone/",type,"/grf_ecdf_",type, "_", n,".pdf"), 
                width = 10, height = 8)

# 3. Apply across all techniques in r0_scores_glm
report_summary_grf <- r0_scores_grf %>%
  as.data.frame() %>%
  map_dfr(evaluate_technique, alpha = alpha, 
          true_scores = true_score_grf, 
          random_scores = r1_grf, .id = "Method")

saveRDS(report_summary_grf, file = paste0("inst/images/tryzone/",type,"/report_grf_",type, "_", n,".rds"))
# proportion of inclusion
results_grf <- t(rbind(#grf.ncs_accuracy = proportions_test_grf, 
                     oracle_accuracy = proportions_test, 
                     oracular_spv = oracular_pv, 
                     r_bar = report_summary_grf$r_bar 
                     #, problematic.cases = prop_grf
                     )) %>% as.data.frame()

df_long_grf <- results_grf %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Method")%>% 
  pivot_longer(cols = c(
    #"grf.ncs_accuracy", 
    "oracle_accuracy", "oracular_spv", "r_bar"#, "problematic.cases"
    ), 
               names_to = "Metric", 
               values_to = "Value")

df_long_grf$Method <- factor(df_long_grf$Method)
method_levels <- levels(df_long_grf$Method)

method_signs_grf <- df_long_grf %>%
  filter(Metric == "r_bar") %>%
  select(Method, Value) %>%
  mutate(Method = factor(Method, levels = method_levels)) %>%
  arrange(Method) %>%
  mutate(
    Color = case_when(
      Value > 1 ~ "red",
      Value < 0 ~ "green",
      TRUE      ~ "black"
    )
  )

y_label_colors_grf <- setNames(method_signs_grf$Color, 
                               method_signs_grf$Method)

# 2. Plot with automatic legend mapping
ggplot(df_long_grf) +
  geom_vline(xintercept = c(0, 1)) +
  geom_point(aes(x = Value,y = Method), size = 3, shape = 3, color = "steelblue") +
  facet_grid(cols = vars(Metric), scales = "free_x") +
  scale_x_continuous(limits = function(x) {
    # If the metric's data goes negative, force c(-1, 1), otherwise c(0, 1)
    #if (min(x, na.rm = TRUE) < -1 ) {
    #  return(c(-1, 1))
    #} else 
    if((max(x, na.rm = TRUE) <=1)&(min(x, na.rm = TRUE)>=0)) {
      return(c(0, 1))
    }else{
      return(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
    }}) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(color= y_label_colors_grf))+
  labs(
    x = "",
    y = "Noisy label generation technique",
    title = "Benchmarking Noisy Label Generation Using GRF-Derived Nonconformity Scores"
  )

ggplot2::ggsave(filename = paste0("inst/images/tryzone/",type,"/grf_benchmark_",type, "_", n,".pdf"), 
                width = 12, height = 8)

if(type=="normal"){
  beta_low_vec  <- c(10,10,1,3,2)  # 1,2 high
  beta_high_vec <- c(2,1,10,10,4)  # 3,4 low 
  w <- stats::plogis(X[,1] + X[,2] - 0.5)
  beta <- (1 - w) * matrix(beta_low_vec,  nrow=nrow(X), ncol=5, byrow=TRUE) +
    w * matrix(beta_high_vec, nrow=nrow(X), ncol=5, byrow=TRUE)
  #epsilon <- matrix(stats::rnorm(nrow(X) * 5), nrow=nrow(X), ncol=5)
  
  treatment_assignment <- beta #+ epsilon
  
  probs <- exp(treatment_assignment - apply(treatment_assignment, 1, max))
  expit_treatment <- probs / rowSums(probs)
  expit_treatment <- as.data.frame(expit_treatment)
  colnames(expit_treatment) <- levels_A
}else{
  beta_high <- c(4, 10, 4, 2, 4)  # treatment 2 and then a bit of 1 and 4
  beta_medium <- c(10, 2, 2, 10, 15) # treatment 4, 1 and 5 and a bit of the others 
  beta_low <- c(4, 2, 10, 4, 2)  # treatment 3 and then a bit of 1 and 4
  z_axis <- X[,1] - X[,2]
  w_low  <- stats::plogis(5*(z_axis - 0.25))
  w_high  <- stats::plogis(5*(-z_axis - 0.5))
  w_medium <- pmax(0, 1 - (w_high + w_low))
  total_w  <- w_low + w_medium + w_high
  
  beta <- (w_low/total_w) * matrix(beta_low,  nrow=nrow(X), ncol=5, byrow=TRUE) +
    (w_high/total_w) * matrix(beta_high, nrow=nrow(X), ncol=5, byrow=TRUE) + 
    (w_medium/total_w) * matrix(beta_medium, nrow=nrow(X), ncol=5, byrow=TRUE)
  
  #epsilon <- matrix(stats::rnorm(nrow(X) * 5), nrow=nrow(X), ncol=5)
  treatment_assignment <- beta #+ epsilon
  probs <- exp(treatment_assignment - apply(treatment_assignment, 1, max))
  expit_treatment <- probs / rowSums(probs)
  expit_treatment <- as.data.frame(expit_treatment)
  colnames(expit_treatment) <- levels_A
}


grf.ps <- stats::predict(g.reg.train, X[SL.out$folds[[3]],])$predictions %>% 
  as.data.frame()
rf.ps <- stats::predict(rf.reg.train, X[SL.out$folds[[3]],], type = "prob")%>% 
  as.data.frame()
glm.ps <- stats::predict(glm.reg.train, 
                         newx = X[SL.out$folds[[3]],], 
                         type="response")[,,1] %>% as.data.frame()

combined_df <- bind_rows(
  grf.ps = grf.ps, 
  rf.ps = rf.ps,
  glm.ps = glm.ps, 
  .id = "Dataset"
) %>%
  pivot_longer(
    cols = -Dataset,
    names_to = "Column",
    values_to = "Value"
  )

ggplot(combined_df, aes(x = Value, colour = Dataset)) +
  geom_density(alpha = 0.3) +
  ggplot2::ylim(c(0,10))+
  geom_density(
    data = expit_treatment %>% 
      pivot_longer(
        cols = everything(),
        names_to = "Column",
        values_to = "Value"),
    aes(x = Value),
    color = "black",
    linetype = "dashed",
    inherit.aes = FALSE) +
  facet_wrap(~ Column) +
  theme_minimal() +
  labs(
    title = "Propensity score comparison",
    x = "Value",
    y = "Density")

ggplot2::ggsave(filename = paste0("inst/images/tryzone/",type,"/propensity_score_",type, "_", n,".pdf"), 
                width = 10, height = 8)

################################################################################
##########################       Cardinality      ##############################
################################################################################
quantiles_glm <- apply(r0_scores_glm[,1:ncol(SL.out$doptFactorPredict_test)], 2, 
                       function(x)stats::quantile(x, probs = 1 - alpha, na.rm = TRUE))

set_valued_policies_glm <- apply(data.frame(1:ncol(SL.out$doptFactorPredict_test)), 1,
                                 function(i){
                                   binary_confidence_set<- ifelse(margin_po_glm_new<quantiles_glm[i], 1, 0)
                                   idx <- which(binary_confidence_set  != 0, arr.ind = TRUE)
                                   split(idx[, "col"], factor(idx[, "row"],
                                                              levels = seq_len(nrow(binary_confidence_set))))})

mean_cardinality_glm <- lapply(set_valued_policies_glm, function(x)apply(do.call(rbind,lapply(x, length)), 2, mean))

quantiles_grf <- apply(r0_scores_glm[,1:ncol(SL.out$doptFactorPredict_test)], 2, 
                       function(x)stats::quantile(x,probs = 1-alpha, na.rm = TRUE))

set_valued_policies_grf <- apply(data.frame(1:ncol(SL.out$doptFactorPredict_test)), 1,
                                 function(i){
                                   binary_confidence_set<- ifelse(margin_po_grf_new<quantiles_grf[i], 1, 0)
                                   idx <- which(binary_confidence_set  != 0, arr.ind = TRUE)
                                   split(idx[, "col"], factor(idx[, "row"],
                                                              levels = seq_len(nrow(binary_confidence_set))))})
mean_cardinality_grf <- lapply(set_valued_policies_grf, function(x)apply(do.call(rbind,lapply(x, length)), 2, mean))
card_pv_grf <- cbind(do.call(rbind, mean_cardinality_grf), results_grf$oracular_spv[1:ncol(SL.out$doptFactorPredict_new)], 
                     results_grf$r_bar[1:ncol(SL.out$doptFactorPredict_new)])

ggplot()+
  geom_point(aes(x=card_pv_grf[,1], y = card_pv_grf[,2], colour = ifelse(card_pv_grf[,3]<0,1,0)))
################################################################################
############################       Heatmap      ################################
################################################################################
true_df <- as.data.frame(do.call(rbind, optimal_policy_test))%>%
  mutate(Row = row_number()) %>%
  pivot_longer(
    cols = -Row, 
    values_to = "Levels"
  ) %>%
  mutate(Levels = factor(Levels, levels=levels_A))

true_heatmap <- ggplot(true_df
                       %>% filter(`Row`<=10), 
                       aes(x = name, y = `Row`, fill = `Levels`)) +
  geom_tile(linewidth = 0.1)+
  scale_fill_viridis_d(option = "viridis", drop = FALSE) +
  scale_color_discrete()+
  scale_y_reverse() +
  xlab("Ground truth")+
  theme(
    axis.text.y = element_blank(),
    legend.position = "none", 
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(b = 44, r=-5)
  )

# 2. Combine and reshape
heatmap_data <- SL.out$doptFactorPredict_test %>% 
  as.data.frame() %>%
  mutate(Row = row_number()) %>%
  pivot_longer(
    cols = -Row,
    names_to = "Policy learning method",
    values_to = "Levels"
  ) %>%
  mutate(Levels = factor(Levels, levels = levels_A))

heatmap_pl <- ggplot(heatmap_data %>% 
                       filter(`Row`<=10), aes(x = `Policy learning method`, 
                         y = `Row`, fill = `Levels`)) +
  geom_tile(linewidth = 0.1)+
  scale_fill_manual(values = viridisLite::viridis(m), 
                    drop=FALSE, 
                    limits = levels_A) +
  scale_y_reverse() +
  ylab("")+
  theme(
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    plot.margin =  margin(l=-10)
  )

heatmap_both <- grid.arrange(true_heatmap, heatmap_pl, ncol=2, widths = c(1, 4))
ggplot2::ggsave(heatmap_both, filename = paste0("inst/images/tryzone/",type,"/heatmap_recommendations_",type, "_", n,".pdf"), 
                width = 10, height = 8)


