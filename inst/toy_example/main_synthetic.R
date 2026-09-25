# ── Set working directory  ──────────────────────────────────────────────────
root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning"
setwd(root.path)
if (!dir.exists(root.path)) {
  warning(sprintf("The directory '%s' does not exist. Creating it...", root.path))
  dir.create(root.path, recursive = TRUE)
}

subdirs <- c("images", "predictions")
for (subdir in subdirs) {
  subdir_path <- file.path(root.path, "inst", subdir)

  if (!dir.exists(subdir_path)) {
    message(sprintf("Creating subdirectory: %s", subdir))
    dir.create(subdir_path, recursive = TRUE)
  }
  if (subdir == "images") {
    for (img_folder in c(6000, 12000, 18000)) {
      img_path <- file.path(subdir_path, img_folder)
      if (!dir.exists(img_path)) {
        message(sprintf("Creating image subfolder: %s", img_folder))
        dir.create(img_path, recursive = TRUE)
      }
    }
  }
}

message("Directory check complete.")


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
seed <- 2026
set.seed(seed)
VFolds <- 3 # folds to split data
synthetic_scenario <- TRUE
type <- "normal" # additional name for images (here: type of synthetic scenario)

name_learner <- "new"


is_RCT <- ifelse(type=="normal", FALSE, TRUE)
RCT_file<- ifelse(is_RCT==TRUE,"RCT/", "non_RCT/")
n_samples <- c(6000, 12000, 18000)
ncov <- 4
n_test <- 100

random_rate <- seq(0,1,0.1) # random rates to test
n_rate <- length(random_rate) # number of random rates to test
alphas <- seq(0,1,0.05) # number of confidence levels to test

# ── Simulations for varying sample sizes  ─────────────────────────────────────
for (n in n_samples){
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
    #if(type=="normal"){c("X1","X2", "X3", "X4")}else{c("X1","X2")}# name for covariates in dataset
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

  family = ifelse(max(Y) <= 1 & min(Y) >= 0, "binomial", "gaussian") # in [0,1] or beyond
  SL.out$family = family
  ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new)))

  # ── 0) Divide data into three even sets ─────────────────────────────────────
  # ── Noisy label generation, scoring model & calibration ─────────────────────
  SL.out$folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y,
                                        cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds,
                                                                                          shuffle = TRUE))

  train1 <- SL.out$df_obs[SL.out$folds[[1]],] # generate noisy labels
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
  SL.out$g.reg.train_spv <- grf::probability_forest(X = X[SL.out$folds[[1]],], 
                                                    Y =  A[SL.out$folds[[1]]] %>% as.factor())
  
  #randomForest::randomForest(x = train1[,covariates_name],y = train1[,treatment_name])
  ## 1.3) Estimate A* (OTR) using experts
  # Training performed on train1
  # Two predictions:
  # (i) on test
  # (ii) on SL.out$df_new_sample
  X_train <- X[SL.out$folds[[1]],] 
  A_train <- A[SL.out$folds[[1]]]
  Y_train <- Y[SL.out$folds[[1]]]
  
  source("inst/toy_example/train_policies.R")

  # Generate noisy calibration labels
  unweighted_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_test,
                                             weights =rep(1/numalgs, numalgs),
                                             df_pred = test,
                                             levels = as.numeric(levels_A))
  unweighted_cal <- apply(apply(unweighted_probs, 1, function(x){
    rmultinom(1,1,prob=x)}), 2, which.max)

  # ── 2) Nonconformity score model (i.e. s(X,A)) ───────────────────────
  # Training performed on train2
  # Two predictions:
  # (i) on test
  # (ii) on SL.out$df_new_sample 
    SL.out$QAW.reg.train = grf::regression_forest(X = cbind(
      X[SL.out$folds[[2]],], A[SL.out$folds[[2]]]), 
      Y = Y[SL.out$folds[[2]]], seed = seed)
    
    potential_outcomes_test <- do.call(cbind,lapply(1:m, function(val) {
      new_data <- cbind(X[SL.out$folds[[3]],], factor(val, levels=levels_A) %>% as.numeric())
      stats::predict(SL.out$QAW.reg.train, newdata = new_data)$predictions}))
    
    potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
      new_data <- cbind(X_new, factor(val, levels=levels_A) %>% as.numeric())
      stats::predict(SL.out$QAW.reg.train, newdata = new_data)$predictions}))
  
  SL.out$g.reg.train <- grf::probability_forest(X = train2[,covariates_name],
                                                Y = train2[,treatment_name])

  # 2.2) Predict nonconformity scores (margin score)
  # Nonconformity scores on calibration data
  margin_po <-  margin_score(potential_outcomes_test)

  # Nonconformity scores on oracular data
  SL.out$true_score <- margin_po[cbind(1:nrow(test), SL.out$true_cal)]

  # Nonconformity scores on new data (used to generate sets)
  SL.out$new_scores <- margin_score(potential_outcomes_new)  # score for all potential outcomes from new data

  # Oracular nonconformity scores
  if(type=="normal"){
    margin_true_test <-margin_score(mu_P0_normal(test[,covariates_name]))
    SL.out$true_marginal_scores_new <- margin_score(mu_P0_normal(SL.out$df_new_sample[,covariates_name]))
  }else{
    margin_true_test <- margin_score(mu_P0_simplex_complicated(test[,covariates_name]))
    SL.out$true_marginal_scores_new <- margin_score(mu_P0_simplex_complicated(SL.out$df_new_sample[,covariates_name]))
  }
  SL.out$true_score_true <- margin_true_test[cbind(1:nrow(test), SL.out$true_cal)]

  # Randomness injection
  source("inst/randomness_injection.R")
  # Save results
  saveRDS(object = SL.out, file = paste0("inst/predictions/", type,"_",n,"_", name_learner,".rds"))
  
  # Create result table
  #source("inst/toy_example/table.R")
  # Evaluate set-valued policies
  source("inst/toy_example/metrics_synthetic.R")
}
# Generate plots
source("inst/toy_example/figures_synthetic.R")
