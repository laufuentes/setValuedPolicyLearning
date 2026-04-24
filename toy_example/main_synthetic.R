# ── Set working directory  ──────────────────────────────────────────────────
root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /CPolicySets"
setwd(root.path)
if (!dir.exists(root.path)) {
  warning(sprintf("The directory '%s' does not exist. Creating it...", root.path))
  dir.create(root.path, recursive = TRUE)
}

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

# ── Load functions from src files  ────────────────────────────────────────────
source("src/synthetic_data.R")
source("src/non-conformity-scores.R")
source("src/label-estimation.R")
source("src/utils.R")
source("src/evaluation.R")

# ── General parameters  ──────────────────────────────────────────────────────────────
seed <- 2026
set.seed(seed)
VFolds <- 3 # folds to split data 
synthetic_scenario <- TRUE 
type <- "complex" # additional name for images (here: type of synthetic scenario)
is_RCT <- FALSE
RCT_file<- ifelse(is_RCT==TRUE,"RCT/", "non_RCT/")
n_samples <- c(6000, 12000, 18000)
ncov <- 4
n_test <- 1 

random_rate <- seq(0,1,0.1) # random rates to test 
n_rate <- length(random_rate) # number of random rates to test 
alphas <- seq(0,1,0.05) # number of confidence levels to test 

subdirs <- c("images", "predictions")
for (subdir in subdirs) {
  subdir_path <- file.path(root.path, subdir)
  
  if (!dir.exists(subdir_path)) {
    message(sprintf("Creating subdirectory: %s", subdir))
    dir.create(subdir_path, recursive = TRUE)
  }
  if (subdir == "images") {
    for (img_folder in n_samples) {
      img_path <- file.path(subdir_path, img_folder)
      if (!dir.exists(img_path)) {
        message(sprintf("Creating image subfolder: %s", img_folder))
        dir.create(img_path, recursive = TRUE)
      }
    }
  }
}

message("Directory check complete.")


for (n in n_samples){
  SL.out<- list() # list where results will be saved 
  
  # Synthetic data generation 
  ## Training observations 
  exp <- generate_data(n, ncov = ncov, type=type, is_RCT=is_RCT, seed = seed)  
  SL.out$df_obs <- exp[[1]] # extract observational data 
  df_complete <- exp[[2]] # extract complete data (unavailable in real scenarios)
  SL.out$optimal_policy <- exp[[3]] # extract optimal policy (unavailable in real scenarios)
  SL.out$potential_outcomes <- df_complete %>% select(starts_with("Potential_outcomes."))
  
  ### Test observations 
  exp_new_sample <- generate_data(n/2, ncov=ncov, type=type) # generate test observations 
  SL.out$df_new_sample <- exp_new_sample[[1]]  # extract observational data 
  SL.out$optimal_policy_new <- exp_new_sample[[3]] # extract optimal policy (unavailable in real scenarios)
  SL.out$potential_outcomes <- exp_new_sample[[2]] %>% select(starts_with("Potential_outcomes."))
  
  # ── Define data parameters  ───────────────────────────────────────────────────
  covariates_name <- c("X1","X2", "X3", "X4") # name for covariates in dataset
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

  # ── 0) Divide data into three even sets ───────────────────────────────────────
  # ── Noisy label generation, scoring model & calibration ───────────────────────
  SL.out$folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y,
                                        cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds, 
                                                                                          shuffle = TRUE))
  
  train1 <- SL.out$df_obs[SL.out$folds[[1]],] # generate noisy labels 
  train2 <-  SL.out$df_obs[SL.out$folds[[2]],] # score model and nuisances
  test <-  SL.out$df_obs[SL.out$folds[[3]],] # calibration 
  optimal_policy_test <- SL.out$optimal_policy[SL.out$folds[[3]]]
  true_potential_outcomes_test <- df_complete[SL.out$folds[[3]],] %>% select(starts_with("Potential_outcomes."))
  
  # ── 1) Generate noisy labels (i.e. estimates of (A*,X)) ───────────────────────
  ## 1.0) True calibration samples (only for synthetic data)
  
  ## Pseudo-label generation
  # options : c("polle", "personalized", "mix")
  expert_technique <- "polle" # indicator of estimation procedure 
  #### Define the library of experts (Q)
  if(expert_technique %in% c("polle", "mix")){
    q_learners <- list() %>%
      add_qlearner(name = "drql_lm", type = "drql",   q_func = "q_glm", 
                   action_name=treatment_name, covariates=covariates_name)
    
    #### Define the propensity score learner (G)
    g_learner <- glearner(m, g_func = "g_rf", sl_library = NULL, num.trees = 500)  
  }
  
  SL.out$true_cal<- apply(data.frame(1:nrow(test)),1,function(x){
    el <- optimal_policy_test[[x]]
    length_el <- length(el)
    if(length_el==1){el}else{el[sample(length_el,1)]}}) 
  
  ## 1.1) Generate random labels (i.e. A_rd)
  A_rd <- apply(data.frame(1:nrow(test)),1,function(i)sample(as.numeric(levels_A),size=1))
  
  ## 1.2) Estimate A* using experts 
  if(expert_technique!= "personalized"){
    SL.init1 = expert_fit_predict(train1, test, new = SL.out$df_new_sample, 
                                  covariates = covariates_name, 
                                  treatment_name=treatment_name, outcome_name=outcome_name,
                                  qlearners_list = q_learners, g_model=g_learner)
    
    SL.out$libraryNames <- SL.init1[["learner_names"]] # expert names 
    numalgs <- length(SL.out$libraryNames) # number of experts
    
    # predictions on calibration
    SL.out$doptFactorPredict_test <- SL.init1[["expert_policies"]]  
    # predictions on new data
    SL.out$doptFactorPredict_new <- SL.init1[["expert_policies_new"]]  
    
    if(expert_technique=="mix"){
      additional_methods <- NULL # replace with additional learning methods 
      SL.out$libraryNames <- c(SL.out$libraryNames, additional_methods)
      other_preds_test <- NULL # replace: predictions on calibration set 
      SL.out$doptFactorPredict_test <- cbind(SL.out$doptFactorPredict_test, 
                                             other_preds_test)
      other_preds_new <- # replace: predictions on new data
        SL.out$doptFactorPredict_new <- cbind(SL.out$doptFactorPredict_test, 
                                              other_preds_new)
    }
  } else{
    # Define experts with other learners 
    SL.out$libraryNames <- NULL # replace with learner names 
    numalgs <- length(SL.out$libraryNames)
    
    SL.out$doptFactorPredict_test <- NULL # replace with predictions on calibration 
    SL.out$doptFactorPredict_new <- NULL # replace with predictions on new
  }
  
  #### 1.2.a) Unweighted calibration samples
  unweighted_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_test,
                                             weights =rep(1/numalgs, numalgs), 
                                             df_pred = test, 
                                             levels = as.numeric(levels_A))
  unweighted_cal <- apply(apply(unweighted_probs, 1, function(x){rmultinom(1,1,prob=x)}),2,which.max)
  
  # 2) Nonconformity score model (on train2)
  # 2.1) Train nuisances 
  SL.library.nuisance <- c(
    "SL.mean", 
    "SL.glm", 
    "SL.xgboost", 
    "SL.ksvm",
    "SL.ranger" 
  )
  SL.out$QAW.reg.train = SuperLearner::SuperLearner( # Outcome model
    Y=train2[,outcome_name], X = train2[,c(covariates_name,treatment_name)],
    SL.library=SL.library.nuisance, family = SL.out$family) 
  
  SL.out$g.reg.train <- randomForest::randomForest(x = train2[,covariates_name], 
                                                   y = train2[,treatment_name])
  
  # 2.2) Predict nonconformity scores (margin score)
  #### Calibration data 
  potential_outcomes_test <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- test[, c(covariates_name, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(SL.out$QAW.reg.train, newdata = new_data)$pred}))
  
  margin_po <-  margin_score(potential_outcomes_test) # score for all potential outcomes from test
  #### New data 
  potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- SL.out$df_new[, c(covariates_name, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(SL.out$QAW.reg.train, newdata = new_data)$pred}))
  
  SL.out$new_scores <- margin_score(potential_outcomes_new)  # score for all potential outcomes from new data
  
  SL.out$true_score <- margin_po[cbind(1:nrow(test), SL.out$true_cal)] 
  
  if(type=="normal"){
    margin_true_test <-margin_score(mu_P0_normal(test[,covariates_name]))
    SL.out$true_marginal_scores_new <- margin_score(mu_P0_normal(SL.out$df_new_sample[,covariates_name]))
  }else{
    margin_true_test <- margin_score(mu_P0_simplex_complicated(test[,covariates_name]))
    SL.out$true_marginal_scores_new <- margin_score(mu_P0_simplex_complicated(SL.out$df_new_sample[,covariates_name]))
  }
  SL.out$true_score_true <- margin_true_test[cbind(1:nrow(test), SL.out$true_cal)]

  source("noise_injection.R")
  saveRDS(object = SL.out, file = paste0("predictions/", type,"_",n,".rds"))
  source("toy_example/table.R")
  source("toy_example/metrics_synthetic.R")
}
source("toy_example/figures_synthetic.R")
