# ── Set working directory  ──────────────────────────────────────────────────
root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /CPolicySets"
setwd(root.path)

if (!dir.exists(root.path)) {
  warning(sprintf("The directory '%s' does not exist. Creating it...", root.path))
  dir.create(root.path, recursive = TRUE)
}

subdirs <- c("images", "predictions")

for (subdir in subdirs) {
  subdir_path <- file.path(root.path, subdir)
  
  if (!dir.exists(subdir_path)) {
    message(sprintf("Creating subdirectory: %s", subdir))
    dir.create(subdir_path, recursive = TRUE)
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
synthetic_scenario <- FALSE 
type <- #file name 

random_rate <- seq(0,0.6,0.1) # random rates to test 
n_rate <- length(random_rate) # number of random rates to test 
n_test <- 100 # number of repetitions for boxplots
alphas <- seq(0,0.5,0.05) # number of confidence levels to test 

SL.out<- list() # where results will be saved 

# ── Load data  ────────────────────────────────────────────────────────────────
all_data <- NULL #LOAD DATA 
  
SL.out$df_obs <- NULL # subsample of all data for training CP (ex. 3/4 of data)
n<- nrow(SL.out$df_obs)  # number of observations for training CP.
  
SL.out$df_new_sample <- NULL # subsample of all data for generating prediction sets (here 1/4 remaining)
  
SL.out$optimal_policy <- NULL
SL.out$optimal_policy_new <- NULL

# ── Define data parameters  ───────────────────────────────────────────────────
covariates_name <- # name for covariates in dataset
X <- SL.out$df_obs[,covariates_name] %>% as.matrix() 
X_new <- SL.out$df_new_sample[,covariates_name] %>% as.matrix()

treatment_name <-  # name of treatment indicator in dataset
A <- SL.out$df_obs[,treatment_name]
A_new <- SL.out$df_new_sample[,treatment_name]

levels_A <- levels(A) # treatment levels 
m <- length(levels(A)) # number of treatment levels 

outcome_name <-  # name of outcome in dataset
Y <- SL.out$df_obs[,outcome_name] 
Y_new <- SL.out$df_new_sample[,outcome_name] 

family = ifelse(max(Y) <= 1 & min(Y) >= 0, "binomial", "gaussian") # in [0,1] or beyond 
SL.out$family = family 
ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new))) 

# ── Define libraries for estimation procedures ────────────────────────────────
## Pseudo-label generation

# options : c("polle", "personalized", "mix")
expert_technique <- "polle" # indicator of estimation procedure 
#### Define the library of experts (Q)
if(expert_technique %in% c("polle", "mix")){
  q_learners <- list() %>%
    add_qlearner(name = "drql_lm", type = "drql", q_func = "q_glm", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "drql_rf", type = "drql", q_func = "q_rf", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "ql_lm", type = "ql",   q_func = "q_glm", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "ql_rf", type = "ql",   q_func = "q_rf", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "drql_xgb", type = "drql", q_func = "q_sl", 
                 sl_library = "SL.xgboost", 
                 action_name=treatment_name, covariates=covariates_name) %>% 
    add_qlearner(name = "ql_xgb", type = "ql",   q_func = "q_sl", 
                 sl_library = "SL.xgboost", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "ql_ranger", type = "ql",   q_func = "q_sl", 
                 sl_library = "SL.ranger", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "ptl_xgb", type = "ptl",  q_func = "q_sl", 
                 sl_library = "SL.xgboost", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "ptl_lm", type = "ptl",  q_func = "q_glm", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "ptl_rf", type = "ptl",  q_func = "q_rf", 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "policytree", type = "policytree", 
                 depth = 2, 
                 action_name=treatment_name, covariates=covariates_name) %>%
    add_qlearner(name = "hybrid_tree", type = "policytree", 
                 depth = 3, hybrid = TRUE, 
                 action_name=treatment_name, covariates=covariates_name)
  
  #### Define the propensity score learner (G)
  g_learner <- glearner(m, g_func = "g_rf", sl_library = NULL, num.trees = 500)  
}

## Nonconformity score learner (nuisance)
SL.library.nuisance <- c(
  "SL.mean", 
  "SL.glm", 
  "SL.xgboost", 
  "SL.ranger", 
  "SL.ksvm" 
)

# ── 0) Divide data into three even sets ───────────────────────────────────────
# ── Noisy label generation, scoring model & calibration ───────────────────────
SL.out$folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y,
                               cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds, 
                                                                                 shuffle = TRUE))

train1 <- SL.out$df_obs[SL.out$folds[[1]],] # generate noisy labels 
train2 <-  SL.out$df_obs[SL.out$folds[[2]],] # score model and nuisances
test <-  SL.out$df_obs[SL.out$folds[[3]],] # calibration 

# ── 1) Generate noisy labels (i.e. estimates of (A*,X)) ───────────────────────
## 1.1) Generate random labels (i.e. A_rd)
A_rd <- apply(data.frame(1:nrow(test)),1,
              function(i)sample(as.numeric(levels_A),size=1))

## 1.2) Estimate A* using experts 
if(expert_technique!= "personalized"){
  SL.init1 = expert_fit_predict(train1, test, new = SL.out$df_new_sample, 
                                covariates = covariates_name, 
                                treatment_name=treatment_name, 
                                outcome_name=outcome_name,
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
unweighted_cal <- apply(
  apply(unweighted_probs, 1, function(x){rmultinom(1, 1, prob=x)}), 
                        2, which.max)

# 2) Nonconformity score model (on train2)
# 2.1) Train nuisances 
## Outcome model
SL.out$QAW.reg.train = SuperLearner::SuperLearner( 
  Y=train2[,outcome_name], X = train2[,c(covariates_name,treatment_name)],
  SL.library=SL.library.nuisance, family = SL.out$family) 

## Propensity score
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

source("noise_injection.R")
saveRDS(object = SL.out, file = paste0("predictions/", type, ".rds"))
source("plots.R")
