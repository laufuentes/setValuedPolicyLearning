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
synthetic_scenario <- FALSE
type <- #file name

random_rate <- seq(0,0.6,0.1) # random rates to test
n_rate <- length(random_rate) # number of random rates to test
n_test <- 1 # number of repetitions for boxplots
alphas <- seq(0,0.5,0.05) # number of confidence levels to test

SL.out<- list() # where results will be saved

# ── Load data  ────────────────────────────────────────────────────────────────
all_data <- NULL #LOAD DATA

SL.out$df_obs <- NULL # subsample of all data for training CP (ex. 3/4 of data)
n<- nrow(SL.out$df_obs)  # number of observations for training CP.

SL.out$df_new_sample <- NULL # subsample of all data for generating prediction sets (here 1/4 remaining)

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

second_outcome_name <- # name of adverse event outcome
Xi <- SL.out$df_obs[,second_outcome_name]
Xi_new <- SL.out$df_new_sample[,second_outcome_name]
ab_xi <- c(min(c(Xi,Xi_new)),max(c(Xi,Xi_new)))

# ── Define libraries for estimation procedures ────────────────────────────────
## Noisy label generation
### options : c("polle", "personalized", "mix")
expert_technique <- "personalized" # indicator of estimation procedure

## Nonconformity score learner (nuisance)
SL.library.nuisance <- c(
  "SL.mean",
  "SL.glm",
  "SL.xgboost",
  "SL.ranger",
  "SL.ksvm")

# ── 0) Divide data into three even sets ───────────────────────────────────────
# ── Noisy label generation, scoring model & calibration ───────────────────────
SL.out$folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y,
                                      cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds,
                                                                                        shuffle = TRUE))

train1 <- SL.out$df_obs[SL.out$folds[[1]],] # generate noisy labels
train2 <-  SL.out$df_obs[SL.out$folds[[2]],] # score model and nuisances
test <-  SL.out$df_obs[SL.out$folds[[3]],] # calibration

# ── 1) Black-box label generation (i.e. estimates of (X,A*)) ───────────────────────
## 1.1) Generate random labels (i.e. A_rd)
A_rd <- apply(data.frame(1:nrow(test)),1,
              function(i)sample(as.numeric(levels_A),size=1))

## 1.2) Estimate A* (OTR) using experts
# Training performed on train1
# Two predictions:
# (i) on test
# (ii) on SL.out$df_new_sample

# Define experts with other learners
SL.out$libraryNames <- NULL # macf
numalgs <- length(SL.out$libraryNames)

# Note: predictions should keep numeric levels consistent
# with SL.out$df_obs[,treatment_name] (i.e. numeric)
#source("external_estimation_file.R")
pred_calibration <- NULL # replace with predictions on calibration
pred_new_data <- NULL # replace with predictions on new
SL.out$doptFactorPredict_test <- array(as.numeric(pred_calibration),
                                       dim = c(length(pred_calibration), numalgs),
                                       dimnames = list(NULL, SL.out$libraryNames))
SL.out$doptFactorPredict_new <- array(as.numeric(pred_new_data),
                                      dim = c(length(pred_new_data), numalgs),
                                      dimnames = list(NULL, SL.out$libraryNames))
# Generate noisy calibration labels
unweighted_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_test,
                                           weights =rep(1/numalgs, numalgs),
                                           df_pred = test,
                                           levels = as.numeric(levels_A))
unweighted_cal <- apply(
  apply(unweighted_probs, 1, function(x){rmultinom(1, 1, prob=x)}), 2, which.max)

# ── 2) Nonconformity score model (i.e. s(X,A)) ───────────────────────
# Training performed on train2
# Two predictions:
# (i) on test
# (ii) on SL.out$df_new_sample

# 2.1) Train nuisances
# Primary outcome (Y) model
SL.out$QAW.reg.train = SuperLearner::SuperLearner(
  Y=train2[,outcome_name], X = train2[,c(covariates_name,treatment_name)],
  SL.library=SL.library.nuisance, family = SL.out$family)

# Second outcome (Xi) model
SL.out$QxiAW.reg.train = SuperLearner::SuperLearner(
  Y=train2[,second_outcome_name], X = train2[,c(covariates_name,treatment_name)],
  SL.library=SL.library.nuisance, family = SL.out$family)

## Propensity score  (G-model)
SL.out$g.reg.train <- randomForest::randomForest(x = train2[,covariates_name],
                                                 y = train2[,treatment_name])

# 2.2) Compute nonconformity scores (margin score)
# Nonconformity scores on calibration data
potential_outcomes_test <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- test[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  SuperLearner::predict.SuperLearner(SL.out$QAW.reg.train, newdata = new_data)$pred}))

margin_po <-  margin_score(potential_outcomes_test) # score for all potential outcomes from test

# Nonconformity scores on new data (used to generate sets)
potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- SL.out$df_new[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  SuperLearner::predict.SuperLearner(SL.out$QAW.reg.train, newdata = new_data)$pred}))

SL.out$new_scores <- margin_score(potential_outcomes_new)  # score for all potential outcomes from new data

# Randomness injection
source("inst/randomness_injection.R")
# Save results
saveRDS(object = SL.out, file = paste0("inst/predictions/", type, ".rds"))

# Evaluate set-valued policies and generate plot
source("inst/ivf_example/plots_ivf.R")
