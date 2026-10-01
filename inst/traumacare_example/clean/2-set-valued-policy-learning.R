root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning/"
setwd(root.path)
seed <- 2026
set.seed(seed)

names <- readRDS("inst/traumacare_example/intermediate/preprocessing.rds")
source("inst/libraries.R")
source("R/utils.R")
source("R/evaluation.R")

# ── Set environment for (conformal) policy learning ───────────────────────────
SL.out<- list() # where results will be saved
type <- "traumacare" #file name
random_rate <- c(0, 0.1, 0.25, 0.5)
n_rate <- length(random_rate)
alpha <- 0.1

# ── Load data  ────────────────────────────────────────────────────────────────
train_imp_onehot <- read.csv("inst/traumacare_example/intermediate/train_imp_one_hot.csv")
SL.out$df_obs <- train_imp_onehot # use train_imp (if want to use no-one-hot encoded version)
n<- nrow(SL.out$df_obs)  # number of observations for training CP.

test_imp_onehot <- read.csv("inst/traumacare_example/intermediate/test_imp_one_hot.csv")
SL.out$df_new_sample <- test_imp_onehot # use test_imp (if want to use no-one-hot encoded version)

# ── Define data parameters  ───────────────────────────────────────────────────
# baseline covariates
covariates_name <- names$covariates_name # use covariate_name (if want to use no-one-hot encoded version ) 
X <- SL.out$df_obs[, covariates_name] %>% 
  as.matrix() %>% 
  apply(2, as.numeric) # Covariates for training data 

X_new <- SL.out$df_new_sample[,covariates_name] %>% 
  as.matrix() %>% 
  apply(2, as.numeric) # Covariates for test data 

# treatment
treatment_name <- names$treatment_name
A <- SL.out$df_obs[,treatment_name] # Treatment vector for training data 

levels_A <- levels(A %>% as.factor()) # treatment levels ("0" and "1")
m <- length(levels_A) # number of treatment levels (2)

# outcome 
outcome_name <- names$outcome_name
Y <- SL.out$df_obs[,outcome_name] 
ab <- c(min(c(Y)),max(c(Y))) 

# ── Split train into three  ───────────────────────────────────────────────────
shuffled_indices <- sample(1:n)
cut1 <- round(0.5 * n)
cut2 <- round(0.2 * n)
custom_folds <- list(
  Fold1 = shuffled_indices[1:cut1],
  Fold2 = shuffled_indices[(cut1 + 1): (cut1 + cut2)],
  Fold3 = shuffled_indices[(cut1 + cut2 + 1):n])

SL.out$folds <- SuperLearner::CVFolds(n, id = NULL, Y = Y,
                                      cvControl = SuperLearner::SuperLearner.CV.control(V = 3L,
                                                                                        validRows = custom_folds))
# train the set-valued policy learning methods
train <- SL.out$df_obs[SL.out$folds[[1]],] 
n_train <- nrow(train)
# alpha selection
pseudo.test.predict <-  SL.out$df_obs[SL.out$folds[[2]],]  # pseudo.test
#SL.out$folds_pseudo <- SuperLearner::CVFolds(nrow(pseudo.test), id = NULL,Y = Y[SL.out$folds[[2]]],
#                                      cvControl = SuperLearner::SuperLearner.CV.control(V = 2L))
#pseudo.test.predict <- pseudo.test[SL.out$folds_pseudo[[1]],]    
#pseudo.test.nuisance <- pseudo.test[SL.out$folds_pseudo[[2]],] 
# alpha evaluation
evaluation <-  SL.out$df_obs[SL.out$folds[[3]],] 

# Set-valued policy learning ───────────────────────────────────────────────────
####### Conformal policy learning ##############################################      
# ── 0) Divide data into three even sets ───────────────────────────────────────
# ── Noisy label generation, scoring model & calibration ───────────────────────
set.seed(seed)
shuffled_indices <- sample(1:n_train)
cut1 <- cut2 <- round(0.4 * n_train)
custom_folds_conformal <- list(
  Fold1 = shuffled_indices[1:cut1],
  Fold2 = shuffled_indices[(cut1 + 1): (cut1 + cut2)],
  Fold3 = shuffled_indices[(cut1 + cut2 + 1):n_train])

set.seed(seed)
SL.out$folds_conformal <- SuperLearner::CVFolds(n_train, id = NULL,Y = Y[SL.out$folds[[1]]],
                                                cvControl = SuperLearner::SuperLearner.CV.control(V = 3L,
                                                                                                  validRows = custom_folds_conformal))

train1 <- SL.out$df_obs[SL.out$folds_conformal[[1]],] # generate noisy labels
train2 <-  SL.out$df_obs[SL.out$folds_conformal[[2]],] # score model and nuisances
calibration <-  SL.out$df_obs[SL.out$folds_conformal[[3]],] # calibration

# ── 1) Black-box label generation (i.e. estimates of (X,A*)) ──────────────────
# (also trains the nuisance for GLB) 
source("inst/traumacare_example/clean/2.b-Train_policies.R")

# Generate noisy calibration labels
# using an aggregation of policies 
unweighted_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_test,
                                           weights =rep(1/numalgs, numalgs),
                                           df_pred = calibration,
                                           levels = as.numeric(levels_A)) # probability distribution

colnames(unweighted_probs) <- levels_A

# # aggregation of multiple policies 
SL.out$unweighted_cal <- apply(
  apply(unweighted_probs, 1, 
        function(x){rmultinom(1, 1, prob=x)}), 2, which.max)-1  # generate calibration samples 

# # single policy 
SL.out$single_policy_cal <- SL.out$doptFactorPredict_test[,c("proba.forest","ql.SL")] 

noisy_labels <- cbind(unweighted = SL.out$unweighted_cal, SL.out$single_policy_cal)

A_rd <- apply(data.frame(1:nrow(calibration)),1,function(i)sample(as.numeric(levels_A),size=1))

perturbed_noisy_labels <- apply(noisy_labels, 2, function(policy) {
  sapply(random_rate, function(rate) {
    mix_factor <- rbinom(nrow(calibration), 1, prob = rate)
    mix_factor*A_rd  + (1 - mix_factor) * policy
  }) |> as.data.frame() |> setNames(random_rate)})

SL.out$perturbed_noisy_labels <- perturbed_noisy_labels


# ── 2) Nonconformity score model (i.e. s(X,A))  ───────────────────────────────
# Training performed on train2
# Two predictions:
# (i) on calibration (calibration)
# (ii) on test set (SL.out$df_new_sample)

# 2.1) Train nuisance  ─────────────────────────────────────────────────────────
## Outcome model (Q-model)  
QAW.reg.train_conformal = grf::probability_forest(
  X = cbind(X[SL.out$folds_conformal[[2]],],A[SL.out$folds_conformal[[2]]]), 
  Y = Y[SL.out$folds_conformal[[2]]] %>% as.factor())

# 2.2) Predict nonconformity scores (margin score)  ────────────────────────────
# Nonconformity scores on calibration data
potential_outcomes_cal <- do.call(cbind,lapply(0:1, function(val) {
  new_data <- cbind(X[SL.out$folds_conformal[[3]],], val)
  stats::predict(QAW.reg.train_conformal, newdata = new_data)$predictions[,2]}))

margin_po <-  margin_score(potential_outcomes_cal) # score for all potential outcomes from calibration

r1_score <- data.frame(random = margin_po[cbind(1:nrow(calibration), A_rd+1)])
ecdf_data <- apply(cbind(SL.out$doptFactorPredict_test, 
                         unweighted=SL.out$unweighted_cal), 2, 
                   function(x)margin_po[cbind(1:nrow(calibration),x+1)])|> 
  as.data.frame() |>
  pivot_longer(cols=everything(), names_to = "methods", values_to = "value")

ggplot(ecdf_data, aes(x = value, colour = methods)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  geom_hline(yintercept = 1-alpha, colour = "red") +
  stat_ecdf(
    data = r1_score,
    aes(x = random, colour = "Random labels"),
    linetype = "dashed", colour="gray",
    linewidth = 1.2
  ) +
  labs(y = "ECDF", x = "Value", colour = "Method")
ggplot2::ggsave(filename = paste0("inst/traumacare_example/images/ecdf_proba_forest.pdf"), 
                width = 10, height = 8)

# Nonconformity scores on new data (used to generate sets)
potential_outcomes_new <- do.call(cbind,lapply(0:1, function(val) {
  new_data <- cbind(X_new, val)
  stats::predict(QAW.reg.train_conformal, newdata = new_data)$predictions[,2]}))

potential_outcomes_pseudo <- do.call(cbind,lapply(0:1, function(val) {
  new_data <- cbind(X[SL.out$folds[[2]],], val) # folds_pseudo[[1]]
  stats::predict(QAW.reg.train_conformal, newdata = new_data)$predictions[,2]}))


# # ── 3) Calibration: generate the nonconformity score density  ─────────────────
# densities based on noisy labels 
# noisy labels: aggregation of multiple policies 

margin_po_perturbed <- lapply(perturbed_noisy_labels, function(df) {
  apply(df, 2, function(x){
    margin_po[cbind(1:nrow(calibration), x+1)]})})

SL.out$calibration_scores <- margin_po_perturbed

# densities used for prediction  
SL.out$new_scores <- margin_score(potential_outcomes_new)  # final test
SL.out$pseudo_scores <- margin_score(potential_outcomes_pseudo)  # pseudo


saveRDS(SL.out, file = "inst/traumacare_example/intermediate/SL.out.rds")

# 3) Train the nuisances for alpha evaluation ────────────────────────────────────────
# Final evaluation 
## Outcome model (Q-model)  
QAW.reg.train.alpha.grf = grf::probability_forest(
  X = cbind(evaluation[,covariates_name],
            evaluation[,treatment_name]), 
  Y = evaluation[,outcome_name] %>% as.factor())

Q.all.pseudo.grf <-  do.call(cbind,lapply(0:1, function(val) {
  new_data <- cbind(pseudo.test.predict[,covariates_name], val)
  stats::predict(QAW.reg.train.alpha.grf, newdata = new_data)$predictions[,2]}))


# QAW.reg.train.alpha.glm = stats::glm(formula = formula(paste(outcome_name, "~", paste(c(covariates_name,treatment_name), collapse = "+"))),
#                                      data = evaluation, family = "binomial")
# 
# Q.all.pseudo.glm <-  do.call(cbind,lapply(0:1, function(val) {
#   new_data <- pseudo.test.predict[,c(covariates_name, treatment_name)]
#   new_data[,treatment_name] <- val
#   stats::predict(QAW.reg.train.alpha.glm, newdata = new_data, type="response")}))


## Propensity score (G-model)
g.reg.train.alpha.grf <- grf::probability_forest(X = evaluation[,covariates_name],
                                                 Y = evaluation[,treatment_name]%>% 
                                                   as.factor())

gAW.pred.pseudo.grf <- stats::predict(g.reg.train.alpha.grf, newdata = pseudo.test.predict[, covariates_name])$predictions
hist(gAW.pred.pseudo.grf[,2])

# g.reg.train.alpha.glm = stats::glm(formula = formula(paste(treatment_name, "~", paste(covariates_name, collapse = "+"))),
#                                    data = evaluation, family = "binomial")
# 
# ps<- stats::predict(g.reg.train.alpha.glm, newdata = pseudo.test.predict[, covariates_name], type="response") %>% as.numeric()
# gAW.pred.pseudo.glm <- cbind(1-ps, ps)
# hist(ps)

source("inst/traumacare_example/clean/3-Evaluation_with_r.R")
