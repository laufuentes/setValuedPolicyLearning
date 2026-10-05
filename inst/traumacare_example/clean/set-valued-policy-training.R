root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning/"
setwd(root.path)
seed <- 2026
set.seed(seed)

names <- readRDS("inst/traumacare_example/intermediate/preprocessing.rds")


# ── Load functions from R folder  ────────────────────────────────────────────
source("inst/libraries.R")
source("R/utils.R")
source("R/evaluation.R")
source("inst/traumacare_example/clean/train_policies.R")

# ── General parameters  ───────────────────────────────────────────────────────
random_rate <- c(0, 0.1, 0.25, 0.5)
n_rate <- length(random_rate)
alpha <- 0.1
z <- qnorm(1 - alpha/2)

# ── Load data  ──────────────────────────────────────────────
### Train samples
train_imp_onehot <- read.csv("inst/traumacare_example/intermediate/train_imp_one_hot.csv")
df_obs <- train_imp_onehot # use train_imp (if want to use no-one-hot encoded version)
n<- nrow(df_obs)  # number of observations. 

### Test samples (for final predictions)
test_imp_onehot <- read.csv("inst/traumacare_example/intermediate/test_imp_one_hot.csv")
df_new_sample <- test_imp_onehot # use test_imp (if want to use no-one-hot encoded version)

# ── Define data parameters  ─────────────────────────────────────────────────
# Baseline covariates
covariates_name <- names$covariates_name # use covariate_name (if want to use no-one-hot encoded version ) 
X <- df_obs[, covariates_name] %>% 
  as.matrix() %>% 
  apply(2, as.numeric) # Covariates for training data 

X_new <- df_new_sample[,covariates_name] %>% 
  as.matrix() %>% 
  apply(2, as.numeric) # Covariates for test data 

# Treatment
treatment_name <- names$treatment_name
A <- df_obs[,treatment_name] %>% as.factor() # Treatment vector for training data 
levels_A <- levels(A) # treatment levels
m <- length(levels(A)) # number of treatment levels

# Outcome
outcome_name <- names$outcome_name
Y <- df_obs[,outcome_name] 
ab <- c(min(Y),max(Y))

# ── Split data into 3 folds  ────────────────────────────────────────────────
# Fold 1: Set-valued policy learning
# Fold 2: Predict set-valued policies (r selection)
# Fold 3: Train nuisances to evaluate on Fold2
set.seed(seed)
shuffled_indices <- sample(1:nrow(df_obs))
cut1 <- round(0.5 * nrow(df_obs))
cut2 <- round(0.2 * nrow(df_obs))
custom_folds <- list(
  Fold1 = shuffled_indices[1:cut1],
  Fold2 = shuffled_indices[(cut1 + 1): (cut1 + cut2)],
  Fold3 = shuffled_indices[(cut1 + cut2 + 1):nrow(df_obs)])
  
folds <- SuperLearner::CVFolds(nrow(df_obs), id = NULL, Y = df_obs[,outcome_name] ,
                                        cvControl = SuperLearner::SuperLearner.CV.control(V = 3L,
                                                                                          validRows = custom_folds))
  
train <- df_obs[folds[[1]],] # train set-valued policies 
n_train <- nrow(train)
pseudo.test.predict <-  df_obs[folds[[2]],] # predict set-valued policies 
evaluation <-  df_obs[folds[[3]],] # train nuisances for set-valued policy evaluation
  
# Train nuisances for evaluation
QAW.reg.train.r = grf::probability_forest(
X = cbind(evaluation[,covariates_name],
              evaluation[,treatment_name]), 
Y = evaluation[,outcome_name] %>% as.factor())
  
Q.all.pseudo.r <-  do.call(cbind,lapply(0:1, function(val) {
    new_data <- cbind(pseudo.test.predict[,covariates_name], val)
    stats::predict(QAW.reg.train.r, newdata = new_data)$predictions[,2]}))
  
g.reg.train.r <- grf::probability_forest(X = evaluation[,covariates_name],
                                           Y = evaluation[,treatment_name]%>% 
                                             as.factor())

gAW.pred.pseudo.r <- stats::predict(g.reg.train.r, 
                                        newdata = pseudo.test.predict[, covariates_name])$predictions
  
# ── SET_VALUED POLICY LEARNING ────────────────────────────────────────────────

# ── 1. CONFORMAL SET_VALUED POLICY LEARNING ───────────────────────────────────
# ── 0) Divide data into three even sets ───────────────────────────────────────
# ── Noisy label generation, scoring model & calibration ───────────────────────
set.seed(seed)
shuffled_indices <- sample(1:n_train)
cut1 <- cut2 <- round(0.4 * n_train)
custom_folds_conformal <- list(
    Fold1 = shuffled_indices[1:cut1],
    Fold2 = shuffled_indices[(cut1 + 1): (cut1 + cut2)],
    Fold3 = shuffled_indices[(cut1 + cut2 + 1):n_train])
  
folds_conformal <- SuperLearner::CVFolds(n_train, id = NULL, 
                                           Y = train[,outcome_name],
                                           cvControl = SuperLearner::SuperLearner.CV.control(V = 3L, validRows = custom_folds_conformal))
  
train1 <- train[folds_conformal[[1]],] # generate noisy labels
train2 <-  train[folds_conformal[[2]],] # score model and nuisances
calibration <-  train[folds_conformal[[3]],] # calibration
  
# ── 1) Black-box label generation (i.e. estimates of (X,A*)) ───────────────────────
## 1.1) Generate random labels (i.e. A_rd)
A_rd <- sample(as.numeric(levels_A), size = nrow(calibration), replace = TRUE)
  
## 1.2) Estimate A* (OTR) using experts
# Training performed on train1
trained_policies_results <- train_policies(train_b = df_obs, train1 = train1, 
                 calibration = calibration, pseudo_test = pseudo.test.predict,
                 seed = seed) # GLB trained inside
  
saveRDS(trained_policies_results,
          file="inst/traumacare_example/images/trained_policies_results.rds")

# Extract results from training
doptFactorPredict_cal <- trained_policies_results$doptFactorPredict_cal 
doptFactorPredict_pseudo <- trained_policies_results$doptFactorPredict_pseudo
doptFactorPredict_new_naive <- trained_policies_results$doptFactorPredict_new_naive 
doptFactorPredict_pseudo_naive <- trained_policies_results$doptFactorPredict_pseudo_naive
model.glb.pf <- trained_policies_results$model.glb.pf 
model.glb.glm <- trained_policies_results$model.glb.glm
selected_methods <- trained_policies_results$selected_methods 
unweighted.naive_new <- trained_policies_results$unweighted.naive
unweighted.pseudo.naive <- trained_policies_results$unweighted.pseudo.naive
numalgs <- ncol(doptFactorPredict_cal) # number of baselines
  
# Unweighted aggregation of baselines
unweighted_probs <- weighted_probs_experts(fitted_experts = doptFactorPredict_cal,
                                             weights =rep(1/numalgs, numalgs),
                                             df_pred = calibration,
                                             levels = as.numeric(levels_A))
  
unweighted_aggregation <- apply(apply(unweighted_probs, 1, function(x){
    rmultinom(1,1,prob=x)}), 2, which.max)
  
# Extract baseline policies for conformal procedure
policy_cal <- doptFactorPredict_cal[,selected_methods]
  
# Randomness injection
policy_cal_r <- matrix(0, nrow = nrow(calibration), ncol = n_rate)
for (i in seq_len(n_rate)){
    r <- random_rate[i]
    mix_factor<- stats::rbinom(nrow(calibration), 1, prob=r) # R ~ Ber(r)
    policy_cal_r[,i] <- mix_factor * A_rd + (1 - mix_factor) * policy_cal
}
  
selected_methods <- c(selected_methods,
                        lapply(random_rate, function(r) paste0(selected_methods, " (r=", r, ")")))
# ── 2) Train nonconformity score model (i.e. s(X,A)) ───────────────────────
# Training performed on train2
# Two predictions:
# (i) on calibration
# (ii) on df_new_sample
QAW.reg.train_conformal = grf::probability_forest(
    X = cbind(train2[,c(covariates_name,treatment_name)]), 
    Y = train2[,outcome_name] %>% as.factor())
  
# ── 3) Calibration step  ────────────────────────────────────────────────────
potential_outcomes_cal <- do.call(cbind,lapply(0:1, function(val) {
    new_data <- calibration[,c(covariates_name,treatment_name)]
    new_data[,treatment_name] <- val 
    stats::predict(QAW.reg.train_conformal, newdata = new_data)$predictions[,2]}))
  
# Compute margin score on calibration data 
margin_po <-  margin_score(potential_outcomes_cal)
  
# Extract scores for different label types
r0_scores_policy <- apply(policy_cal_r,2,function(x){
  margin_po[cbind(1:nrow(calibration), x)]})  # Baseline policies and perturbed labels 
r0_scores_aggregation <- margin_po[cbind(1:nrow(calibration), 
                                           unweighted_aggregation)] # aggregation of baselines
r1_score <- margin_po[cbind(seq_len(nrow(calibration)), A_rd)] # random labels 
r_true <- margin_po[cbind(seq_len(nrow(calibration)), true_cal)] # true labels
   
potential_outcomes_pseudo <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- pseudo.test.predict[, c(covariates_name, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(QAW.reg.train, newdata = new_data)$pred}))
  
# Compute margin score on pseudo test data (for r selection)
margin_po_pseudo <-  margin_score(potential_outcomes_pseudo)
  
potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- df_new_sample[, c(covariates_name, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(QAW.reg.train, newdata = new_data)$pred}))
  
# Compute margin score on new data (for final prediction)
margin_po_new <-  margin_score(potential_outcomes_new)
  
# Build and evaluate conformal set-valued policies   ─────────────────────────
quantiles_and_svp <- list()
# Baseline policy-based noisy labels
conf_policy <- apply(r0_scores_policy, 2, function(x){
  quant <- stats::quantile(x, (1-alpha))
  binary_to_confidence_set(margin_po_pseudo < quant)
}) # evaluation
  
results.policy <- lapply(conf_policy, function(x){
table.evaluation.real(x, prop_score_new = gAW.pred.pseudo.r, 
                          potential_outcomes = Q.all.pseudo.r,
                          df_new_sample = df_new_sample,
                          levels_A = levels_A, 
                          treatment_name = treatment_name, 
                          outcome_name = outcome_name)})
  
# Aggregation-based noisy labels
quant_agg <- quantile(r0_scores_aggregation, 1 - alpha)
conf_set_agg <- binary_to_confidence_set(margin_po_new < quant_agg)
list.set.valued.policies[["conf.agg"]] <- conf_set_agg
results_agg <- table.evaluation.real(conf_set_agg, 
                                       prop_score_new = gAW.pred.pseudo.r, 
                                       potential_outcomes = Q.all.pseudo.r,
                                       df_new_sample = df_new_sample,
                                       levels_A = levels_A, 
                                       treatment_name = treatment_name, 
                                       outcome_name = outcome_name)
  
# ── 2.GREATEST LOWER BOUND (GLB) ──────────────────────────────────────────────
# ── Using regression forest for estimation ──────────────────────────────────
lowers <- uppers <- matrix(0, nrow=nrow(df_new_sample), ncol=m)
for (l in as.numeric(levels_A)){
    data_l <- data.frame(df_new_sample[,covariates_name], Treatment=l)
    pred <- stats::predict(glb.model.grf, newdata = data_l, estimate.variance = TRUE)
    se <- sqrt(pred$variance.estimates)
    lowers[,l] <- pred$predictions - z * se
    uppers[,l] <- pred$predictions + z * se}
uppest_lrw_bound <- apply(lowers, 1, max)
conf_set_grf <- binary_to_confidence_set(uppers >= uppest_lrw_bound)
list.set.valued.policies[["glb.grf"]] <- conf_set_grf
results_glb_grf <- table.evaluation.real(conf_set_grf, 
                                           prop_score_new = gAW.pred.pseudo.r, 
                                           potential_outcomes = Q.all.pseudo.r,
                                           df_new_sample = df_new_sample,
                                           levels_A = levels_A, 
                                           treatment_name = treatment_name, 
                                           outcome_name = outcome_name)
  
# ── Using linear model with interactions for estimation ─────────────────────
lowers <- uppers <- matrix(0, nrow=nrow(df_new_sample), ncol=m)
for (l in as.numeric(levels_A)){
    data_l <- data.frame(df_new_sample[,covariates_name], A=factor(l, levels = levels_A))
    pred <- stats::predict(glb.model.lm, newdata = data_l, se.fit = TRUE)
    se <- pred$se.fit
    lowers[,l] <- (pred$fit - z * se) %>% as.numeric()
    uppers[,l] <- (pred$fit + z * se) %>% as.numeric()}
uppest_lrw_bound <- apply(lowers, 1, max)
conf_set_lm <- binary_to_confidence_set(uppers>=uppest_lrw_bound)
list.set.valued.policies[["glb.glm"]] <- conf_set_lm
results_glb_lm <- table.evaluation.real(conf_set_lm, 
                                          prop_score_new = gAW.pred.pseudo.r, 
                                          potential_outcomes = Q.all.pseudo.r,
                                          df_new_sample = df_new_sample,
                                          levels_A = levels_A, 
                                          treatment_name = treatment_name, 
                                          outcome_name = outcome_name)
  
  
# ── Collect results for this bootstrap iteration ────────────────────────────
  indices <- c(cardinality = 1, spv_data = 2)
  c(lapply(indices, function(i) {
      do.call(cbind, c(
        list(
          results_glb_grf[[i]],
          results_glb_lm[[i]],
          results.policy.agg[[i]],
          results_agg[[i]]
        ),
        lapply(results.single.policy, `[[`, i),
        lapply(results_policy, `[[`, i)
      ))
    }),
    list(quantiles_and_svp, 
      selected_methods = selected_methods, 
      doptFactorPredict_new_naive = doptFactorPredict_new_naive,
      margin_po_new = margin_po_new, list.set.valued.policies))