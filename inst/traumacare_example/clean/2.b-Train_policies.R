set.seed(seed)

# ── Conformal procedure training  ─────────────────────────────────────────────

X_train <- X[SL.out$folds_conformal[[1]],]
A_train <- A[SL.out$folds_conformal[[1]]]
Y_train <- Y[SL.out$folds_conformal[[1]]]
  
pred_calibration <- list()  # replace with predictions on calibration 
pred_new_data <- list() # replace with predictions on new
pred_pseudo_data <- list()
cat("training experts (conformal prediction)...")

## Probability forest  ─────────────────────────────────────────────────────────
proba.forest <- probability_forest(X=cbind(X_train, A_train), 
                                   Y = Y_train %>% as.factor())

# create predictions for counterfactuals 
potential_outcomes_calibration <- do.call(cbind,
                                          lapply(levels_A %>% as.numeric(), 
                                                 function(val) {
                                                   new_data <- calibration[, c(covariates_name, treatment_name)]
                                                   new_data[,treatment_name] <- val
                                                   stats::predict(proba.forest, newdata = new_data)$predictions[,2]}))

pred_calibration[["proba.forest"]] <- apply(potential_outcomes_calibration, 1, function(x)which.max(x) -1 ) %>% as.numeric()

potential_outcomes_new <- do.call(cbind, 
                                  lapply(levels_A %>% as.numeric(), 
                                         function(val) {
                                           new_data <- SL.out$df_new_sample[,covariates_name]
                                           new_data[,treatment_name] <- val
                                           stats::predict(proba.forest, newdata = new_data)$predictions[,2]}))

pred_new_data[["proba.forest"]] <- apply(potential_outcomes_new, 1, function(x)which.max(x) -1 ) %>% as.numeric()

potential_outcomes_pseudo <- do.call(cbind, 
                                  lapply(levels_A %>% as.numeric(), 
                                         function(val) {
                                           new_data <- pseudo.test.predict[,covariates_name]
                                           new_data[,treatment_name] <- val
                                           stats::predict(proba.forest, newdata = new_data)$predictions[,2]}))

pred_pseudo_data[["proba.forest"]] <- apply(potential_outcomes_pseudo, 1, function(x)which.max(x) -1 ) %>% as.numeric()

## Multi arm causal forest (MACF required for training trees)  ─────────────────
multi.forest <- grf::multi_arm_causal_forest(X = X_train, 
                                             Y = Y_train, 
                                             W = A_train %>% as.factor())

preds_cal_macf <- predict(multi.forest, calibration[,covariates_name])$predictions 
pred_calibration[["MACF"]] <- ifelse(preds_cal_macf[,,1] >0, 1, 0) %>% as.numeric()

preds_new_macf <- predict(multi.forest, X_new)$predictions 
pred_new_data[["MACF"]] <- ifelse(preds_new_macf[,,1] >0, 1, 0)%>% as.numeric()

preds_pseudo_macf <- predict(multi.forest, pseudo.test.predict[,covariates_name])$predictions 
pred_pseudo_data[["MACF"]] <- ifelse(preds_pseudo_macf[,,1] >0, 1, 0)%>% as.numeric()

## policytree  ─────────────────────────────────────────────────────────────────
forest <- grf::causal_forest(X = X_train, 
                             Y = Y_train, 
                             W = as.numeric(A_train))
DR.scores <- policytree::double_robust_scores(forest)
tree <- policytree::policy_tree(X_train, Gamma= DR.scores)

pred_calibration[["Tree"]] <- stats::predict(tree, 
                                             newdata = calibration[,covariates_name]) -1  # calibration  
pred_new_data[["Tree"]] <- stats::predict(tree, newdata = X_new) -1  # test  
pred_pseudo_data[["Tree"]] <- stats::predict(tree, newdata = pseudo.test.predict[,covariates_name]) -1  # test  


# hybrid policytree  ───────────────────────────────────────────────────────────
hybrid_tree <- policytree::hybrid_policy_tree(X_train, Gamma=DR.scores)

pred_calibration[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, 
                                                    newdata = calibration[,covariates_name]) - 1 # calibration
pred_new_data[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = X_new) -1 # test 
pred_pseudo_data[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = pseudo.test.predict[,covariates_name]) -1 # test 

# Many experts with SuperLearner (Q-learning)   ────────────────────────────────
# Note: You can write a custom function to pre-resample X and Y within a SuperLearner template,
# or use the built-in sub-sampling wrappers if available.
SL.library_cond <- c("SL.randomForest", "SL.mean", "SL.gam", 
                      "SL.glm", "SL.xgboost")
 
# train model
data_A1 <- train1[train1[[treatment_name]] == 1, ]
data_A0 <- train1[train1[[treatment_name]] == 0, ]
 
SL_A1 <- SuperLearner::SuperLearner(
   Y = Y_train[train1[[treatment_name]] == 1] %>% as.numeric(),
   X = data_A1[, covariates_name] %>% as.data.frame(),
   SL.library = SL.library_cond,
   family = binomial())
 
SL_A0 <- SuperLearner::SuperLearner(
   Y = Y_train[train1[[treatment_name]] == 0] %>% as.numeric(),
   X = data_A0[, covariates_name] %>% as.data.frame(),
   SL.library = SL.library_cond,
   family = binomial())

cal_preds1 <- predict(SL_A1,
                      newdata = calibration[, covariates_name] %>%
                        as.data.frame())

cal_preds0 <- predict(SL_A0, 
                      newdata = calibration[, covariates_name] %>%
                        as.data.frame())
pred_calibration[["ql.SL"]] <- ifelse(cal_preds1$pred > cal_preds0$pred, 1,0) %>% as.numeric()

new_preds1 <- stats::predict(SL_A1,
                             newdata = SL.out$df_new_sample[, covariates_name] %>%
                               as.data.frame())

new_preds0 <- predict(SL_A0,
                       newdata = SL.out$df_new_sample[, covariates_name] %>%
                         as.data.frame())
 
pred_new_data[["ql.SL"]] <-  ifelse(new_preds1$pred > new_preds0$pred, 1, 0) %>% as.numeric()

pseudo_preds1 <- stats::predict(SL_A1,
                             newdata = pseudo.test.predict[, covariates_name] %>%
                               as.data.frame())

pseudo_preds0 <- predict(SL_A0,
                      newdata = pseudo.test.predict[, covariates_name] %>%
                        as.data.frame())

pred_pseudo_data[["ql.SL"]] <-  ifelse(pseudo_preds1$pred > pseudo_preds0$pred, 1, 0) %>% as.numeric()

# Q-learning: glm with interactions ────────────────────────────────────────────
ql.glm = stats::glm(formula = formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name)),
                  data = train1, family = "binomial")

potential_outcomes_cal.ql.glm <- do.call(cbind,lapply(levels_A %>% as.numeric(), 
                                                function(val) {
                                                  new_data <- calibration[, c(covariates_name, treatment_name)]
                                                  new_data[,treatment_name] <- val
                                                  stats::predict(ql.glm, newdata = new_data, type = "response")}))

pred_calibration[["ql.lm.interact"]] <- apply(potential_outcomes_cal.ql.glm, 1, function(x)which.max(x)-1) %>% as.numeric()

potential_outcomes_new.ql.glm <- do.call(cbind,lapply(levels_A %>% as.numeric(), 
                                               function(val) {
                                                 new_data <- SL.out$df_new_sample[, covariates_name]
                                                 new_data[,treatment_name] <- val
                                                 stats::predict(ql.glm, newdata = new_data, type = "response")}))

pred_new_data[["ql.glm.interact"]] <-  apply(potential_outcomes_new.ql.glm, 1, function(x)which.max(x)-1) %>% as.numeric()

potential_outcomes_pseudo.ql.glm <- do.call(cbind,lapply(levels_A %>% as.numeric(), 
                                                      function(val) {
                                                        new_data <- pseudo.test.predict[, covariates_name]
                                                        new_data[,treatment_name] <- val
                                                        stats::predict(ql.glm, newdata = new_data, type = "response")}))

pred_pseudo_data[["ql.glm.interact"]] <-  apply(potential_outcomes_pseudo.ql.glm, 1, function(x)which.max(x)-1) %>% as.numeric()


SL.out$libraryNames <- names(pred_calibration)
numalgs <- length(SL.out$libraryNames)

SL.out$doptFactorPredict_test <- do.call(cbind, pred_calibration) %>% as.array()
colnames(SL.out$doptFactorPredict_test) <- SL.out$libraryNames

SL.out$doptFactorPredict_new <- do.call(cbind, pred_new_data) %>% as.array()
colnames(SL.out$doptFactorPredict_new) <- SL.out$libraryNames

SL.out$doptFactorPredict_pseudo <- do.call(cbind, pred_pseudo_data) %>% as.array()
colnames(SL.out$doptFactorPredict_pseudo) <- SL.out$libraryNames

# ── Naive procedure training  ─────────────────────────────────────────────
training_data <- rbind(train1, train2, calibration)

X_train <- X[SL.out$folds[[1]],]
A_train <- A[SL.out$folds[[1]]]
Y_train <- Y[SL.out$folds[[1]]]
cat("training experts (naive version)...")

pred_new_data_naive <- list() # replace with predictions on new

## Probability forest (also GLB model)  ────────────────────────────────────────
SL.out$model.glb <- probability_forest(X=cbind(X_train, A_train), 
                                   Y = Y_train %>% as.factor())

# create predictions for counterfactuals 
potential_outcomes_new_naive <- do.call(cbind, 
                                  lapply(levels_A %>% as.numeric(), 
                                         function(val) {
                                           new_data <- SL.out$df_new_sample[,covariates_name]
                                           new_data[,treatment_name] <- val
                                           stats::predict(SL.out$model.glb, newdata = new_data)$predictions[,2]}))

pred_new_data_naive[["proba.forest"]] <- apply(potential_outcomes_new_naive, 1, function(x)which.max(x) -1 ) %>% as.numeric()

## Multi arm causal forest (MACF required for training trees)  ─────────────────
multi.forest <- grf::multi_arm_causal_forest(X = X_train, 
                                             Y = Y_train, 
                                             W = A_train %>% as.factor())

preds_new_macf_naive <- predict(multi.forest, X_new)$predictions 
pred_new_data_naive[["MACF"]] <- ifelse(preds_new_macf_naive[,,1] >0, 1, 0) %>% 
  as.numeric()

## policytree  ─────────────────────────────────────────────────────────────────
forest <- grf::causal_forest(X = X_train, 
                             Y = Y_train, 
                             W = as.numeric(A_train))
DR.scores <- policytree::double_robust_scores(forest)
tree_naive <- policytree::policy_tree(X_train, Gamma= DR.scores)

pred_new_data_naive[["Tree"]] <- stats::predict(tree_naive, newdata = X_new) -1  # test  


# hybrid policytree  ───────────────────────────────────────────────────────────
hybrid_tree_naive <- policytree::hybrid_policy_tree(X_train, Gamma=DR.scores)

pred_new_data_naive[["Hybrid_Tree"]] <- stats::predict(hybrid_tree_naive, newdata = X_new) -1 # test 


# Many experts with SuperLearner (Q-learning)   ────────────────────────────────
# train model
data_A1 <- training_data[training_data[[treatment_name]] == 1, ]
data_A0 <- training_data[training_data[[treatment_name]] == 0, ]

SL_A1 <- SuperLearner::SuperLearner(
  Y = Y_train[training_data[[treatment_name]] == 1] %>% as.numeric(),
  X = data_A1[, covariates_name] %>% as.data.frame(),
  SL.library = SL.library_cond,
  family = binomial())

SL_A0 <- SuperLearner::SuperLearner(
  Y = Y_train[training_data[[treatment_name]] == 0] %>% as.numeric(),
  X = data_A0[, covariates_name] %>% as.data.frame(),
  SL.library = SL.library_cond,
  family = binomial())

new_preds1_naive <- stats::predict(SL_A1,
                             newdata = SL.out$df_new_sample[, covariates_name] %>%
                               as.data.frame())

new_preds0_naive <- predict(SL_A0,
                      newdata = SL.out$df_new_sample[, covariates_name] %>%
                        as.data.frame())

pred_new_data_naive[["ql.SL"]] <-  ifelse(new_preds1_naive$pred > new_preds0_naive$pred, 1, 0) %>% as.numeric()

# Q-learning: glm with interactions ────────────────────────────────────────────
ql.glm = stats::glm(formula = formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name)),
                    data = training_data, family = "binomial")


potential_outcomes_new <- do.call(cbind,lapply(levels_A %>% as.numeric(), 
                                               function(val) {
                                                 new_data <- SL.out$df_new_sample[, covariates_name]
                                                 new_data[,treatment_name] <- val
                                                 stats::predict(ql.glm, newdata = new_data, type = "response")}))

pred_new_data_naive[["ql.glm.interact"]] <-  apply(potential_outcomes_new, 1, function(x)which.max(x)-1) %>% as.numeric()

  
doptFactorPredict_new_naive <- do.call(cbind, pred_new_data_naive) %>% as.array()
colnames(doptFactorPredict_new_naive) <- SL.out$libraryNames
  
# Generate the distribution of unweighted experts 
unweighted_probs_naive <- weighted_probs_experts(fitted_experts = doptFactorPredict_new_naive,
                                             weights =rep(1/numalgs, numalgs),
                                             df_pred = SL.out$df_new_sample,
                                             levels = as.numeric(levels_A))
  
colnames(unweighted_probs_naive) <- as.numeric(levels_A)
write.csv((doptFactorPredict_new_naive) %>% 
              as.data.frame() %>% 
              mutate("SUBJECT_REF"= SL.out$df_new_sample$SUBJECT_REF), 
          file=paste0("inst/traumacare_example/intermediate/Experts_aggregation_", outcome_name, ".csv"))
