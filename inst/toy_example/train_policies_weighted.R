set.seed(seed)

pred_calibration <- list()  # replace with predictions on calibration 
pred_new_data <- list() # replace with predictions on new

cat("training experts (conformal prediction)...")

X_train <- X[SL.out$folds[[1]],]  #X_train <- X[SL.out$folds[[1]],] 
A_train <- A[SL.out$folds[[1]]]   #A_train <- A[SL.out$folds[[1]]]
Y_train <- Y[SL.out$folds[[1]]]   #Y_train <- Y[SL.out$folds[[1]]]


## Q-learning: lm interaction  ─────────────────────────────────────────────────
ql.lm = stats::lm(formula = formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name)),
                  data = train1)

potential_outcomes_test <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- test[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  stats::predict(ql.lm, newdata = new_data)}))

pred_calibration[["ql.lm.interact"]] <- apply(potential_outcomes_test, 1, which.max) %>% as.numeric()

potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- SL.out$df_new_sample[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  stats::predict(ql.lm, newdata = new_data)}))
pred_new_data[["ql.lm.interact"]] <- apply(potential_outcomes_new, 1, which.max)

## Q-learning: regression forest  ───────────────────────────────────────────────
ql.reg.forest <- grf::regression_forest(X = cbind(X_train, A_train), 
                                        Y = Y_train, seed = seed)
potential_outcomes_test <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- cbind(X[SL.out$folds[[3]],], factor(val, levels=levels_A) %>% as.numeric())
  stats::predict(ql.reg.forest, newdata = new_data)$predictions}))

pred_calibration[["ql.reg.forest"]] <- apply(potential_outcomes_test, 1, which.max)

potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- cbind(X_new, factor(val, levels=levels_A) %>% as.numeric())
  stats::predict(ql.reg.forest, newdata = new_data)$predictions}))

pred_new_data[["ql.reg.forest"]] <- apply(potential_outcomes_new, 1, which.max)

## MACF ────────────────────────────────────────────────────────────────────────
forest <- grf::multi_arm_causal_forest(X = X_train, 
                                       Y = Y_train, 
                                       W = A_train %>% as.factor()) 
cate_matrix <- predict(forest, X[SL.out$folds[[3]], ])$predictions[, , 1]
cate_matrix_all <- cbind(0, cate_matrix)
colnames(cate_matrix_all) <- levels_A
pred_calibration[["MACF"]] <- apply(cate_matrix_all, 1, which.max)

cate_new_matrix <- predict(forest, X_new)$predictions[, , 1]
cate_new_matrix_all <- cbind(0, cate_new_matrix)
colnames(cate_new_matrix_all) <- levels_A
pred_new_data[["MACF"]] <- apply(cate_new_matrix_all, 1, which.max) 

## policytree  ─────────────────────────────────────────────────────────────────
DR.scores <- policytree::double_robust_scores(forest)
tree <- policytree::policy_tree(X_train, Gamma= DR.scores)

pred_calibration[["Tree"]] <- stats::predict(tree, 
                                             newdata = X[SL.out$folds[[3]],])
pred_new_data[["Tree"]] <- stats::predict(tree, X_new)  

# hybrid policytree  ─────────────────────────────────────────────────────────
hybrid_tree <- policytree::hybrid_policy_tree(X_train, Gamma=DR.scores, depth = 3)

pred_calibration[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, 
                                                    newdata = X[SL.out$folds[[3]],]) # calibration
pred_new_data[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = X_new) # test 

# Q-learning: SuperLearner   ───────────────────────────────────────────────────
SL.library_cond <- c("SL.randomForest", "SL.ksvm", "SL.mean",
                     "SL.glm", "SL.xgboost")

QL_mod <- SuperLearner::SuperLearner( 
  Y=train1[,outcome_name], X = train1[,c(covariates_name,treatment_name)],
  SL.library=SL.library_cond, family = "gaussian")

potential_outcomes_test <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- test[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  SuperLearner::predict.SuperLearner(QL_mod, newdata = new_data)$pred}))

pred_calibration[["ql.SL"]] <- apply(potential_outcomes_test, 1, which.max)

potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- SL.out$df_new_sample[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  SuperLearner::predict.SuperLearner(QL_mod, newdata = new_data)$pred}))

pred_new_data[["ql.SL"]] <- apply(potential_outcomes_new, 1, which.max)

# POLLE package  ───────────────────────────────────────────────────────────────
train_polle <- data.table::as.data.table(train1)
train_polle[[treatment_name]] <- as.factor(train_polle[[treatment_name]])
pd_train <- polle::policy_data(train_polle, action = treatment_name,
                               covariates = covariates_name, utility = outcome_name)

test_polle  <- data.table::as.data.table(test)
test_polle[[treatment_name]]  <- as.factor(test_polle[[treatment_name]])
pd_test  <- polle::policy_data(test_polle,  action = treatment_name, 
                               covariates = covariates_name, utility = outcome_name)

new_polle <- data.table::as.data.table(SL.out$df_new_sample)
new_polle[[treatment_name]] <- as.factor(new_polle[[treatment_name]])
pd_new <- polle::policy_data(new_polle, action = treatment_name, 
                             covariates = covariates_name, 
                             utility = outcome_name)

formula_obj <- stats::reformulate(covariates_name)

# lm 
l_obj.drql.glm <- polle::policy_learn(
  type    = "drql",
  control = polle::control_drql(qv_models = polle::q_glm(formula = formula_obj, family = stats::gaussian())), 
  cross_fit_g_models = FALSE)

po.drql.lm <- l_obj.drql.glm(
  policy_data = pd_train, 
  q_models = q_glm(),
  g_models = g_rf())

pred_calibration[["drql.lm"]] <- as.numeric(polle::get_policy(po.drql.lm)(pd_test)$d)
pred_new_data[["drql.lm"]] <- as.numeric(polle::get_policy(po.drql.lm)(pd_new)$d)

l_obj.drql.rf <- polle::policy_learn(
  type    = "drql",
  control = polle::control_drql(qv_models = polle::q_rf(formula = formula_obj)), 
  cross_fit_g_models = FALSE)

po.drql.rf <- l_obj.drql.rf(
  policy_data = pd_train, 
  q_models = q_rf(),
  g_models = g_rf())

pred_calibration[["drql.rf"]] <- as.numeric(polle::get_policy(po.drql.rf)(pd_test)$d)
pred_new_data[["drql.rf"]] <- as.numeric(polle::get_policy(po.drql.rf)(pd_new)$d)

l_obj.drql.xgboost <- polle::policy_learn(
  type    = "drql",
  control = polle::control_drql(qv_models = polle::q_sl(formula = formula_obj, SL.library = "SL.xgboost")), 
  cross_fit_g_models = FALSE)

po.drql.xgboost <- l_obj.drql.xgboost(
  policy_data = pd_train, 
  q_models = q_sl(),
  g_models = g_rf())

pred_calibration[["drql.xgboost"]] <- as.numeric(polle::get_policy(po.drql.xgboost)(pd_test)$d)
pred_new_data[["drql.xgboost"]] <- as.numeric(polle::get_policy(po.drql.xgboost)(pd_new)$d)

l_obj.drql.ksvm <- polle::policy_learn(
  type    = "drql",
  control = polle::control_drql(qv_models = polle::q_sl(formula = formula_obj, 
                                                        SL.library = "SL.ksvm")), 
  cross_fit_g_models = FALSE)

po.drql.ksvm <- l_obj.drql.ksvm(
  policy_data = pd_train, 
  q_models = q_sl(),
  g_models = g_rf())

pred_calibration[["drql.ksvm"]] <- as.numeric(polle::get_policy(po.drql.ksvm)(pd_test)$d)
pred_new_data[["drql.ksvm"]] <- as.numeric(polle::get_policy(po.drql.ksvm)(pd_new)$d)


libraryNames <- names(pred_calibration)
numalgs <- length(libraryNames)

doptFactorPredict_test <- do.call(cbind, pred_calibration) %>% as.array()
colnames(doptFactorPredict_test) <- libraryNames

doptFactorPredict_new <- do.call(cbind, pred_new_data) %>% as.array()
colnames(doptFactorPredict_new) <- libraryNames

libraryNames <- names(pred_calibration)
numalgs <- length(libraryNames)

doptFactorPredict_test_weighted <- do.call(cbind, pred_calibration) %>% as.array()
doptFactorPredict_new_weighted <- do.call(cbind, pred_new_data) %>% as.array()

colnames(doptFactorPredict_test_weighted) <- libraryNames
colnames(doptFactorPredict_new_weighted) <- libraryNames


mod_y <-grf::regression_forest(X = cbind(
  X[SL.out$folds[[4]],], A[SL.out$folds[[4]]]), 
  Y = Y[SL.out$folds[[4]]], seed = seed)

Q.all.actions <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- cbind(X[SL.out$folds[[3]],], factor(val, levels=levels_A) %>% as.numeric())
  stats::predict(mod_y, newdata = new_data)$predictions}))

mod_ps <-  grf::probability_forest(X = SL.out$df_obs[SL.out$folds[[4]],covariates_name],
                                   Y = SL.out$df_obs[SL.out$folds[[4]],treatment_name])

gAW.pred <- stats::predict(mod_ps,
                           newdata = test[, covariates_name, drop = FALSE],
                           type = "prob")$pred

weights <- apply(doptFactorPredict_test_weighted, 2, function(x)set_policy_value(x, test,
                                                                                 covariates = covariates_name,
                                                                                 treatment_name = treatment_name,
                                                                                 outcome_name = outcome_name, gAW.pred.spv =NULL,
                                                                                 Q_all_actions = Q.all.actions, 
                                                                                 gAW.pred = gAW.pred, ab, n_test = 1, levels_A))

# Generate the distribution of unweighted experts 
weighted_probs <- weighted_probs_experts(fitted_experts = doptFactorPredict_test_weighted,
                                         weights = weights/sum(weights),
                                         df_pred = test,
                                         levels = as.numeric(levels_A))

weighted_aggregation <- apply(apply(weighted_probs, 1, function(x){
  rmultinom(1,1,prob=x)}),2, which.max)

weighted_probs_new <- weighted_probs_experts(fitted_experts = doptFactorPredict_new_weighted,
                                         weights = weights/sum(weights),
                                         df_pred = SL.out$df_new_sample,
                                         levels = as.numeric(levels_A))

weighted_aggregation_new <- apply(apply(weighted_probs_new, 1, function(x){
  rmultinom(1,1,prob=x)}),2, which.max)
