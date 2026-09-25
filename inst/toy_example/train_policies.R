set.seed(seed)

pred_calibration <- list()  # replace with predictions on calibration 
pred_new_data <- list() # replace with predictions on new

cat("training experts (conformal prediction)...")


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


SL.out$libraryNames <- names(pred_calibration)
numalgs <- length(SL.out$libraryNames)

SL.out$doptFactorPredict_test <- do.call(cbind, pred_calibration) %>% as.array()
colnames(SL.out$doptFactorPredict_test) <- SL.out$libraryNames

SL.out$doptFactorPredict_new <- do.call(cbind, pred_new_data) %>% as.array()
colnames(SL.out$doptFactorPredict_new) <- SL.out$libraryNames

###############################################################################  
############################  Naive method #################################### 
############################################################################### 
training_data <- rbind(train1, train2, test)
X_train <- X
A_train <- A
Y_train <- Y
cat("training experts (naive version)...")

pred_new_data_naive <- list() # replace with predictions on new

## Q-learning: lm interaction ──────────────────────────────────────────────────
ql.lm = stats::lm(formula = formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name)),
                  data = training_data)

potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- SL.out$df_new_sample[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  stats::predict(ql.lm, newdata = new_data)}))
pred_new_data_naive[["ql.lm.interact"]] <- apply(potential_outcomes_new, 1, which.max) %>% as.numeric() %>% as.numeric()

## Q-learning: regression forest  ──────────────────────────────────────────────
ql.reg.forest <- grf::regression_forest(X = cbind(X_train, A_train), 
                              Y = Y_train, seed = seed)
potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- cbind(X_new, factor(val, levels=levels_A) %>% as.numeric())
  stats::predict(ql.reg.forest, newdata = new_data)$predictions}))
pred_new_data_naive[["ql.reg.forest"]] <- apply(potential_outcomes_new, 1, which.max)

## MACF ────────────────────────────────────────────────────────────────────────
forest <- grf::multi_arm_causal_forest(X = X_train, 
                                       Y = Y_train, 
                                       W = A_train %>% as.factor())
  
cate_new_matrix <- predict(forest, X_new)$predictions[, , 1]
cate_new_matrix_all <- cbind(0, cate_new_matrix)
colnames(cate_new_matrix_all) <- levels_A
pred_new_data_naive[["MACF"]] <- apply(cate_new_matrix_all, 1, which.max)
  
## policytree  ─────────────────────────────────────────────────────────────────
DR.scores <- policytree::double_robust_scores(forest)
tree <- policytree::policy_tree(X_train, Gamma= DR.scores)

pred_new_data_naive[["Tree"]] <- stats::predict(tree, X_new)  # test  

# hybrid policytree  ───────────────────────────────────────────────────────────
hybrid_tree <- policytree::hybrid_policy_tree(X_train, Gamma=DR.scores)

pred_new_data_naive[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = X_new) # test 


# Q-learning: SuperLearner   ───────────────────────────────────────────────────

QL_mod <- SuperLearner::SuperLearner( # Outcome model
  Y=training_data[,outcome_name], X = training_data[,c(covariates_name,treatment_name)],
  SL.library=SL.library_cond, family = "gaussian")

potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- SL.out$df_new_sample[, c(covariates_name, treatment_name)]
  new_data[,treatment_name] <- factor(val, levels=levels_A)
  SuperLearner::predict.SuperLearner(QL_mod, newdata = new_data)$pred}))

pred_new_data_naive[["ql.SL"]] <- apply(potential_outcomes_new, 1, which.max)

# POLLE package  ───────────────────────────────────────────────────────────────
train_polle <- data.table::as.data.table(training_data)
train_polle[[treatment_name]] <- as.factor(train_polle[[treatment_name]])
pd_train <- polle::policy_data(train_polle, action = treatment_name,
                               covariates = covariates_name, utility = outcome_name)

po.drql.lm <- l_obj.drql.glm(
  policy_data = pd_train, 
  q_models = q_glm(),
  g_models = g_rf())

pred_new_data_naive[["drql.lm"]] <- as.numeric(polle::get_policy(po.drql.lm)(pd_new)$d)

po.drql.rf <- l_obj.drql.rf(
  policy_data = pd_train, 
  q_models = q_rf(),
  g_models = g_rf())

pred_new_data_naive[["drql.rf"]] <- as.numeric(polle::get_policy(po.drql.rf)(pd_new)$d)

po.drql.xgboost <- l_obj.drql.xgboost(
  policy_data = pd_train, 
  q_models = q_sl(),
  g_models = g_rf())

pred_new_data_naive[["drql.xgboost"]] <- as.numeric(polle::get_policy(po.drql.xgboost)(pd_new)$d)

po.drql.ksvm <- l_obj.drql.ksvm(
  policy_data = pd_train, 
  q_models = q_sl(),
  g_models = g_rf())

pred_new_data_naive[["drql.ksvm"]] <- as.numeric(polle::get_policy(po.drql.ksvm)(pd_new)$d)

libraryNames_naive <- names(pred_new_data_naive)
numalgs_naive <- length(libraryNames_naive)

doptFactorPredict_new_naive <- do.call(cbind, pred_new_data_naive) %>% as.array()
colnames(doptFactorPredict_new_naive) <- libraryNames_naive

# Generate the distribution of unweighted experts 
unweighted_probs_naive <- weighted_probs_experts(fitted_experts = doptFactorPredict_new_naive,
                                           weights =rep(1/numalgs_naive, numalgs_naive),
                                           df_pred = SL.out$df_new_sample,
                                           levels = as.numeric(levels_A))


SL.out$unweighted_new_naive <- apply(apply(unweighted_probs_naive, 1, function(x){
  rmultinom(1,1,prob=x)}),2, which.max)
