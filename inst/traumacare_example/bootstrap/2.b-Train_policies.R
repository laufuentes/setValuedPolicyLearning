train_policies <- function(train_b, train1, calibration, pseudo_test, seed) {
  set.seed(seed)
  # Extract base datasets
  X_b <- train_b[, covariates_name, drop = FALSE]
  A_b <- train_b[, treatment_name]
  Y_b <- train_b[, outcome_name]
  
  X_train <- train1[, covariates_name, drop = FALSE]
  A_train <- train1[, treatment_name]   
  Y_train <- train1[, outcome_name]  
  
  m_levels <- 1:m
  df_new <- df_new_sample
  
  pred_calibration <- list()
  pred_pseudo_data <- list()

  ## Probability forest  ─────────────────────────────────────────────────────────
  proba.forest <- probability_forest(X=cbind(X_train, A_train), 
                                     Y = Y_train %>% as.factor())
  po_calibration <- do.call(cbind, lapply(levels_A %>% as.numeric(), 
                                          function(val) {
                                                     new_data <- calibration[, c(covariates_name, treatment_name)]
                                                     new_data[,treatment_name] <- val
                                                     stats::predict(proba.forest, newdata = new_data)$predictions[,2]}))
  
  pred_calibration[["proba.forest"]] <- apply(po_calibration, 1, function(x)which.max(x) -1 ) %>% as.numeric()
  
  po_pseudo <- do.call(cbind, lapply(levels_A %>% as.numeric(), 
                                              function(val) {
                                                new_data <- pseudo_test[,covariates_name]
                                                new_data[,treatment_name] <- val
                                                stats::predict(proba.forest, newdata = new_data)$predictions[,2]}))
  
  pred_pseudo_data[["proba.forest"]] <- apply(po_pseudo, 1, function(x)which.max(x) -1 ) %>% as.numeric()
  
  ## ── 2. MACF ─────────────────────────────────────────────────────────────────
  multi.forest <- grf::multi_arm_causal_forest(X = X_train, 
                                               Y = Y_train, 
                                               W = A_train %>% as.factor())
  
  preds_cal_macf <- predict(multi.forest, calibration[,covariates_name])$predictions 
  pred_calibration[["MACF"]] <- ifelse(preds_cal_macf[,,1] >0, 1, 0) %>% as.numeric()
  
  preds_pseudo_macf <- predict(multi.forest, pseudo_test[,covariates_name])$predictions 
  pred_pseudo_data[["MACF"]] <- ifelse(preds_pseudo_macf[,,1] >0, 1, 0)%>% as.numeric()
  
  ## ── 3. Policytree & Hybrid ──────────────────────────────────────────────────
  forest <- grf::causal_forest(X = X_train, 
                               Y = Y_train, 
                               W = as.numeric(A_train))
  DR.scores <- policytree::double_robust_scores(forest)
  tree <- policytree::policy_tree(X_train, Gamma= DR.scores)
  
  pred_calibration[["Tree"]] <- stats::predict(tree, newdata = calibration[,covariates_name]) -1  
  pred_pseudo_data[["Tree"]] <- stats::predict(tree, newdata = pseudo_test[,covariates_name]) -1 
  
  hybrid_tree <- policytree::hybrid_policy_tree(X_train, Gamma=DR.scores)
  pred_calibration[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, 
                                                      newdata = calibration[,covariates_name]) - 1 # calibration
  pred_pseudo_data[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = pseudo_test[,covariates_name]) -1 # test 
  
  ## ── 5. Q-learning: SuperLearner ────────────────────────────────────────────
  SL.library_cond <- c("SL.randomForest", "SL.mean", "SL.gam", 
                       "SL.glm", "SL.xgboost")
  
  # train model
  data_A1 <- train1[train1[[treatment_name]] == 1, ]
  data_A0 <- train1[train1[[treatment_name]] == 0, ]
  
  SL_A1 <- SuperLearner::SuperLearner(
    Y = Y_train[train1[[treatment_name]] == 1] %>% as.numeric(),
    X = data_A1[, covariates_name] %>% as.data.frame(),
    SL.library = SL.library_cond, family = binomial())
  
  SL_A0 <- SuperLearner::SuperLearner(
    Y = Y_train[train1[[treatment_name]] == 0] %>% as.numeric(),
    X = data_A0[, covariates_name] %>% as.data.frame(),
    SL.library = SL.library_cond, family = binomial())
  
  cal_preds1 <- predict(SL_A1,newdata = calibration[, covariates_name] %>% as.data.frame())
  cal_preds0 <- predict(SL_A0, newdata = calibration[, covariates_name] %>% as.data.frame())
  pred_calibration[["ql.SL"]] <- ifelse(cal_preds1$pred > cal_preds0$pred, 1,0) %>% as.numeric()
  
  pseudo_preds1 <- stats::predict(SL_A1,newdata = pseudo_test[, covariates_name] %>% as.data.frame())
  pseudo_preds0 <- predict(SL_A0, newdata = pseudo_test[, covariates_name] %>% as.data.frame())
  pred_pseudo_data[["ql.SL"]] <-  ifelse(pseudo_preds1$pred > pseudo_preds0$pred, 1, 0) %>% as.numeric()
  
  ## ── 6. Q-learning : linear model with interactions ─────────────────────────
  ql.glm.interact = stats::glm(formula = formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name)),
                               data = train1, family = "binomial")
  
  po.glm <- do.call(cbind,lapply(levels_A %>% as.numeric(), 
                                 function(val) {
                                   new_data <- calibration[, c(covariates_name, treatment_name)]
                                   new_data[,treatment_name] <- val
                                   stats::predict(ql.glm.interact, newdata = new_data, type = "response")}))
  pred_calibration[["ql.glm.interact"]] <- apply(po.glm, 1, function(x)which.max(x)-1) %>% as.numeric()
  po_pseudo.glm <- do.call(cbind,lapply(levels_A %>% as.numeric(), 
                                           function(val) {
                                             new_data <- pseudo_test[, covariates_name]
                                             new_data[,treatment_name] <- val
                                             stats::predict(ql.glm.interact, newdata = new_data, type = "response")}))
  
  pred_pseudo_data[["ql.glm.interact"]] <-  apply(po_pseudo.glm, 1, function(x)which.max(x)-1) %>% as.numeric()
  
  libraryNames <- names(pred_calibration)
  numalgs <- length(libraryNames)
  
  doptFactorPredict_cal <- do.call(cbind, pred_calibration) %>% as.array()
  colnames(doptFactorPredict_cal) <- libraryNames
  
  doptFactorPredict_pseudo <- do.call(cbind, pred_pseudo_data) %>% as.array()
  colnames(doptFactorPredict_pseudo) <- libraryNames
  
  ## ── 7. Naive method & GLB ──────────────────────────────────────────────────
  cat("training experts (naive version)...")
  pred_pseudo_data_naive <- pred_new_data_naive <- list()
  
  ## Probability forest (also GLB model)  ────────────────────────────────────────
  model.glb.pf <- probability_forest(X=cbind(X_b, A_b), 
                                         Y = Y_b %>% as.factor())
  po_naive <- do.call(cbind, 
                      lapply(levels_A %>% as.numeric(), 
                             function(val) {
                               new_data <- df_new_sample[,covariates_name]
                               new_data[,treatment_name] <- val
                               stats::predict(model.glb.pf, 
                                              newdata = new_data)$predictions[,2]}))
  
  pred_new_data_naive[["proba.forest"]] <- apply(po_naive, 1, function(x)which.max(x) -1 ) %>% as.numeric()
  
  po_pseudo_naive <- do.call(cbind, 
                             lapply(levels_A %>% as.numeric(), 
                                    function(val) {
                                      new_data <-  pseudo_test[,covariates_name]
                                      new_data[,treatment_name] <- val
                                      stats::predict(model.glb.pf, newdata = new_data)$predictions[,2]}))
  
  pred_pseudo_data_naive[["proba.forest"]] <- apply(po_pseudo_naive, 1, function(x)which.max(x) -1 ) %>% as.numeric()
  
  ## ── 2. MACF ─────────────────────────────────────────────────────────────────
  multi.forest <- grf::multi_arm_causal_forest(X = X_b, 
                                               Y = Y_b, 
                                               W = A_b %>% as.factor())
  
  preds_new_macf_naive <- predict(multi.forest, X_new)$predictions 
  pred_new_data_naive[["MACF"]] <- ifelse(preds_new_macf_naive[,,1] >0, 1, 0) %>% as.numeric()
  
  preds_pseudo_macf_naive <- predict(multi.forest, pseudo_test[,covariates_name])$predictions 
  pred_pseudo_data_naive[["MACF"]] <- ifelse(preds_pseudo_macf_naive[,,1] >0, 1, 0) %>% as.numeric()
  
  ## ── 3. Policytree & Hybrid ──────────────────────────────────────────────────
  forest <- grf::causal_forest(X = X_b, Y = Y_b, 
                               W = as.numeric(A_b))
  DR.scores <- policytree::double_robust_scores(forest)
  tree_naive <- policytree::policy_tree(X_b, Gamma= DR.scores)
  
  pred_new_data_naive[["Tree"]] <- stats::predict(tree_naive, newdata = X_new) -1  # test  
  pred_pseudo_data_naive[["Tree"]] <- stats::predict(tree_naive, 
                                                     newdata = pseudo_test[,covariates_name]) -1  # test  
  hybrid_tree_naive <- policytree::hybrid_policy_tree(X_b, Gamma=DR.scores)
  
  pred_new_data_naive[["Hybrid_Tree"]] <- stats::predict(hybrid_tree_naive, newdata = X_new) -1 # test 
  pred_pseudo_data_naive[["Hybrid_Tree"]] <- stats::predict(hybrid_tree_naive, 
                                                            newdata = pseudo_test[,covariates_name]) -1 # test 
  ## ── 5. Q-learning: SuperLearner ────────────────────────────────────────────
  data_A1 <- train_b[train_b[[treatment_name]] == 1, ]
  data_A0 <- train_b[train_b[[treatment_name]] == 0, ]
  
  SL_A1 <- SuperLearner::SuperLearner(
    Y = Y_b[train_b[[treatment_name]] == 1] %>% as.numeric(),
    X = data_A1[, covariates_name] %>% as.data.frame(),
    SL.library = SL.library_cond, family = binomial())
  
  SL_A0 <- SuperLearner::SuperLearner(
    Y = Y_b[train_b[[treatment_name]] == 0] %>% as.numeric(),
    X = data_A0[, covariates_name] %>% as.data.frame(),
    SL.library = SL.library_cond, family = binomial())
  
  new_preds1_naive <- stats::predict(SL_A1,
                                     newdata = df_new_sample[, covariates_name] %>%
                                       as.data.frame())
  
  new_preds0_naive <- predict(SL_A0,
                              newdata = df_new_sample[, covariates_name] %>%
                                as.data.frame())
  
  pred_new_data_naive[["ql.SL"]] <-  ifelse(new_preds1_naive$pred > new_preds0_naive$pred, 1, 0) %>% as.numeric()
  
  pseudo_preds1_naive <- stats::predict(SL_A1,
                                     newdata = pseudo_test[, covariates_name] %>%
                                       as.data.frame())
  pseudo_preds0_naive <- predict(SL_A0,
                              newdata = pseudo_test[, covariates_name] %>%
                                as.data.frame())
  
  pred_pseudo_data_naive[["ql.SL"]] <-  ifelse(pseudo_preds1_naive$pred > pseudo_preds0_naive$pred, 1, 0) %>% as.numeric()
  
  ## ── 6. Q-learning : generalized linear model with interactions ─────────────
  model.glb.glm = stats::glm(formula = formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name)),
                      data = train_b, family = "binomial")
  
  
  po_new <- do.call(cbind,lapply(levels_A %>% as.numeric(), 
                                                 function(val) {
                                                   new_data <- df_new_sample[, covariates_name]
                                                   new_data[,treatment_name] <- val
                                                   stats::predict(ql.glm, newdata = new_data, type = "response")}))
  
  pred_new_data_naive[["ql.glm.interact"]] <-  apply(po_new, 1, function(x)which.max(x)-1) %>% as.numeric()
  
  po_naive <- do.call(cbind,lapply(levels_A %>% as.numeric(), 
                                   function(val) {
                                     new_data <- pseudo_test[, covariates_name]
                                     new_data[,treatment_name] <- val
                                     stats::predict(ql.glm, newdata = new_data, type = "response")}))
  
  pred_pseudo_data_naive[["ql.glm.interact"]] <-  apply(po_naive, 1, function(x)which.max(x)-1) %>% as.numeric()
  
  
  doptFactorPredict_new_naive <- do.call(cbind, pred_new_data_naive) %>% as.array()
  colnames(doptFactorPredict_new_naive) <- libraryNames
  
  doptFactorPredict_pseudo_naive <- do.call(cbind, pred_pseudo_data_naive) %>% as.array()
  colnames(doptFactorPredict_pseudo_naive) <- libraryNames
  
  unweighted_probs_naive <- weighted_probs_experts(fitted_experts = doptFactorPredict_new_naive,
                                                   weights =rep(1/numalgs, numalgs),
                                                   df_pred = df_new_sample,
                                                   levels = as.numeric(levels_A))
  
  unweighted.naive <- apply(
    apply(unweighted_probs_naive, 1, 
          function(x){rmultinom(1, 1, prob=x)}), 2, which.max)-1
  
  unweighted_probs_naive_pseudo <- weighted_probs_experts(fitted_experts = doptFactorPredict_pseudo_naive,
                                                   weights =rep(1/numalgs, numalgs),
                                                   df_pred = pseudo_test,
                                                   levels = as.numeric(levels_A))
  
  unweighted.pseudo.naive <- apply(
    apply(unweighted_probs_naive_pseudo, 1, 
          function(x){rmultinom(1, 1, prob=x)}), 2, which.max)-1
  
  colnames(unweighted_probs_naive) <- as.numeric(levels_A)
  return(list(
    doptFactorPredict_cal = doptFactorPredict_cal, 
    doptFactorPredict_pseudo = doptFactorPredict_pseudo,
    doptFactorPredict_new_naive = doptFactorPredict_new_naive, 
    doptFactorPredict_pseudo_naive = doptFactorPredict_pseudo_naive,
    glb.model.grf = glb.model.grf, 
    glb.model.lm = glb.model.lm, 
    selected_methods = selected_methods, 
    unweighted.naive_new = unweighted.naive,
    unweighted.pseudo.naive=unweighted.pseudo.naive))
}