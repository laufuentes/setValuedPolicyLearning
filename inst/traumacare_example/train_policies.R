train_policies <- function(train_b, train1, calibration, pseudo_test, seed, SL.library_cond=c("SL.mean", "SL.glm")) {
  set.seed(seed)
  cat("training conformal...")
  # Extract base datasets
  X_b <- train_b[, covariates_name, drop = FALSE]
  A_b <- train_b[, treatment_name]
  Y_b <- train_b[, outcome_name]
  
  X_train <- train1[, covariates_name, drop = FALSE]
  A_train <- train1[, treatment_name]   
  Y_train <- train1[, outcome_name]  
  
  numeric_levels_A <- as.numeric(levels_A)
  pred_calibration <- list()
  pred_pseudo_data  <- list()
  
  ## 1. Probability forest  ────────────────────────────────────────────────────
  proba.forest <- probability_forest(
    X = cbind(X_train, A_train), 
    Y = as.factor(Y_train))
  
  cal_template <- calibration[, c(covariates_name, treatment_name)]
  po_calibration <- sapply(numeric_levels_A, function(val) {
    cal_template[[treatment_name]] <- val
    stats::predict(proba.forest, newdata = cal_template)$predictions[, 2]})
  
  pred_calibration[["proba.forest"]] <- max.col(po_calibration) - 1
  
  pseudo_template <- pseudo_test[, covariates_name, drop = FALSE]
  po_pseudo <- sapply(numeric_levels_A, function(val) {
    pseudo_template[[treatment_name]] <- val
    stats::predict(proba.forest, newdata = pseudo_template)$predictions[, 2]
  })
  pred_pseudo_data[["proba.forest"]] <- max.col(po_pseudo) - 1
  
  ## ── 2. MACF ─────────────────────────────────────────────────────────────────
  multi.forest <- grf::multi_arm_causal_forest(X = X_train, Y = Y_train, 
                                               W = as.factor(A_train))
  
  preds_cal_macf <- predict(multi.forest, calibration[, covariates_name])$predictions 
  pred_calibration[["MACF"]] <- as.numeric(preds_cal_macf[, , 1] > 0)
  
  preds_pseudo_macf <- predict(multi.forest, pseudo_test[, covariates_name])$predictions 
  pred_pseudo_data[["MACF"]] <- as.numeric(preds_pseudo_macf[, , 1] > 0)
  
  ## ── 3. Policytree & Hybrid ──────────────────────────────────────────────────
  forest <- grf::causal_forest(X = X_train, Y = Y_train, W = as.numeric(A_train))
  DR.scores <- policytree::double_robust_scores(forest)
  
  tree <- policytree::policy_tree(X_train, Gamma = DR.scores)
  pred_calibration[["Tree"]] <- stats::predict(tree, newdata = calibration[, covariates_name]) - 1  
  pred_pseudo_data[["Tree"]] <- stats::predict(tree, newdata = pseudo_test[, covariates_name]) - 1 
  
  hybrid_tree <- policytree::hybrid_policy_tree(X_train, Gamma = DR.scores)
  pred_calibration[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = calibration[, covariates_name]) - 1 
  pred_pseudo_data[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = pseudo_test[, covariates_name]) - 1
  
  ## ── 4. Q-learning: SuperLearner ────────────────────────────────────────────
  idx_A1 <- train1[[treatment_name]] == 1
  idx_A0 <- train1[[treatment_name]] == 0
  
  SL_A1 <- SuperLearner::SuperLearner(
    Y = as.numeric(Y_train[idx_A1]),
    X = train1[idx_A1, covariates_name, drop = FALSE],
    SL.library = SL.library_cond, family = binomial())
  
  SL_A0 <- SuperLearner::SuperLearner(
    Y = as.numeric(Y_train[idx_A0]),
    X = train1[idx_A0, covariates_name, drop = FALSE],
    SL.library = SL.library_cond, family = binomial())
  
  cal_df <- as.data.frame(calibration[, covariates_name, drop = FALSE])
  cal_preds1 <- predict(SL_A1, newdata = cal_df)$pred
  cal_preds0 <- predict(SL_A0, newdata = cal_df)$pred
  pred_calibration[["ql.SL"]] <- as.numeric(cal_preds1 > cal_preds0)
  
  pseudo_df <- as.data.frame(pseudo_test[, covariates_name, drop = FALSE])
  pseudo_preds1 <- predict(SL_A1, newdata = pseudo_df)$pred
  pseudo_preds0 <- predict(SL_A0, newdata = pseudo_df)$pred
  pred_pseudo_data[["ql.SL"]] <- as.numeric(pseudo_preds1 > pseudo_preds0)
  
  ## ── 5. Q-learning : linear model with interactions ─────────────────────────
  form_str <- formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name))
  ql.glm.interact <- stats::glm(formula = form_str, data = train1, family = "binomial")
  
  po.glm <- sapply(numeric_levels_A, function(val) {
    cal_template[[treatment_name]] <- val
    stats::predict(ql.glm.interact, newdata = cal_template, type = "response")
  })
  pred_calibration[["ql.glm.interact"]] <- max.col(po.glm) - 1
  
  po_pseudo.glm <- sapply(numeric_levels_A, function(val) {
    pseudo_template[[treatment_name]] <- val
    stats::predict(ql.glm.interact, newdata = pseudo_template, type = "response")
  })
  pred_pseudo_data[["ql.glm.interact"]] <- max.col(po_pseudo.glm) - 1
  
  # Group expert predictions
  libraryNames <- names(pred_calibration)
  numalgs <- length(libraryNames)
  
  doptFactorPredict_cal <- do.call(cbind, pred_calibration)
  doptFactorPredict_pseudo <- do.call(cbind, pred_pseudo_data)
  
  ## ── 6. Naive Method & GLB Models ───────────────────────────────────────────
  cat("training baselines & GLB...")
  
  pred_pseudo_data_naive <- list()
  pred_new_data_naive    <- list()
  
  # Probability Forest (GLB Model)
  model.glb.pf <- probability_forest(
    X = cbind(X_b, A_b), 
    Y = as.factor(Y_b))
  
  df_new_template <- df_new_sample[, covariates_name, drop = FALSE]
  po_naive <- sapply(numeric_levels_A, function(val) {
    df_new_template[[treatment_name]] <- val
    stats::predict(model.glb.pf, newdata = df_new_template)$predictions[, 2]
  })
  pred_new_data_naive[["proba.forest"]] <- max.col(po_naive) - 1
  
  po_pseudo_naive <- sapply(numeric_levels_A, function(val) {
    pseudo_template[[treatment_name]] <- val
    stats::predict(model.glb.pf, newdata = pseudo_template)$predictions[, 2]
  })
  pred_pseudo_data_naive[["proba.forest"]] <- max.col(po_pseudo_naive) - 1
  
  ## ── 2. MACF ─────────────────────────────────────────────────────────────────
  multi.forest_b <- grf::multi_arm_causal_forest(
    X = X_b, Y = Y_b, W = as.factor(A_b)
  )
  X_new <- df_new_sample[, covariates_name, drop = FALSE]
  
  preds_new_macf_naive <- predict(multi.forest_b, X_new)$predictions 
  pred_new_data_naive[["MACF"]] <- as.numeric(preds_new_macf_naive[, , 1] > 0)
  
  preds_pseudo_macf_naive <- predict(multi.forest_b, pseudo_test[, covariates_name])$predictions 
  pred_pseudo_data_naive[["MACF"]] <- as.numeric(preds_pseudo_macf_naive[, , 1] > 0)
  
  ## ── 3. Policytree & Hybrid ──────────────────────────────────────────────────
  forest_b <- grf::causal_forest(X = X_b, Y = Y_b, W = as.numeric(A_b))
  DR.scores_b <- policytree::double_robust_scores(forest_b)
  
  tree_naive <- policytree::policy_tree(X_b, Gamma = DR.scores_b)
  pred_new_data_naive[["Tree"]] <- stats::predict(tree_naive, newdata = X_new) - 1   
  pred_pseudo_data_naive[["Tree"]] <- stats::predict(tree_naive, newdata = pseudo_test[, covariates_name]) - 1   
  
  hybrid_tree_naive <- policytree::hybrid_policy_tree(X_b, Gamma = DR.scores_b)
  pred_new_data_naive[["Hybrid_Tree"]] <- stats::predict(hybrid_tree_naive, newdata = X_new) - 1 
  pred_pseudo_data_naive[["Hybrid_Tree"]] <- stats::predict(hybrid_tree_naive, newdata = pseudo_test[, covariates_name]) - 1
  
  ## ── 4. Q-learning: SuperLearner ────────────────────────────────────────────
  idx_b1 <- train_b[[treatment_name]] == 1
  idx_b0 <- train_b[[treatment_name]] == 0
  
  SL_A1_b <- SuperLearner::SuperLearner(
    Y = as.numeric(Y_b[idx_b1]),
    X = train_b[idx_b1, covariates_name, drop = FALSE],
    SL.library = SL.library_cond, family = binomial())
  
  SL_A0_b <- SuperLearner::SuperLearner(
    Y = as.numeric(Y_b[idx_b0]),
    X = train_b[idx_b0, covariates_name, drop = FALSE],
    SL.library = SL.library_cond, family = binomial())
  
  df_new_df <- as.data.frame(df_new_sample[, covariates_name, drop = FALSE])
  new_preds1_naive <- predict(SL_A1_b, newdata = df_new_df)$pred
  new_preds0_naive <- predict(SL_A0_b, newdata = df_new_df)$pred
  pred_new_data_naive[["ql.SL"]] <- as.numeric(new_preds1_naive > new_preds0_naive)
  
  pseudo_preds1_naive <- predict(SL_A1_b, newdata = pseudo_df)$pred
  pseudo_preds0_naive <- predict(SL_A0_b, newdata = pseudo_df)$pred
  pred_pseudo_data_naive[["ql.SL"]] <- as.numeric(pseudo_preds1_naive > pseudo_preds0_naive)
  
  ## ── 6. Q-learning : generalized linear model with interactions ─────────────
  model.glb.glm <- stats::glm(formula = form_str, data = train_b, family = "binomial")
  
  po.glm <- sapply(numeric_levels_A, function(val) {
    df_new_template[[treatment_name]] <- val
    stats::predict(model.glb.glm, newdata = df_new_template, type = "response")
  })
  pred_new_data_naive[["ql.glm.interact"]] <- max.col(po.glm) - 1
  
  po_pseudo.glm <- sapply(numeric_levels_A, function(val) {
    pseudo_template[[treatment_name]] <- val
    stats::predict(model.glb.glm, newdata = pseudo_template, type = "response")
  })
  pred_pseudo_data_naive[["ql.glm.interact"]] <- max.col(po_pseudo.glm) - 1
  
  doptFactorPredict_new_naive <- do.call(cbind, pred_new_data_naive)
  doptFactorPredict_pseudo_naive <- do.call(cbind, pred_pseudo_data_naive)
  
  unweighted_probs_naive <- weighted_probs_experts(fitted_experts = doptFactorPredict_new_naive,
                                                   weights =rep(1/numalgs, numalgs),
                                                   df_pred = df_new_sample,
                                                   levels =levels_A)
  
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
  
  colnames(unweighted_probs_naive)  <- colnames(unweighted_probs_naive_pseudo)  <- as.numeric(levels_A)
  return(list(
    doptFactorPredict_cal = doptFactorPredict_cal, 
    doptFactorPredict_pseudo = doptFactorPredict_pseudo,
    doptFactorPredict_new_naive = doptFactorPredict_new_naive, 
    doptFactorPredict_pseudo_naive = doptFactorPredict_pseudo_naive,
    model.glb.pf = model.glb.pf, 
    model.glb.glm = model.glb.glm, 
    unweighted.naive_new = unweighted.naive,
    unweighted.pseudo.naive=unweighted.pseudo.naive))
}