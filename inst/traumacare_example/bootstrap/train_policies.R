train_policies <- function(train_b, train1, calibration, pseudo_test, seed) {
  set.seed(seed)
  
  message("training experts (conformal prediction)...\n")
    get_best_action <- function(model, base_df, m_levels, covariates, treatment, pred_fun = stats::predict) {
    n <- nrow(base_df)
    # Expand data once across all treatment levels
    expanded_df <- base_df[rep(seq_len(n), times = length(m_levels)), c(covariates, treatment), drop = FALSE]
    expanded_df[[treatment]] <- factor(rep(m_levels, each = n), levels = levels_A)
    
    preds <- pred_fun(model, expanded_df)
    if (is.list(preds) && "pred" %in% names(preds)) preds <- preds$pred
    if (is.list(preds) && "predictions" %in% names(preds)) preds <- preds$predictions
    pred_mat <- matrix(preds, nrow = n, ncol = length(m_levels))
    max.col(pred_mat, ties.method = "first")
  }
  
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
  
  ## ── 1. Probability forest
  proba.forest <- grf::probability_forest(X=cbind(X_train, A_train), 
                                     Y = Y_train %>% as.factor(), seed = seed)
  
  pf_pred <- function(mod, data) {
    X_mat <- cbind(as.matrix(data[, covariates_name, treatment_name]))
    X_mat[, treatment_name] <- as.numeric(data[[treatment_name]])
    stats::predict(mod, newdata = X_mat)$predictions[,2]
  }
  
  pred_calibration[["proba.forest"]] <- get_best_action(proba.forest, calibration, 
                                                        m_levels, covariates_name, 
                                                        treatment_name, 
                                                        pred_fun = pf_pred, Factor=FALSE)
  
  pred_pseudo_data[["proba.forest"]] <- get_best_action(proba.forest, pseudo_test, 
                                                        m_levels, covariates_name, treatment_name, pred_fun = pf_pred)

  ## ── 2. MACF ─────────────────────────────────────────────────────────────────
  forest <- grf::multi_arm_causal_forest(
    X = X_train, Y = Y_train, W = as.factor(A_train))
  
  cate_cal <- predict(forest, calibration[, covariates_name])$predictions[, , 1]
  pred_calibration[["MACF"]] <- max.col(cbind(0, cate_cal), ties.method = "first")
  
  cate_cal_p <- predict(forest, pseudo_test[, covariates_name])$predictions[, , 1]
  pred_pseudo_data[["MACF"]] <- max.col(cbind(0, cate_cal_p), ties.method = "first")
  
  ## ── 3. Policytree & Hybrid ──────────────────────────────────────────────────
  DR.scores <- policytree::double_robust_scores(forest)
  
  tree <- policytree::policy_tree(X_train, Gamma = DR.scores)
  pred_calibration[["Tree"]] <- stats::predict(tree, newdata = calibration[, covariates_name])
  pred_pseudo_data[["Tree"]] <- stats::predict(tree, newdata = pseudo_test[, covariates_name])

  hybrid_tree <- policytree::hybrid_policy_tree(X_train, Gamma = DR.scores, depth = 3)
  pred_calibration[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = calibration[, covariates_name])
  pred_pseudo_data[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = pseudo_test[, covariates_name])

  ## ── 5. Q-learning: SuperLearner ────────────────────────────────────────────
  SL.library_cond <- c("SL.randomForest", "SL.ksvm", "SL.mean", "SL.glm", "SL.xgboost")
  
  QL_mod <- SuperLearner::SuperLearner( 
    Y = Y_train, X = train1[, c(covariates_name, treatment_name)],
    SL.library = SL.library_cond, family = "binomial")
  
  sl_pred <- function(mod, data) SuperLearner::predict.SuperLearner(mod, newdata = data)$pred
  
  pred_calibration[["ql.SL"]] <- get_best_action(QL_mod, calibration, m_levels, covariates_name, treatment_name, pred_fun = sl_pred)
  pred_pseudo_data[["ql.SL"]] <- get_best_action(QL_mod, pseudo_test, m_levels, covariates_name, treatment_name, pred_fun = sl_pred)
  
  ## ── 6. Q-learning with linear model  ───────────────────────────────────────
  f_lm <- stats::as.formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name))
  ql.lm <- stats::glm(formula = f_lm, data = train1, family = "binomial")
  
  pred_calibration[["ql.lm.interact"]] <- get_best_action(ql.lm, calibration, m_levels, covariates_name, treatment_name)
  pred_pseudo_data[["ql.lm.interact"]] <- get_best_action(ql.lm, pseudo_test, m_levels, covariates_name, treatment_name)
  
  libraryNames <- names(pred_calibration)
  doptFactorPredict_cal <- do.call(cbind, pred_calibration)
  doptFactorPredict_pseudo <- do.call(cbind, pred_pseudo_data)

  ## ── 7. Naive method & GLB ──────────────────────────────────────────────────
  message("training experts (naive version & GLB)...\n")
  
  pred_new_data_naive <- list()
  
  glb.model.lm <- stats::glm(formula = f_lm, data = train_b, family = "binomial")
  pred_new_data_naive[["ql.lm.interact"]] <- get_best_action(glb.model.lm, df_new, m_levels, covariates_name, treatment_name)
  
  glb.model.pf <- grf::probability_forest(X = cbind(X_b, as.numeric(A_b)), Y = Y_b %>% as.factor(), seed = seed)
  pred_new_data_naive[["proba.forest"]] <- get_best_action(glb.model.pf, df_new, m_levels, covariates_name, treatment_name, pred_fun = pf_pred)
  
  forest_b <- grf::multi_arm_causal_forest(X = X_b, Y = Y_b, W = as.factor(A_b))
  cate_new_b <- predict(forest_b, X_new)$predictions[, , 1]
  pred_new_data_naive[["MACF"]] <- max.col(cbind(0, cate_new_b), ties.method = "first")
  
  DR.scores_b <- policytree::double_robust_scores(forest_b)
  tree_b <- policytree::policy_tree(X_b, Gamma = DR.scores_b)
  pred_new_data_naive[["Tree"]] <- stats::predict(tree_b, X_new)
  
  hybrid_tree_b <- policytree::hybrid_policy_tree(X_b, Gamma = DR.scores_b)
  pred_new_data_naive[["Hybrid_Tree"]] <- stats::predict(hybrid_tree_b, newdata = X_new)
  
  QL_mod_b <- SuperLearner::SuperLearner(
    Y = Y_b, X = train_b[, c(covariates_name, treatment_name)],
    SL.library = SL.library_cond, family = "binomial"
  )
  pred_new_data_naive[["ql.SL"]] <- get_best_action(QL_mod_b, df_new, m_levels, covariates_name, treatment_name, pred_fun = sl_pred)
  
  doptFactorPredict_new_naive <- do.call(cbind, pred_new_data_naive)
  numalgs_naive <- ncol(doptFactorPredict_new_naive)
  
  selected_methods <- c("MACF")
  
  
  unweighted_probs_naive <- weighted_probs_experts(
    fitted_experts = doptFactorPredict_new_naive,
    weights = rep(1 / numalgs_naive, numalgs_naive),
    df_pred = df_new,
    levels = as.numeric(levels_A))
  
  # Vectorized Multinomial Sampling replace apply loop
  unweighted_new_naive <- max.col(
    t(apply(unweighted_probs_naive, 1, function(p) stats::rmultinom(1, 1, prob = p))), 
    ties.method = "first")
  
  results.policy.agg <- table.evaluation.real(unweighted_new_naive, 
                                              prop_score_new = gAW.pred.pseudo.r, 
                                              potential_outcomes = Q.all.pseudo.r,
                                              df_new_sample = df_new_sample,
                                              levels_A = levels_A, 
                                              treatment_name = treatment_name, 
                                              outcome_name = outcome_name)
  
  return(list(
    doptFactorPredict_cal = doptFactorPredict_cal, 
    doptFactorPredict_new_naive = doptFactorPredict_new_naive, 
    glb.model.grf = glb.model.grf, 
    glb.model.lm = glb.model.lm, 
    selected_methods = selected_methods,
    results.policy = results.policy, 
    results.policy.agg = results.policy.agg))
}