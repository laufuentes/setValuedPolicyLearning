train_policies <- function(train_b, train1, calibration, seed) {
  set.seed(seed)
  
  message("training experts (conformal prediction)...\n")
  
  get_best_action <- function(model, base_df, m_levels, covariates, treatment, pred_fun = stats::predict) {
    n <- nrow(base_df)
    expanded_df <- base_df[rep(seq_len(n), times = length(m_levels)), c(covariates, treatment), drop = FALSE]
    expanded_df[[treatment]] <- factor(rep(m_levels, each = n), levels = levels_A)
    preds <- pred_fun(model, expanded_df)
    if (is.list(preds) && "pred" %in% names(preds)) preds <- preds$pred
    if (is.list(preds) && "predictions" %in% names(preds)) preds <- preds$predictions
    pred_mat <- matrix(preds, nrow = n, ncol = length(m_levels))
    max.col(pred_mat, ties.method = "first")
  }
  
  # Extract datasets
  # For conformal procedure (subset of complete boostrap sample)
  X_train <- train1[, covariates_name, drop = FALSE]
  A_train <- train1[, treatment_name]   
  Y_train <- train1[, outcome_name]   
  
  # For standard policy learning & GLB (complete boostrap sample)
  X_b <- train_b[, covariates_name, drop = FALSE]
  A_b <- train_b[, treatment_name]
  Y_b <- train_b[, outcome_name]
  
  m_levels <- 1:m
  df_new <- SL.out$df_new_sample
  
  pred_calibration <- list()

  ## ── 1. Q-learning: lm 
  ## ── 1.1 Q-learning: without interactions ───────────────────────────────────
  f_lm_simple <- stats::as.formula(paste(outcome_name, "~", paste(c(covariates_name,treatment_name), collapse = "+")))
  ql.lm_simple <- stats::lm(formula = f_lm_simple, data = train1)
  
  pred_calibration[["ql.lm"]] <- get_best_action(ql.lm_simple, calibration, m_levels, covariates_name, treatment_name)
  
  ## ── 1.2 Q-learning: with interactions ──────────────────────────────────────
  f_lm <- stats::as.formula(paste(outcome_name, "~ (", paste(covariates_name, collapse = "+"), ")*", treatment_name))
  ql.lm <- stats::lm(formula = f_lm, data = train1)
  
  pred_calibration[["ql.lm.interact"]] <- get_best_action(ql.lm, calibration, m_levels, covariates_name, treatment_name)

  ## ── 2. Q-learning: regression forest ───────────────────────────────────────
  ql.reg.forest <- grf::regression_forest(
    X = cbind(X_train, as.numeric(A_train)), 
    Y = Y_train, 
    seed = seed
  )
  
  rf_pred <- function(mod, data) {
    X_mat <- cbind(as.matrix(data[, covariates_name]), as.numeric(data[[treatment_name]]))
    stats::predict(mod, newdata = X_mat)$predictions
  }
  
  pred_calibration[["ql.reg.forest"]] <- get_best_action(ql.reg.forest, calibration, m_levels, covariates_name, treatment_name, pred_fun = rf_pred)

  ## ── 3. MACF ────────────────────────────────────────────────────────────────
  forest <- grf::multi_arm_causal_forest(
    X = X_train, Y = Y_train, W = as.factor(A_train)
  )
  
  cate_cal <- predict(forest, calibration[, covariates_name])$predictions[, , 1]
  pred_calibration[["MACF"]] <- max.col(cbind(0, cate_cal), ties.method = "first")
  
  ## ── 4. Policytree & Hybrid ──────────────────────────────────────────────────
  DR.scores <- policytree::double_robust_scores(forest)
  
  tree <- policytree::policy_tree(X_train, Gamma = DR.scores)
  pred_calibration[["Tree"]] <- stats::predict(tree, newdata = calibration[, covariates_name])

  hybrid_tree <- policytree::hybrid_policy_tree(X_train, Gamma = DR.scores, depth = 3)
  pred_calibration[["Hybrid_Tree"]] <- stats::predict(hybrid_tree, newdata = calibration[, covariates_name])

  ## ── 5. Q-learning: SuperLearner ────────────────────────────────────────────
  SL.library_cond <- c("SL.randomForest", "SL.ksvm", "SL.mean", "SL.glm", "SL.xgboost")
  
  QL_mod <- SuperLearner::SuperLearner( 
    Y = Y_train, X = train1[, c(covariates_name, treatment_name)],
    SL.library = SL.library_cond, family = "gaussian"
  )
  
  sl_pred <- function(mod, data) SuperLearner::predict.SuperLearner(mod, newdata = data)$pred
  
  pred_calibration[["ql.SL"]] <- get_best_action(QL_mod, calibration, m_levels, covariates_name, treatment_name, pred_fun = sl_pred)

  ## ── 6. POLLE Package ────────────────────────────────────────────────────────
  # Pre-convert datasets once
  prep_polle <- function(df) {
    dt <- data.table::as.data.table(df)
    dt[[treatment_name]] <- as.factor(dt[[treatment_name]])
    polle::policy_data(dt, action = treatment_name, covariates = covariates_name, utility = outcome_name)
  }
  
  pd_train <- prep_polle(train1)
  pd_cal   <- prep_polle(calibration)
  pd_new   <- prep_polle(df_new)
  
  formula_obj <- stats::reformulate(covariates_name)
  
  polle_configs <- list(
    drql.lm      = list(qv = polle::q_glm(formula = formula_obj, family = stats::gaussian()), q = q_glm()),
    drql.rf      = list(qv = polle::q_rf(formula = formula_obj), q = q_rf()),
    drql.xgboost = list(qv = polle::q_sl(formula = formula_obj, SL.library = "SL.xgboost"), q = q_sl()),
    drql.ksvm    = list(qv = polle::q_sl(formula = formula_obj, SL.library = "SL.ksvm"), q = q_sl())
  )
  
  for (name in names(polle_configs)) {
    cfg <- polle_configs[[name]]
    l_obj <- polle::policy_learn(type = "drql", control = polle::control_drql(qv_models = cfg$qv), cross_fit_g_models = FALSE)
    po <- l_obj(policy_data = pd_train, q_models = cfg$q, g_models = g_rf())
    
    pred_calibration[[name]] <- as.numeric(polle::get_policy(po)(pd_cal)$d)
  }
  
  SL.out$libraryNames <- names(pred_calibration)
  doptFactorPredict_cal <- do.call(cbind, pred_calibration)

  ## ── 7. Standard policy learning & GLB ──────────────────────────────────────
  message("training experts (baseline policies & GLB)...\n")
  
  pred_new_data_naive <- list()
  
  ql.lm_simple <- stats::lm(formula = f_lm_simple, data = train_b)
  
  pred_new_data_naive[["ql.lm"]] <- get_best_action(ql.lm_simple, df_new, m_levels, covariates_name, treatment_name)
  
  
  glb.model.lm <- stats::lm(formula = f_lm, data = train_b)
  pred_new_data_naive[["ql.lm.interact"]] <- get_best_action(glb.model.lm, df_new, m_levels, covariates_name, treatment_name)
  
  glb.model.grf <- grf::regression_forest(X = cbind(X_b, as.numeric(A_b)), Y = Y_b, seed = seed)
  pred_new_data_naive[["ql.reg.forest"]] <- get_best_action(glb.model.grf, df_new, m_levels, covariates_name, treatment_name, pred_fun = rf_pred)
  
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
    SL.library = SL.library_cond, family = "gaussian"
  )
  pred_new_data_naive[["ql.SL"]] <- get_best_action(QL_mod_b, df_new, m_levels, covariates_name, treatment_name, pred_fun = sl_pred)
  
  pd_train_b <- prep_polle(train_b)
  
  for (name in names(polle_configs)) {
    cfg <- polle_configs[[name]]
    l_obj <- polle::policy_learn(type = "drql", control = polle::control_drql(qv_models = cfg$qv), cross_fit_g_models = FALSE)
    po <- l_obj(policy_data = pd_train_b, q_models = cfg$q, g_models = g_rf())
    pred_new_data_naive[[name]] <- as.numeric(polle::get_policy(po)(pd_new)$d)
  }
  
  doptFactorPredict_new_naive <- do.call(cbind, pred_new_data_naive)
  numalgs_naive <- ncol(doptFactorPredict_new_naive)
  
  selected_methods <- switch(type,
                             "tree" = c("MACF", "ql.reg.forest", "drql.ksvm", "ql.lm"),
                             "linear" = c("MACF", "drql.lm", "drql.ksvm", "ql.lm"), 
                             "complex" = c("MACF", "drql.lm", "drql.ksvm", "ql.lm"))
  
  results.policy <- lapply(selected_methods, function(method){
    single.naive <- doptFactorPredict_new_naive[, method]
    table.evaluation(single.naive, 
                     optimal_policy_new = SL.out$optimal_policy_new, 
                     prop_score_new = SL.out$prop_score_new, 
                     potential_outcomes = SL.out$potential_outcomes, 
                     df_new_sample = SL.out$df_new_sample,
                     levels_A = levels_A,
                     treatment_name = treatment_name, 
                     outcome_name = outcome_name)})
  
  
  unweighted_probs_naive <- weighted_probs_experts(
    fitted_experts = doptFactorPredict_new_naive,
    weights = rep(1 / numalgs_naive, numalgs_naive),
    df_pred = df_new,
    levels = as.numeric(levels_A)
  )
  
  # Vectorized Multinomial Sampling replace apply loop
  unweighted_new_naive <- max.col(
    t(apply(unweighted_probs_naive, 1, function(p) stats::rmultinom(1, 1, prob = p))), 
    ties.method = "first"
  )
  results.policy.agg <- table.evaluation(unweighted_new_naive, 
                                         optimal_policy_new = SL.out$optimal_policy_new, 
                                         prop_score_new = SL.out$prop_score_new, 
                                         potential_outcomes = SL.out$potential_outcomes, 
                                         df_new_sample = SL.out$df_new_sample,
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
    results.policy.agg = results.policy.agg
  ))
}