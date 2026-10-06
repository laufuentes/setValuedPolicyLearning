results_list <- parallel::mclapply(seq_len(n_bootstrap), function(bootstrap_idx) {
  system2("echo", args = sprintf("'iteration %d'", bootstrap_idx), stderr = "")
  
  # ── Bootstrap sample setup ──────────────────────────────────────────────────
  idx_b <- bootstrap_indices[[bootstrap_idx]]
  train_b <- SL.out$df_obs[idx_b, ]
  X_b <- X[idx_b, ]
  A_b <- A[idx_b]
  Y_b <- Y[idx_b]
  
  optimal_policy_train <- SL.out$optimal_policy[idx_b]
  train_complete <- df_complete[idx_b, ]
  
  set.seed(seed + bootstrap_idx)
  
  # ── CONFORMAL SET_VALUED POLICY LEARNING ────────────────────────────────────
  # ── 0) Divide data into three even sets ─────────────────────────────────────
  # ── Noisy label generation, scoring model & calibration ─────────────────────
  folds <- SuperLearner::CVFolds(nrow(train_b), id = NULL, Y = Y_b,
                                 cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds,
                                                                                          shuffle = TRUE))
  train1 <- train_b[c(folds[[1]],folds[[4]]),] # generate noisy labels
  train2 <-  train_b[folds[[2]],] # score model and nuisances
  calibration <-  train_b[folds[[3]],] # calibration
  
  optimal_policy_cal <- optimal_policy_train[folds[[3]]]
  true_potential_outcomes_cal <-train_complete[folds[[3]],] |> 
    select(starts_with("Potential_outcomes."))
  
  # ── 1) Black-box label generation (i.e. estimates of (X,A*)) ───────────────────────
  ## 1.1) Generate random labels (i.e. A_rd)
  A_rd <- sample(as.numeric(levels_A), size = nrow(calibration), replace = TRUE)
  
  ## 1.2) Sample an optimal treatment per observation to generate oracular (X,A*)
  true_cal<- sapply(optimal_policy_cal, function(el) {
    if (length(el) == 1) el else el[sample(length(el), 1)]})
  
  ## 1.3) Estimate A* (OTR) using experts
  # Training performed on train1
  trained_policies_results <- train_policies(train_b = train_b, train1 = train1, 
                 calibration = calibration, 
                 seed = seed + bootstrap_idx) # GLB trained inside
  # Extract results from training
  doptFactorPredict_cal <- trained_policies_results$doptFactorPredict_cal # predictions on calibration 
  doptFactorPredict_new_naive <- trained_policies_results$doptFactorPredict_new_naive
  glb.model.grf <- trained_policies_results$glb.model.grf # GLB model (grf)
  glb.model.lm <- trained_policies_results$glb.model.lm # GLB models (lm)
  selected_methods <- trained_policies_results$selected_methods # baselines used to train set-valued policies
  results.single.policy <- trained_policies_results$results.policy 
  results.policy.agg <- trained_policies_results$results.policy.agg
  numalgs <- ncol(trained_policies_results$doptFactorPredict_cal) # number of baselines
  
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
  r <- ifelse(type=="tree", 0.15, 0.1)
  mix_factor<- stats::rbinom(nrow(calibration), 1, prob=r) # R ~ Ber(r)
  
  # Create perturbed labels
  noisy_policy <- mix_factor * A_rd + (1 - mix_factor) * policy_cal[, 1:2]
  policy_cal_r <- cbind(policy_cal[,1], noisy_policy[,1], 
                        policy_cal[,2], noisy_policy[,2],
                        policy_cal[, 3:4])
  
  selected_methods <- c(selected_methods[1],
                        paste0(selected_methods[1], paste0(" (r=",r,")")), 
                        selected_methods[2],
                        paste0(selected_methods[2],paste0(" (r=",r,")")), 
                        selected_methods[3:4]) 
  
  # ── 2) Train nonconformity score model (i.e. s(X,A)) ───────────────────────
  # Training performed on train2
  # Two predictions:
  # (i) on calibration
  # (ii) on SL.out$df_new_sample
  SL.library_cond <- c("SL.randomForest", "SL.ksvm", "SL.mean", "SL.glm", "SL.xgboost")
  QAW.reg.train <- SuperLearner::SuperLearner(
    Y = train2[, outcome_name], X = train2[, c(covariates_name, treatment_name)],
    SL.library = SL.library_cond, family = "gaussian")
  
  # ── 3) Calibration step  ────────────────────────────────────────────────────
  potential_outcomes_cal <- do.call(cbind,lapply(1:m, function(val) {
      new_data <- calibration[, c(covariates_name, treatment_name)]
      new_data[,treatment_name] <- factor(val, levels=levels_A)
      SuperLearner::predict.SuperLearner(QAW.reg.train, newdata = new_data)$pred}))
  
  # Compute margin scores on calibration data 
  margin_po <-  margin_score(potential_outcomes_cal)
  
  # Extract scores for different label types
  r0_scores_policy <- apply(policy_cal_r,2,function(x){
    margin_po[cbind(1:nrow(calibration), x)]})  # baseline policies and perturbed labels
  r0_scores_aggregation <- margin_po[cbind(1:nrow(calibration), 
                                           unweighted_aggregation)] # aggregation of baseline policies 
  r1_score <- margin_po[cbind(seq_len(nrow(calibration)), A_rd)] # random labels 
  r_true <- margin_po[cbind(seq_len(nrow(calibration)), true_cal)] # true labels
  
  potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- SL.out$df_new_sample[, c(covariates_name, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(QAW.reg.train, newdata = new_data)$pred}))
  
  # Compute margin score on new data (for final prediction)
  margin_po_new <-  margin_score(potential_outcomes_new)
  
  # Build and evaluate conformal set-valued policies   ─────────────────────────
  # Baseline policy-based noisy labels
  results_policy <- apply(r0_scores_policy,2, function(x){
    quant <- stats::quantile(x, (1-alpha))
    conf_set_policy <- binary_to_confidence_set(margin_po_new < quant)
    table.evaluation(conf_set_policy, 
                     optimal_policy_new = SL.out$optimal_policy_new, 
                     prop_score_new = SL.out$prop_score_new, 
                     potential_outcomes = SL.out$potential_outcomes, 
                     df_new_sample = SL.out$df_new_sample,
                     levels_A = levels_A, 
                     treatment_name = treatment_name, 
                     outcome_name = outcome_name)}) # evaluation 
  
  # Aggregation-based noisy labels
  quant_agg <- quantile(r0_scores_aggregation, 1 - alpha)
  conf_set_agg <- binary_to_confidence_set(margin_po_new < quant_agg)
  results_agg <- table.evaluation(conf_set_agg, 
                                  optimal_policy_new = SL.out$optimal_policy_new, 
                                  prop_score_new = SL.out$prop_score_new, 
                                  potential_outcomes = SL.out$potential_outcomes, 
                                  df_new_sample = SL.out$df_new_sample,
                                  levels_A = levels_A,
                                  treatment_name = treatment_name, 
                                  outcome_name = outcome_name) # evaluation 
  
  # ── GREATEST LOWER BOUND (GLB) ──────────────────────────────────────────────
  # ── Using regression forest for estimation ──────────────────────────────────
  lowers <- uppers <- matrix(0, nrow=nrow(SL.out$df_new), ncol=m)
  for (l in as.numeric(levels_A)){
    data_l <- data.frame(SL.out$df_new[,covariates_name], Treatment=l)
    pred <- stats::predict(glb.model.grf, newdata = data_l, estimate.variance = TRUE)
    se <- sqrt(pred$variance.estimates)
    lowers[,l] <- pred$predictions - z * se
    uppers[,l] <- pred$predictions + z * se
  }
  uppest_lrw_bound <- apply(lowers, 1, max)
  conf_set_grf <- binary_to_confidence_set(uppers >= uppest_lrw_bound)
  results_glb_grf <- table.evaluation(conf_set_grf, 
                                      optimal_policy_new = SL.out$optimal_policy_new, 
                                      prop_score_new = SL.out$prop_score_new, 
                                      potential_outcomes = SL.out$potential_outcomes, 
                                      df_new_sample = SL.out$df_new_sample,
                                      levels_A = levels_A,
                                      treatment_name = treatment_name, 
                                      outcome_name = outcome_name) # evaluation
  
  # ── Using linear model with interactions for estimation ─────────────────────
  lowers <- uppers <- matrix(0, nrow=nrow(SL.out$df_new), ncol=m)
  for (l in as.numeric(levels_A)){
    data_l <- data.frame(SL.out$df_new[,covariates_name], A=factor(l, levels = levels_A))
    pred <- stats::predict(glb.model.lm, newdata = data_l, se.fit = TRUE)
    se <- pred$se.fit
    lowers[,l] <- (pred$fit - z * se) |> as.numeric()
    uppers[,l] <- (pred$fit + z * se) |> as.numeric()
  }
  uppest_lrw_bound <- apply(lowers, 1, max)
  conf_set_lm <- binary_to_confidence_set(uppers>=uppest_lrw_bound)
  results_glb_lm <- table.evaluation(conf_set_lm, 
                                     optimal_policy_new = SL.out$optimal_policy_new, 
                                     prop_score_new = SL.out$prop_score_new, 
                                     potential_outcomes = SL.out$potential_outcomes, 
                                     df_new_sample = SL.out$df_new_sample,
                                     levels_A = levels_A,
                                     treatment_name = treatment_name, 
                                     outcome_name = outcome_name) # evaluation
  
  # ── Collect results for this bootstrap iteration ────────────────────────────
  indices <- c(exact_match = 1, coverage = 3, 
               relaxed_coverage = 4, cardinality = 2, 
               spv_uniform = 5, spv_propensity = 6)
  
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
    list(
      selected_methods = selected_methods, 
      doptFactorPredict_new_naive = doptFactorPredict_new_naive,
      r0_scores_policy = r0_scores_policy,
      r0_scores_agg = r0_scores_aggregation,
      margin_po_new = margin_po_new))}, mc.cores = 4)  