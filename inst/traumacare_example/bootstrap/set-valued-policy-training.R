results_list <- parallel::mclapply(seq_len(n_bootstrap), function(bootstrap_idx) {
  system2("echo", args = sprintf("'iteration %d'", bootstrap_idx), stderr = "")
  # ── Bootstrap sample setup ──────────────────────────────────────────────────
  idx_b <- bootstrap_indices[[bootstrap_idx]]

  train_b <- df_obs[idx_b, ]
  set.seed(seed + bootstrap_idx)
  
  # ── Split data into 3 folds  ────────────────────────────────────────────────
  # Fold 1: Set-valued policy learning
  # Fold 2: Predict set-valued policies (r selection)
  # Fold 3: Train nuisances to evaluate on Fold2
  shuffled_indices <- sample(1:nrow(train_b))
  cut1 <- round(0.5 * nrow(train_b))
  cut2 <- round(0.2 * nrow(train_b))
  custom_folds <- list(
    Fold1 = shuffled_indices[1:cut1],
    Fold2 = shuffled_indices[(cut1 + 1): (cut1 + cut2)],
    Fold3 = shuffled_indices[(cut1 + cut2 + 1):nrow(train_b)])
  
  folds <- SuperLearner::CVFolds(nrow(train_b), id = NULL, Y = train_b[,outcome_name] ,
                                        cvControl = SuperLearner::SuperLearner.CV.control(V = 3L,
                                                                                          validRows = custom_folds))
  
  train <- train_b[folds[[1]],] # train set-valued policies 
  n_train <- nrow(train)
  pseudo.test.predict <-  train_b[folds[[2]],] # predict set-valued policies 
  evaluation <-  train_b[folds[[3]],] # train nuisances for set-valued policy evaluation
  
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
  
  
  # ── CONFORMAL SET_VALUED POLICY LEARNING ────────────────────────────────────
  # ── 0) Divide data into three even sets ─────────────────────────────────────
  # ── Noisy label generation, scoring model & calibration ─────────────────────
  set.seed(seed + bootstrap_idx)
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
  trained_policies_results <- train_policies(train_b = train_b, train1 = train1, 
                 calibration = calibration, 
                 seed = seed + bootstrap_idx) # GLB trained inside
  
  # Extract results from training
  doptFactorPredict_cal <- trained_policies_results$doptFactorPredict_cal # predictions on calibration
  doptFactorPredict_pseudo <- trained_policies_results$doptFactorPredict_pseudo
  doptFactorPredict_new_naive <- trained_policies_results$doptFactorPredict_new_naive
  results.single.policy <- trained_policies_results$results.policy 
  results.policy.agg <- trained_policies_results$results.policy.agg
  glb.model.grf <- trained_policies_results$glb.model.grf # GLB model (grf)
  glb.model.lm <- trained_policies_results$glb.model.lm # GLB models (lm)
  selected_methods <- trained_policies_results$selected_methods  # baselines used to train set-valued policies
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
  SL.library_cond <- c("SL.randomForest", "SL.ksvm", "SL.mean", "SL.glm", "SL.xgboost")
  QAW.reg.train <- SuperLearner::SuperLearner(
    Y = train2[, outcome_name], X = train2[, c(covariates_name, treatment_name)],
    SL.library = SL.library_cond, family = "gaussian")
  
  # ── 3) Calibration step  ────────────────────────────────────────────────────
  potential_outcomes_cal <- do.call(cbind,lapply(1:m, function(val) {
      new_data <- calibration[, c(covariates_name, treatment_name)]
      new_data[,treatment_name] <- factor(val, levels=levels_A)
      SuperLearner::predict.SuperLearner(QAW.reg.train, newdata = new_data)$pred}))
  
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
  
  # Compute margin score on new data (for final prediction)
  margin_po_pseudo <-  margin_score(potential_outcomes_pseudo)
  
  potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- df_new_sample[, c(covariates_name, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(QAW.reg.train, newdata = new_data)$pred}))
  
  # Compute margin score on new data (for final prediction)
  margin_po_new <-  margin_score(potential_outcomes_new)
  
  # Build and evaluate conformal set-valued policies   ─────────────────────────
  list.set.valued.policies <- list()
  # Baseline policy-based noisy labels
  conf_policy <- apply(r0_scores_policy, 2, function(x){
    quant <- stats::quantile(x, (1-alpha))
    binary_to_confidence_set(margin_po_pseudo < quant)
    }) # evaluation
  
  list.set.valued.policies[["conf.r"]] <- conf_policy
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
  
  # ── GREATEST LOWER BOUND (GLB) ──────────────────────────────────────────────
  # ── Using regression forest for estimation ──────────────────────────────────
  lowers <- uppers <- matrix(0, nrow=nrow(df_new_sample), ncol=m)
  for (l in as.numeric(levels_A)){
    data_l <- data.frame(df_new_sample[,covariates_name], Treatment=l)
    pred <- stats::predict(glb.model.grf, newdata = data_l, estimate.variance = TRUE)
    se <- sqrt(pred$variance.estimates)
    lowers[,l] <- pred$predictions - z * se
    uppers[,l] <- pred$predictions + z * se
  }
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
    uppers[,l] <- (pred$fit + z * se) %>% as.numeric()
  }
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
  indices <- c(cardinality = 1, spv_uniform.aipw = 2, spv_propensity.aipw = 3, 
               spv_uniform.tmle = 4, spv_propensity.tmle = 5)
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
      margin_po_new = margin_po_new, list.set.valued.policies))}, mc.cores = 4)  