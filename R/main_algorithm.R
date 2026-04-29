learn_set_valued_policies <- function(X, A, Y, X_test, A_test, Y_test,
                                      random_rate=c(0,0.2,0.5),
                                      alphas = seq(0,0.5,0.05),
                                      SL.library.nuisance = c("SL.mean", "SL.glm",
                                        "SL.xgboost","SL.ranger","SL.ksvm"),
                                      root.path = ".",
                                      ){
  if (!dir.exists(root.path)) {
    warning(sprintf("The directory '%s' does not exist. Creating it...", root.path))
    dir.create(root.path, recursive = TRUE)
  }

  subdirs <- c("images", "predictions")

  for (subdir in subdirs) {
    subdir_path <- file.path(root.path,"inst", subdir)

    if (!dir.exists(subdir_path)) {
      message(sprintf("Creating subdirectory: %s", subdir))
      dir.create(subdir_path, recursive = TRUE)
    }
  }
  message("Directory check complete.")

  # Previous sanity checks:
  # A and A_test must be numeric

  SL.out<- list()

  family = ifelse(max(Y) <= 1 & min(Y) >= 0, "binomial", "gaussian") # in [0,1] or beyond
  SL.out$family = family
  ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new)))

  covariate_name <- NULL
  treatment_name <- NULL
  outcome_name <- NULL
  SL.out$df_obs <- data.frame(X, A, Y)

  SL.out$df_new_sample <- data.frame(X_test, A_test, Y_test)

  levels_A <- levels(A) # treatment levels
  m <- length(levels(A)) # number of treatment levels

  # ── 0) Divide data into three even sets ─────────────────────────────────────
  SL.out$folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y,
                                        cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds,
                                                                                          shuffle = TRUE))

  train1 <- SL.out$df_obs[SL.out$folds[[1]],] # generate noisy labels
  train2 <-  SL.out$df_obs[SL.out$folds[[2]],] # score model and nuisances
  test <-  SL.out$df_obs[SL.out$folds[[3]],] # calibration

  # ── 1) Black-box label generation (i.e. estimates of (X,A*)) ────────────────
  ## 1.1) Generate random labels (i.e. A_rd)
  A_rd <- apply(data.frame(1:nrow(test)),1,
                function(i)sample(as.numeric(levels_A),size=1))

  ## 1.2) Estimate A* (OTR) using experts
  # Training performed on train1
  # Two predictions:
  # (i) on test
  # (ii) on SL.out$df_new_sample
  SL.init1 = expert_fit_predict(train1, test, new = SL.out$df_new_sample,
                                covariates = covariates_name,
                                treatment_name=treatment_name,
                                outcome_name=outcome_name,
                                qlearners_list = q_learners, g_model=g_learner)

  SL.out$libraryNames <- SL.init1[["learner_names"]] # expert names
  numalgs <- length(SL.out$libraryNames) # number of experts

  # predictions on calibration
  SL.out$doptFactorPredict_test <- SL.init1[["expert_policies"]]
  # predictions on new data
  SL.out$doptFactorPredict_new <- SL.init1[["expert_policies_new"]]

  # Generate noisy calibration labels
  unweighted_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_test,
                                             weights =rep(1/numalgs, numalgs),
                                             df_pred = test,
                                             levels = as.numeric(levels_A))

  unweighted_cal <- apply(
    apply(unweighted_probs, 1, function(x){rmultinom(1, 1, prob=x)}),
    2, which.max)

  # ── 2) Nonconformity score model (i.e. s(X,A)) ───────────────────────
  # Training performed on train2
  # Two predictions:
  # (i) on test
  # (ii) on SL.out$df_new_sample
  # 2.1) Train nuisances
  ## Outcome model (Q-model)
  SL.out$QAW.reg.train = SuperLearner::SuperLearner(
    Y=train2[,outcome_name], X = train2[,c(covariates_name,treatment_name)],
    SL.library=SL.library.nuisance, family = SL.out$family)

  ## Propensity score  (G-model)
  SL.out$g.reg.train <- randomForest::randomForest(x = train2[,covariates_name],
                                                   y = train2[,treatment_name])

  # 2.2) Predict nonconformity scores (margin score)
  # Nonconformity scores on calibration data
  potential_outcomes_test <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- test[, c(covariates_name, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(SL.out$QAW.reg.train, newdata = new_data)$pred}))
  margin_po <-  margin_score(potential_outcomes_test) # score for all potential outcomes from test

  # Nonconformity scores on new data (used to generate sets)
  potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- SL.out$df_new[, c(covariates_name, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(SL.out$QAW.reg.train, newdata = new_data)$pred}))
  SL.out$new_scores <- margin_score(potential_outcomes_new)  # score for all potential outcomes from new data

  # ── Generate perturbed labels  ─────────────────────────────────────────────────
  SL.out$rate_cal_labels_unweighted <- SL.out$rate_scores_unweighted_cal<- matrix(0,nrow=nrow(test), ncol=n_rate)
  SL.out$rate_scores_unweighted_true <- matrix(0,nrow=nrow(test), ncol=n_rate)
  for(i in 1:random_rate){
    rate <- random_rate[i] # randomness level r
    mix_factor<- stats::rbinom(nrow(test),1,prob=rate) # R ~ Ber(r)
    # Combine noisy labels with random
    SL.out$rate_cal_labels_unweighted[,i]<- mix_factor*A_rd +  (1-mix_factor)*unweighted_cal
    # Compute associated scores
    SL.out$rate_scores_unweighted_cal[,i] <- margin_po[cbind(1:nrow(test), SL.out$rate_cal_labels_unweighted[,i])]
    if(synthetic_scenario){
      # Compute the true scores
      SL.out$rate_scores_unweighted_true[,i] <- margin_true_test[cbind(1:nrow(test), SL.out$rate_cal_labels_unweighted[,i])]
    }
  }

  data_toghether <- dplyr::bind_rows(
    lapply(1:ncol(SL.out$rate_scores_unweighted_cal), function(i) {
      data.frame(
        value = SL.out$rate_scores_unweighted_cal[,i],
        model = "Estimated score",
        mechanism = "Unweighted",
        type = random_rate[i])}))
  SL.out$rate_cal_labels_behavioral <- margin_po[cbind(1:nrow(test), test[,treatment_name])]

  for(i in 1:length(alphas)){
    alpha <- alphas[i]
    for (r in 1:ncol(SL.out$rate_cal_labels_unweighted)){
      # conformal set-valued policy learning
      quant <- stats::quantile(SL.out$rate_scores_unweighted_cal[, r], (1-alpha))
      binary_confidence_set <-  ifelse(SL.out$new_scores<quant, 1, 0)
      idx <- which(binary_confidence_set  != 0, arr.ind = TRUE)
      confidence_set <- split(idx[, "col"],
                              factor(idx[, "row"],
                                     levels = seq_len(nrow(binary_confidence_set))))
      }}

}
