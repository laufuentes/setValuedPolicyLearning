#' Learn set valued policies
#'
#' Generates set-valued policies using Greatest Lower Bound and conformal 
#' approach based on noisy labels. 
#'
#' @param X A numeric matrix of covariates. Categorical variables must be 
#'   pre-processed (e.g., one-hot encoded).
#' @param A A factor vector representing the observed treatment assignments.
#' @param Y A numeric vector or matrix of primary outcomes.
#' @param X_test A numeric matrix of covariates for the test/evaluation set.
#' @param A_test A factor vector of treatment assignments for the test set.
#' @param Y_test A numeric vector or matrix of outcomes for the test set.
#' @param random_rate A numeric vector of randomness levels for the conformal 
#'   approach. Defaults to `c(0, 0.2, 0.5)`.
#' @param alphas A numeric vector of significance levels (1 - confidence) used 
#'   to construct the set-valued policies.
#' @param VFolds Integer specifying the number of cross-validation folds. 
#'   Currently must be set to `3L`.
#' @param label_generation A list of configurations for noisy label generation 
#'   (via the \code{polle} package). Each list element should contain:
#'   \itemize{
#'     \item \code{type}: Type of learner: \code{"drql"} (Doubly Robust), 
#'       \code{"ql"} (Q-learning), or \code{"plt"} (Policy Tree).
#'     \item \code{q_func}: ML model for the Q-function: \code{"q_glm"}, 
#'       \code{"q_rf"}, \code{"q_sl"}, \code{"q_xgboost"}, or \code{"q_glmnet"}.
#'     \item \code{sl_library}: Character vector of learners (only required 
#'       if \code{q_func = "q_sl"}).
#'   }
#'   Defaults to a list of DR-QL with GLM and SuperLearner (XGBoost).
#' @param SL.library.nuisance Character vector of SuperLearner libraries used 
#'   to estimate nonconformity scores (nuisance functions). 
#'   Defaults to \code{c("SL.mean", "SL.glm", "SL.ranger", "SL.xgboost")}.
#'  
#' @return A list with the following components:
#' \item{SL.out}{A list containing the features and models used during training.}
#' \item{confidence_sets}{A nested list containing the resulting set-valued 
#'   policies: \code{GLB} for the Greatest Lower Bound approach, and 
#'   \code{conformal} for results across different \code{random_rate} values.}
#'
#' @export
learn_set_valued_policies <- function(X, A, Y, X_test, A_test, Y_test,
                                      random_rate=c(0,0.2,0.5),
                                      alphas = seq(0,0.5,0.05), VFolds=3,
                                      label_generation = list(
                                        list(type = "drql", q_func = "q_glm"), 
                                        list(type = "drql", q_func = "q_sl", sl_library = "SL.xgboost")),
                                      SL.library.nuisance = c("SL.mean", "SL.glm","SL.ranger", "SL.xgboost")){

  # ── Data sanity checks ──────────────────────────────────────────────────────
  # Check no missing data
  if (any(is.na(X)) || any(is.na(Y)) || 
      any(is.na(X_test)) || 
      any(is.na(A)) || any(is.na(A_test))) {
    stop(paste("Execution Halted: Missing values (NAs) detected in imput data"))
  }
  # Ensure covariates are numeric
  is_X_numeric <- all(sapply(X, is.numeric)) && all(sapply(X_test, is.numeric))
  if (!is_X_numeric) {
    stop("Covariates (X/X_test) must be numeric. Please check your one-hot encoding.")
  }
  covariate_names <- colnames(X)
  # Verify treatments are factors
  is_A_factor <- is.factor(A) && is.factor(A_test)
  if (!is_A_factor) {
    stop("Treatment variables (A/A_test) must be factors.")
  }
  # Verify treatments have same levels
  if (!identical(levels(A), levels(A_test))) {
    stop("A and A_test have mismatched factor levels.")
  }
  levels_A <- levels(A) # treatment levels
  m <- length(levels(A)) # number of treatment levels
  treatment_name <- "A"
  # Outcome family
  family = ifelse(max(Y) <= 1 & min(Y) >= 0, "binomial", "gaussian") # in [0,1] or beyond
  ab <- c(min(cbind(Y,Y_test)),max(cbind(Y,Y_test)))
  outcome_name <- "Y"
  
  
  # Check that
  alphas_form <- all(alphas>=0 & alphas<=1)
  if (!alphas_form) {
    stop("Confidence levels beyond interval [0,1].")
  }
  
  randomness_form <- all(random_rate>=0 & random_rate<=1)
  if (!randomness_form) {
    stop("Randomness levels beyond interval [0,1].")
  }
  n_rate <- length(random_rate) # number of random rates to test
  
  
  # ── Storage object  ─────────────────────────────────────────────────────────
  SL.out<- list()
  SL.out$feature_names <- list(covariate_names, treatment_name, outcome_name)
  SL.out$ab <- ab
  SL.out$family = family
  SL.out$df_obs <- data.frame(X, A, Y)
  colnames(SL.out$df_obs) <- c(covariate_names, treatment_name, outcome_name)
  n <- nrow(SL.out$df_obs)

  SL.out$df_new_sample <- data.frame(X_test, A_test, Y_test)
  colnames(SL.out$df_new_sample) <- c(covariate_names, treatment_name,
                                      outcome_name)

  # ── 0) Divide data into three even sets ─────────────────────────────────────
  SL.out$folds <- SuperLearner::CVFolds(n, id = NULL,Y = Y,
                                        cvControl = SuperLearner::SuperLearner.CV.control(V = VFolds,
                                                                                          shuffle = TRUE))

  train1 <- SL.out$df_obs[SL.out$folds[[1]],] # generate noisy labels
  train2 <-  SL.out$df_obs[SL.out$folds[[2]],] # score model and nuisances
  test <-  SL.out$df_obs[SL.out$folds[[3]],] # calibration

  # ── 1) Black-box label generation (i.e. estimates of (X,A*)) ────────────────
  ## 1.1) Generate random labels (i.e. A_rd)   ─────────────────────────────────
  A_rd <- apply(data.frame(1:nrow(test)),1,
                function(i)sample(as.numeric(levels_A),size=1))

  ## 1.2) Estimate A* (OTR) using experts  ─────────────────────────────────────
  # Training performed on train1
  # Two predictions:
  # (i) on test
  # (ii) on SL.out$df_new_sample
  # ──   ─────────────────────────────────────────────────────────
  q_learners <- list() 
  for (conf in label_generation) {
    args_list <- list(
      name = paste0(conf$type,"_", conf$q_func), type = conf$type, 
      q_func = conf$q_func, sl_library = conf$sl_library)
    
    q_learners <- q_learners %>% 
      add_qlearner(name = args_list$name, type = args_list$type, 
                   q_func = args_list$q_func, sl_library = args_list$sl_library,
                   action_name = treatment_name, covariates = covariate_names)
  }
  
  g_learner <- glearner(m, g_func = "g_rf", sl_library = NULL, num.trees = 500)
  
  SL.init1 = expert_fit_predict(train1, test, new = SL.out$df_new_sample,
                                covariates = covariate_names,
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
  # 2.1) Train nuisances   ─────────────────────────────────────────────────────
  ## Outcome model (Q-model)
  SL.out$QAW.reg.train = SuperLearner::SuperLearner(
    Y=train2[,outcome_name], X = train2[,c(covariate_names,treatment_name)],
    SL.library=SL.library.nuisance, family = SL.out$family)

  ## Propensity score  (G-model)
  SL.out$g.reg.train <- randomForest::randomForest(x = train2[,covariate_names],
                                                   y = train2[,treatment_name])

  # 2.2) Predict nonconformity scores (margin score) ───────────────────────────
  # Nonconformity scores on calibration data
  potential_outcomes_test <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- test[, c(covariate_names, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(SL.out$QAW.reg.train, newdata = new_data)$pred}))
  margin_po <-  margin_score(potential_outcomes_test) # score for all potential outcomes from test

  # Nonconformity scores on new data (used to generate sets)
  potential_outcomes_new <- do.call(cbind,lapply(1:m, function(val) {
    new_data <- SL.out$df_new[, c(covariate_names, treatment_name)]
    new_data[,treatment_name] <- factor(val, levels=levels_A)
    SuperLearner::predict.SuperLearner(SL.out$QAW.reg.train, newdata = new_data)$pred}))
  SL.out$new_scores <- margin_score(potential_outcomes_new)  # score for all potential outcomes from new data

  # ── 3) Noisy calibration  ───────────────────────────────────────────────────
  # 3.1) Generate perturbed labels  ────────────────────────────────────────────
  SL.out$rate_cal_labels_unweighted <- SL.out$rate_scores_unweighted_cal<- matrix(0,nrow=nrow(test), ncol=n_rate)
  SL.out$rate_scores_unweighted_true <- matrix(0,nrow=nrow(test), ncol=n_rate)
  for(i in 1:length(random_rate)){
    rate <- random_rate[i] # randomness level r
    mix_factor<- stats::rbinom(nrow(test),1,prob=rate) # R ~ Ber(r)
    # Combine noisy labels with random
    SL.out$rate_cal_labels_unweighted[,i]<- mix_factor*A_rd +  (1-mix_factor)*unweighted_cal
    # Compute associated scores
    SL.out$rate_scores_unweighted_cal[,i] <- margin_po[cbind(1:nrow(test), SL.out$rate_cal_labels_unweighted[,i])]
  }
  SL.out$rate_cal_labels_behavioral <- margin_po[cbind(1:nrow(test), test[,treatment_name])]
  
  # ─── Save labels  ───────────────────────────────────────────────────────────
  data_toghether <- dplyr::bind_rows(
    lapply(1:ncol(SL.out$rate_scores_unweighted_cal), function(i) {
      data.frame(
        value = SL.out$rate_scores_unweighted_cal[,i],
        model = "Estimated score",
        mechanism = "Unweighted",
        type = random_rate[i])}))
  
  SL.out$data_toghether <- data_toghether
  
  # 3.2) Generate set-valued policies  ─────────────────────────────────────────
  confidence_sets <- list()
  confidence_sets[["GLB"]] <- list()
  confidence_sets[["conformal"]] <- list()
  
  treatment<- matrix(0, nrow=nrow(rbind(train1, train2)), ncol=m)
  treatment[cbind(1:nrow(treatment), c(train1[, treatment_name],
                                       train2[, treatment_name]))] <- 1
  model <- grf::regression_forest(
    X = cbind(rbind(train1, train2)[, covariate_names], treatment),
    Y = rbind(train1, train2)[, outcome_name])
  
  for(i in 1:length(alphas)){
    alpha <- alphas[i]
    alpha_key <- paste0("alpha=",alpha)
    confidence_sets[["conformal"]][[alpha_key]] <- list()
    for (r in 1:ncol(SL.out$rate_cal_labels_unweighted)){
      r_key <- paste0("r=",random_rate[r])
      # conformal set-valued policy learning
      quant <- stats::quantile(SL.out$rate_scores_unweighted_cal[, r], (1-alpha))
      binary_confidence_set <-  ifelse(SL.out$new_scores<quant, 1, 0)
      idx <- which(binary_confidence_set  != 0, arr.ind = TRUE)
      confidence_sets[["conformal"]][[alpha_key]][[r_key]]  <- heatmap_treatments(
        split(idx[, "col"], factor(idx[, "row"],levels = seq_len(nrow(binary_confidence_set)))), levels_A) %>% as.matrix()
    }
    lowers <- uppers <- lowers_test <- uppers_test <- matrix(0,
                                                             nrow=nrow(SL.out$df_new),
                                                             ncol=m)
    z <- stats::qnorm(1 - alpha/2)
    for (l in as.numeric(levels_A)){
      treatment_l <- matrix(0, nrow=nrow(SL.out$df_new), ncol=m)
      treatment_l[,l] <- 1
      data_l <- data.frame(SL.out$df_new[,covariate_names], Treatment=treatment_l)
      pred <- stats::predict(model, newdata = data_l, estimate.variance = TRUE)
      se <- sqrt(pred$variance.estimates)
      lowers[,l] <- pred$predictions - z * se
      uppers[,l] <- pred$predictions + z * se
    }
    uppest_lrw_bound <- apply(lowers, 1, max)
    C_set_binary_naive <- ifelse(uppers>=uppest_lrw_bound, 1, 0)
    indices_naive <- which(C_set_binary_naive != 0, arr.ind = TRUE)
    confidence_sets[["GLB"]][[alpha_key]] <- heatmap_treatments(
      split(indices_naive[, "col"], indices_naive[, "row"]), levels_A) %>% as.matrix()
  }
  return(list(SL.out=SL.out, confidence_sets=confidence_sets))
}
