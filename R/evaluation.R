#' Average relaxed coverage
#'
#' Computes the average proportion of the true set that is contained within
#' the predicted set. This is a normalized measure of recall across multiple
#' observations.
#'
#' @param true_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param pred_set A `list` of numeric or character vectors representing the predicted sets.
#'
#' @return A numeric value representing the mean coverage (proportion of
#'   intersected elements over true set size) across all observations.
#' @export
#'
#' @examples
#' true <- list(c(1, 2), c(3))
#' pred <- list(c(1), c(3, 4))
#' coverage_relaxed(true, pred)
coverage_relaxed <- function(true_set, pred_set){
  if(!(is.list(pred_set) & is.list(true_set))){
    msg_list <- paste("Sets are not lists")
    warning(msg_list)
  }
  if(!(length(pred_set)== length(true_set))){
    msg_length_list <- paste("Sets of different size")
    warning(msg_length_list)
  }
  n <- length(true_set)
  coverage <- matrix(0, nrow=n)
  for(i in 1:n){
    intersection <- length(intersect(true_set[[i]], pred_set[[i]]))
    coverage[i] <- ifelse(intersection>0,1,0)
  }
  mean(coverage)
}

#' Strict coverage for single observation
#'
#' Returns one if true set is strictly contained within
#' the predicted set zero otherwose.
#'
#' @param true_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param pred_set A `list` of numeric or character vectors representing the predicted sets.
#' 
#' @return A numeric value representing the mean coverage (proportion of
#'   intersected elements over true set size) across all observations.
#' @export
#'
#' @examples
#' true <- list(c(1, 2), c(3))
#' pred <- list(c(1), c(3, 4))
#' coverage_strict(true, pred)

coverage_strict_single <- function(true_set, pred_set) {
  if (all(true_set %in% pred_set)) return(1) else return(0)
}

#' Strict coverage
#'
#' Computes the average proportion of the true set that is strictly contained within
#' the predicted set. This is a normalized measure of recall across multiple
#' observations.
#'
#' @param true_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param pred_set A `list` of numeric or character vectors representing the predicted sets.
#'
#' @return A numeric value representing the mean coverage (proportion of
#'   intersected elements over true set size) across all observations.
#' @export
#'
#' @examples
#' true <- list(c(1, 2), c(3))
#' pred <- list(c(1), c(3, 4))
#' coverage_strict(true, pred)
coverage_strict <- function(true_set, pred_set){
  if(!(is.list(pred_set) & is.list(true_set))){
    msg_list <- paste("Sets are not lists")
    warning(msg_list)
  }
  if(!(length(pred_set)== length(true_set))){
    msg_length_list <- paste("Sets of different size")
    warning(msg_length_list)
  }
  n <- length(true_set)
  coverage <- matrix(0, nrow=n)
  for(i in 1:n){
    coverage[i] <- coverage_strict_single(true_set = true_set[[i]], 
                                          pred_set= pred_set[[i]])
  }
  mean(coverage)
}

#' Plug-in set-policy value
#'
#' Computes the plug-in set-policy value of a set-valued policy using 
#' estimated conditional means.
#'
#' @param test_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param test Data frame used for prediction (calibration/test set).
#' @param Q.all.actions A data frame containing the estimated conditional means.
#' @param gAX.pred  A data frame containing the estimated propensity scores
#' @param levels Vector of possible treatment/action levels. Defaults to `1:4`.
#' @param zero_indexed A logical indicator for binary treatments. 
#'
#' @return A list with two numeric values representing the uniform set-policy 
#' value and propensity set-policy value.
#' @export
set_policy_value_plug_in <- function(test_set, test, Q.all.actions,
                                    gAX.pred, levels= 1:4, zero_indexed = FALSE) {

  n <- nrow(test)
  m<- length(levels)
  row_idx <- seq_len(n)
  col_offset <- (0:(m - 1)) * n

  test_set <- if (is.list(test_set)) test_set else as.list(test_set)
  if (zero_indexed) test_set <- lapply(test_set, function(s) s + 1L)
  
  Ind <- matrix(0, n, m)
  for (i in row_idx) Ind[i, test_set[[i]]] <- 1
  Ind[rowSums(Ind)==0,] <- 1
  
  q_unif <- Ind /rowSums(Ind)
  results <- rowSums(Q.all.actions * q_unif) |> mean()

  g_c <- rowSums(Ind * gAX.pred)
  q_p <- Ind * gAX.pred / g_c
  
  results_non_random <- rowSums(Q.all.actions * q_p) |> mean()
  
  return(list(results,results_non_random))

}

#' TMLE set-policy value
#'
#' Computes the TMLE set-policy value of a set-valued policy using 
#' estimated conditional means and propensity score.
#'
#' @param test_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param Q.all.actions A data frame containing the estimated conditional means.
#' @param gAX.pred  A data frame containing the estimated propensity scores.
#' @param Y A vector containing the outcome variable. 
#' @param A A vector containing the treatment variable.
#' @param ab A vector containing the minimal and maximal value for Y.
#' @param zero_indexed A logical indicator for binary treatments. 
#'
#' @return A list with two numeric values representing the estimated uniform
#'  set-policy value and propensity set-policy value.
#' @export
set_policy_value_tmle <- function(test_set, Q.all.actions, gAX.pred, Y, A, ab, 
                                  levels= 1:4, zero_indexed = FALSE, z= qnorm(0.95)) {
  
  n <- length(Y) 
  m <- ncol(Q.all.actions)
  row_idx <- seq_len(n)
  A_col <- as.integer(setNames(seq_along(levels), levels)[as.character(A)])
  
  Y01 = (Y - min(ab))/diff(ab)
  Qd01 = ((Q.all.actions - min(ab))/diff(ab)) |> as.matrix()
  Qd01[Qd01 > 0.999] <- 0.999
  Qd01[Qd01 < 0.001] <- 0.001

  logit.Qd01 = stats::qlogis(Qd01)
  
  conf <- if (is.list(test_set)) test_set else as.list(test_set)
  if (zero_indexed) conf <- lapply(conf, function(s) s + 1L)   # -> 1-indexed
  
  g <- gAX.pred
  g <- pmax(g, 0.01)
  g <- g / rowSums(g)
  

  tmle <- function(q_mat, H_mat, propensity = FALSE) {
    H_obs <- H_mat[cbind(row_idx, A_col)]
    logit_Q_obs <- logit.Qd01[cbind(row_idx, A_col)]
    update <- stats::glm(
      Y01 ~ -1 + H_obs + offset(logit_Q_obs), family = "quasibinomial")
    e <- update$coefficients
    Qd01_star = stats::plogis(logit.Qd01 + e * H_mat)
    Qd_stars  = Qd01_star * diff(ab) + min(ab)
    m_star <- rowSums(Qd_stars * q_mat)

    phi <- if (propensity){
      m_star + H_obs * (Y - m_star)
    }else{
        m_star + H_obs * (Y - Qd_stars[cbind(row_idx, A_col)])
    }
    
    out <- mean(m_star) 
    se  <- sd(phi) / sqrt(n)
    attr(out, "se")   <- se
    attr(out, "low")  <- out - z * se
    attr(out, "high") <- out + z * se
    out
  }
  
  # uniform SPV
  Ind <- matrix(0, n, m)
  for (i in row_idx) Ind[i, conf[[i]]] <- 1
  Ind[rowSums(Ind)==0,] <- 1
  
  q_unif <- Ind /rowSums(Ind)
  res_unif <- tmle(q_unif, q_unif / g)
  
  # propensity SPV
  g <- gAX.pred
  g_c <- rowSums(Ind * g)
  g_c <- pmax(g_c, 0.01)
  q_p <- Ind * g / g_c
  res_prop <- tmle(q_p, Ind / g_c, propensity = TRUE)
  
  list(results = res_unif, results_non_random = res_prop)
}

#' AIPW set-policy value
#'
#' Computes the AIPW set-policy value of a set-valued policy using 
#' estimated conditional means and propensity score.
#'
#' @param test_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param Q.all.actions A data frame containing the estimated conditional means.
#' @param gAX.pred  A data frame containing the estimated propensity scores.
#' @param Y A vector containing the outcome variable. 
#' @param A A vector containing the treatment variable.
#' @param levels Vector of possible treatment/action levels. Defaults to `1:4`.
#' @param zero_indexed A logical indicator for binary treatments. 
#'
#' @return A list with two numeric values representing the estimated uniform 
#' set-policy value and propensity set-policy value.
#' @export
set_policy_value_aipw <- function(test_set, Q.all.actions, gAX.pred, Y, A, 
                                  levels= 1:4, zero_indexed = FALSE, z= qnorm(0.95)) {
  
  n <- nrow(Q.all.actions)
  m<- length(levels)
  row_idx <- seq_len(n)
  col_offset <- (0:(m - 1)) * n
  A_col <- as.integer(setNames(seq_along(levels), levels)[as.character(A)])
  
  test_set <- if (is.list(test_set)) test_set else as.list(test_set)
  if (zero_indexed){
    test_set <- lapply(test_set, function(s) s + 1L)
  }
  
  Ind <- matrix(0, n, m)
  for (i in row_idx) Ind[i, test_set[[i]]] <- 1
  Ind[rowSums(Ind)==0,] <- 1
  
  ci <- function(phi) {
    out <- mean(phi)
    se  <- sd(phi) / sqrt(n)
    attr(out, "se")   <- se
    attr(out, "low")  <- out - z * se
    attr(out, "high") <- out + z * se
    out}
  
  q_unif <- Ind /rowSums(Ind)
  g_obs <- gAX.pred[cbind(row_idx, A_col)]
  Q.obs <-  Q.all.actions[cbind(row_idx, A_col)]
  
  one_step_unif <- rowSums(Q.all.actions * q_unif) + 
    (q_unif[cbind(row_idx, A_col)]/g_obs)*(Y- Q.obs)
  
  g_c <- rowSums(Ind * gAX.pred)
  q_p    <- Ind * gAX.pred / g_c
  w <- Ind[cbind(row_idx, A_col)] / g_c
  
  m_est_prop <- rowSums(Q.all.actions * q_p)
  one_step_prop <- m_est_prop + w*(Y-m_est_prop)

  return(list(ci(one_step_unif), ci(one_step_prop)))
}

#' Set-policy values for IVF data example
#'
#' Estimates the uniform set-policy value for a primary outcome (Y) and an
#' adverse event (xi). Additionally computes the value for a "minimal treatment"
#' strategy—selecting the lowest available treatment level across all cases.
#' **Note:** This function is uses the package `SL.ODTR` from
#' Montoya, L. M., van der Laan, M. J., Luedtke, A. R., Skeem, J. L., Coyle, J. R.,
#' & Petersen, M. L. (2023). The optimal dynamic treatment rule superlearner:
#' considerations, performance, and application to criminal justice interventions.
#'
#' @param test_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param test Data frame used for prediction (calibration/test set).
#' @param covariates Character vector of covariate names. Defaults to `c("x1", "x2")`.
#' @param treatment_name String indicating the treatment variable. Defaults to "A".
#' @param outcome_name String indicating the outcome variable. Defaults to "Y".
#' @param second_outcome String indicating the second outcome variable. Defaults to "xi".
#' @param mod_y Model to predict the conditional mean outcome.
#' @param mod_xi Model to predict the conditional mean of the second outcome.
#' @param mod_ps Model to predict the propensity score.
#' @param ab Float indicating the largest difference in the outcome.
#' @param ab_xi Float indicating the largest difference in the second outcome.
#' @param n_test Integer indicating the number of policies to sample at random 
#' for the set-valued policy. Defaults to 1.
#' @param levels Vector of possible treatment/action levels. Defaults to `1:5`.
#'
#' @return A list with the estimated set-policy values (random and lowest
#' strategy for Y and xi) of the set-valued policy.
#' @export
ivf_set_policy_values <- function(test_set, test,
                                  covariates = c("x1","x2"),
                                  treatment_name = "A",
                                  outcome_name = "Y",
                                  second_outcome ="xi",
                                  mod_y, mod_xi, mod_ps,
                                  ab, ab_xi, n_test=1, levels) {

  if(!is.list(test_set)){
    test_set <- as.list(test_set)
  }
  n <- nrow(test)
  m <-length(levels)
  row_idx <- seq_len(n)
  col_offset <- (0:(m - 1)) * n
  random_policy <- matrix(NA_integer_, n, n_test)
  lowest_policy <- matrix(NA_integer_, n, n_test)
  for (i in seq_len(n)) {
    allowed <- test_set[[i]]
    if (length(allowed) > 0) {
      random_policy[i, ] <- allowed[sample.int(length(allowed), n_test, replace = TRUE)]
      lowest_policy[i, ] <- min(allowed)
    } else {
        random_policy[i, ] <- sample.int(m, n_test, replace = TRUE)
        lowest_policy[i, ] <- 1
    }}

  gAX.pred <- stats::predict(
    mod_ps,
    newdata = test[, covariates, drop = FALSE],
    type = "prob")
  gAW_bounded <- pmax(gAX.pred, 0.01)

  base_newdata <- test[, covariates, drop = FALSE]

  get_Q <- function(mod) {
    sapply(levels, function(a) {
      newdata_temp <- base_newdata
      newdata_temp[[treatment_name]] <- factor(a, levels = levels)
      stats::predict(mod, newdata = newdata_temp, type = "response")$pred
    })
  }

  Q_all_Y  <- get_Q(mod_y)
  Q_all_xi <- get_Q(mod_xi)

  # 2. Extract common variables
  test_A <- test[,treatment_name]
  Y_vec  <- test[,outcome_name]
  xi_vec <- test[,second_outcome]

  compute_psi <- function(policy_mat, outcome_vec, Q_mat, ab_vec) {
    # Overhead check: Don't fork processes if n_test is 1
    if (n_test <= 1) {
      d <- policy_mat[, 1]
      lin_idx <- row_idx + col_offset[d]
      return(SL.ODTR::tmle.d.fun(A = test_A, Y = outcome_vec, d = d,
                                 Qd = Q_mat[lin_idx], gAW = gAW_bounded[lin_idx],
                                 ab = ab_vec)$psi)
    }

    unlist(parallel::mclapply(seq_len(n_test), function(p) {
      d <- policy_mat[, p]
      lin_idx <- row_idx + col_offset[d]
      SL.ODTR::tmle.d.fun(A = test_A, Y = outcome_vec, d = d,
                          Qd = Q_mat[lin_idx], gAW = gAW_bounded[lin_idx],
                          ab = ab_vec)$psi
    }, mc.cores = parallel::detectCores()))
  }

  list(
    results_random_Y  = compute_psi(random_policy, Y_vec,  Q_all_Y, ab),
    results_random_xi = compute_psi(random_policy, xi_vec, Q_all_xi, ab_xi),
    results_min_Y     = compute_psi(lowest_policy, Y_vec,  Q_all_Y, ab),
    results_min_xi    = compute_psi(lowest_policy, xi_vec, Q_all_xi, ab_xi)
  )
}

#' Margin Nonconformity Score
#'
#' Generates the margin nonconformity scores for a matrix of potential outcomes.
#' The score is calculated as the difference between the maximum potential
#' outcome and all other outcomes, with a specific "margin" calculation for
#' the winning class.
#'
#' @param potential_outcomes Matrix of potential outcomes (observations by treatments).
#'
#' @return A matrix of the same dimensions as `potential_outcomes` containing
#'   the margin nonconformity scores.
#' @export
#'
#' @examples
#' margin_score(matrix(runif(10 * 5), 10, 5))
margin_score <- function(potential_outcomes) {
  which_max <- max.col(potential_outcomes, ties.method = "first")
  row_indices <- seq_len(nrow(potential_outcomes))
  max_vals <- potential_outcomes[cbind(row_indices, which_max)]
  temp_outcomes <- potential_outcomes
  temp_outcomes[cbind(row_indices, which_max)] <- -Inf
  second_max_vals <- apply(temp_outcomes, 1, max)
  score_matrix <- max_vals - potential_outcomes
  score_matrix[cbind(row_indices, which_max)] <- second_max_vals - max_vals
  return(score_matrix)
}

#' Complete evaluation of a set-valued policy in the synthetic settings.
#'
#'
#'
#' @param test_set A `list` of numeric or character vectors representing the 
#' predicted sets.
#' @param optimal_policy_new A `list` of numeric or character vectors representing 
#' the true optimal sets.
#' @param prop_score_new A data frame containing the true propensity scores.
#' @param potential_outcomes A data frame containing the true potential outcomes.
#' @param df_new_sample Data frame used for prediction.
#' @param levels_A Vector of possible treatment/action levels. Defaults to `1:5`.
#' @param treatment_name String indicating the treatment variable. Defaults to "A".
#' @param outcome_name String indicating the outcome variable. Defaults to "Y".
#'
#' @return A list containing different evaluation metrics: exact match, 
#' strict and relaxed coverage, and set-policy values (uniform and propensity). 
#' @export
#'
#' @examples
#' margin_score(matrix(runif(10 * 5), 10, 5))
table.evaluation <- function(test_set, optimal_policy_new,
                             prop_score_new, potential_outcomes, 
                             df_new_sample, levels_A,
                             treatment_name = "A", outcome_name = "Y"){
  exact.matches <- sapply(1:length(test_set), function(i) {
    setequal(test_set[[i]], optimal_policy_new[[i]])|> 
      as.numeric()
  })
  
  cardinality.mean <- sapply(1:length(test_set), function(i) {
    length(test_set[[i]])|> 
      as.numeric()
  }) |> mean()
  
  cov <- sapply(1:length(test_set), function(i) {
    coverage_strict_single(pred_set = test_set[[i]], true_set = optimal_policy_new[[i]])})
  
  cov.relaxed <- coverage_relaxed(true_set = optimal_policy_new, 
                                  pred_set = test_set)
  
  spv.mean <- set_policy_value_plug_in(test_set = test_set, 
                                       test = df_new_sample, 
                                       Q.all.actions = potential_outcomes, 
                                       gAX.pred = prop_score_new, 
                                       levels = levels_A)
  spv.unif <- spv.mean[[1]]
  spv.propensity <- spv.mean[[2]]
  return(list(exact.matches, cardinality.mean, cov, cov.relaxed, spv.unif, spv.propensity))
}

#' Complete evaluation of a set-valued policy using real data.
#'
#'
#'
#' @param test_set A `list` of numeric or character vectors representing the 
#' predicted sets.
#' @param prop_score_new A data frame containing the true propensity scores.
#' @param potential_outcomes A data frame containing the true potential outcomes.
#' @param df_new_sample Data frame used for prediction.
#' @param levels_A Vector of possible treatment/action levels. Defaults to `1:5`.
#' @param ab A vector containing the minimal and maximal value for Y.
#' @param zero_indexed A logical indicator for binary treatments. 
#' @param treatment_name String indicating the treatment variable. Defaults to "A".
#' @param outcome_name String indicating the outcome variable. Defaults to "Y".
#'
#' @return A list containing different evaluation metrics: exact match, 
#' strict and relaxed coverage, and set-policy values (uniform and propensity). 
#' @export
table.evaluation.real <- function(test_set,
                             prop_score_new, potential_outcomes, 
                             df_new_sample, levels_A, ab, zero_indexed=FALSE,
                             treatment_name = "A", outcome_name = "Y"){
  
  cardinality.mean <- sapply(1:length(test_set), function(i) {
    length(test_set[[i]])|> 
      as.numeric()}) |> mean()
  
  spv.aipw <- set_policy_value_aipw(test_set = test_set, 
                                    Y = df_new_sample[,outcome_name], 
                                    A = df_new_sample[,treatment_name],
                                    Q.all.actions = potential_outcomes, 
                                    gAX.pred = prop_score_new,
                                    levels = levels_A, zero_indexed = zero_indexed)
  spv.unif.aipw <- spv.aipw[[1]]
  spv.propensity.aipw <- spv.aipw[[2]]
  
  spv.tmle <- set_policy_value_tmle(test_set = test_set, 
                                    Y = df_new_sample[,outcome_name], 
                                    A = df_new_sample[,treatment_name],
                                    Q.all.actions = potential_outcomes, 
                                    gAX.pred = prop_score_new, ab = ab,
                                    levels = levels_A, zero_indexed = zero_indexed)
  spv.unif.tmle <- spv.tmle[[1]]
  spv.propensity.tmle <- spv.tmle[[2]]
  
  results_df <- rbind(
    data.frame(
      estimator = "AIPW",
      policy = "Uniform",
      value = spv.unif.aipw,       # adjust $value to match your object structure
      lower = attr(spv.unif.aipw, "low"),
      upper = attr(spv.unif.aipw, "high"),
      stringsAsFactors = FALSE),
    data.frame(
      estimator = "AIPW",
      policy = "Propensity",
      value = spv.propensity.aipw,
      lower = attr(spv.propensity.aipw,"low"),
      upper = attr(spv.propensity.aipw,"high"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      estimator = "TMLE",
      policy = "Uniform",
      value = spv.unif.tmle,
      lower = attr(spv.unif.tmle, "low"),
      upper = attr(spv.unif.tmle, "high"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      estimator = "TMLE",
      policy = "Propensity",
      value = spv.propensity.tmle,
      lower = attr(spv.propensity.tmle, "low"),
      upper = attr(spv.propensity.tmle, "high"),
      stringsAsFactors = FALSE
    )
  )
  
  
  
  return(list(cardinality.mean, results_df))
}
