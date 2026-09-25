#' Average uniform coverage
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
#' coverage_unif(true, pred)
coverage_unif <- function(true_set, pred_set) {
  n <- length(true_set)
  total <- 0

  for (i in seq_len(n)) {
    ti <- true_set[[i]]
    if (length(ti) == 0) next

    total <- total +
      length(intersect(ti, pred_set[[i]])) / length(ti)
  }

  total / n
}

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

#' Average width of prediction set
#'
#' Computes the mean width the predicted set.
#'
#' @param pred_set A `list` of numeric or character vectors representing the predicted sets.
#' @param levels Vector of possible treatment levels. Defaults to `1:5`.
#'
#' @return A numeric value representing the mean of the prediction set across all observations.
#' @export
#'
#' @examples
#' pred <- list(c(1), c(3, 4))
#' width(pred)
width <- function(pred_set, levels=1:5){
  if(!is.list(pred_set)){
    coords <- which(pred_set == 1, arr.ind = TRUE)
    pred_set <- split(coords[, "col"], coords[, "row"])
  }
  n <- length(pred_set)
  m_width <- matrix(0, nrow=n)
  for(i in 1:n){
    m_width[i] <- length(pred_set[[i]])
  }
  mean(m_width)
}

#' Oracular set-policy value
#'
#' Computes the oracular set-policy value of a set-valued policy using potential outcomes.
#' **Note:** This function is uses the package `SL.ODTR` from
#' Montoya, L. M., van der Laan, M. J., Luedtke, A. R., Skeem, J. L., Coyle, J. R.,
#' & Petersen, M. L. (2023). The optimal dynamic treatment rule superlearner:
#' considerations, performance, and application to criminal justice interventions.
#'
#' @param test_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param test Data frame used for prediction (calibration/test set).
#' @param test_potential_outcome A data frame containing the potential outcomes.
#' @param covariates Character vector of covariate names. Defaults to `c("x1", "x2")`.
#' @param treatment_name String indicating the treatment variable. Defaults to "A".
#' @param outcome_name String indicating the outcome variable. Defaults to "Y".
#' @param mod_ps Model to predict the propensity score trained on train1, used for modeling doctor's predictions.
#' @param n_test Integer indicating the number of policies to sample at random for the set-valued policy.
#' @param levels Vector of possible treatment/action levels. Defaults to `1:5`.
#'
#' @return A numeric value representing the oracular set-policy value of the set-valued policy.
#' @export
oracular_set_policy_value <- function(test_set, test, test_potential_outcome,
                                    covariates = c("x1","x2"),
                                    treatment_name = "A",
                                    outcome_name = "Y", 
                                    gAW.pred=NULL,
                                    n_test = 1, levels= 1:5) {

  n <- nrow(test)
  m<- length(levels)
  row_idx <- seq_len(n)
  col_offset <- (0:(m - 1)) * n

  if(!is.list(test_set)){
    test_set <- as.list(test_set)
    # coords <- which(test_set == 1, arr.ind = TRUE)
    # test_set <- split(coords[, "col"], coords[, "row"])
  }
  
  ## ---- 1. FAST policy sampling ----
    random_policy <- matrix(NA_integer_, n, n_test)
    
    for (i in seq_len(n)) {
      allowed <- test_set[[i]]
      
      if (length(allowed) > 0) {
        random_policy[i, ] <- allowed[sample.int(length(allowed), n_test, replace = TRUE)]
      } else {
        random_policy[i, ] <- sample.int(m, n_test, replace = TRUE)
      }
    } 
  
  results <- unlist(
    parallel::mclapply(seq_len(n_test), function(p) {
      d <- random_policy[, p]
      mean(test_potential_outcome[cbind(1:n,d)])
    }, mc.cores = parallel::detectCores())
  )
  
  if(!is.null(gAW.pred)){
    #gAW.pred <- stats::predict(
    #  mod_ps,
    #  newdata = test[, covariates, drop = FALSE],
    #  type = "prob")
    gAW_bounded <- pmax(gAW.pred, 0.01)
    
    non_random_policy <- matrix(NA_integer_, n, n_test)
    for (i in seq_len(n)) {
      allowed <- test_set[[i]]
      
      if (length(allowed) > 0) {
        probs_allowed <- gAW_bounded[i, allowed]/sum(gAW_bounded[i, allowed])
        non_random_policy[i, ] <-  allowed[sample.int(length(allowed), n_test, replace = TRUE, prob = probs_allowed)]
      } else {
        non_random_policy[i, ] <- sample.int(m, n_test, replace = TRUE, prob = gAW_bounded[i,])
      }
    } 
    results_non_random <- unlist(
      parallel::mclapply(seq_len(n_test), function(p) {
        d <- non_random_policy[, p]
        mean(test_potential_outcome[cbind(1:n,d)])
      }, mc.cores = parallel::detectCores())
    )
    return(list(results,results_non_random))
  }
  return(results)
}

#' Set-policy value
#'
#' Estimates the uniform set-policy value of a set-valued policy.
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
#' @param mod_ps_spv Model to predict the propensity score trained on train1, used for modeling doctor's predictions.
#' @param mod_y Model to predict the conditional mean outcome (estimate potential outcomes).
#' @param mod_ps Model to predict the propensity score.
#' @param ab Float indicating the largest difference in the outcome.
#' @param n_test Integer indicating the number of policies to sample at random for the set-valued policy.
#' @param levels Vector of possible treatment/action levels. Defaults to `1:5`.
#'
#' @return A numeric value representing the estimated set-policy value of the set-valued policy.
#' @export
set_policy_value <- function(test_set, test,
                             covariates = c("x1", "x2"),
                             treatment_name = "A",
                             outcome_name = "Y",
                             Q_all_actions, gAW.pred, 
                             gAW.pred.spv = NULL, 
                             ab, n_test = 1, levels) {
  n <- nrow(test)
  m <- length(levels)
  row_idx <- seq_len(n)
  
  # Guarantee levels match 1..m indexing
  level_map <- setNames(seq_along(levels), levels)
  
  if (!is.list(test_set)) {
    test_set <- as.list(test_set)
  }
  
  A_vec <- test[[treatment_name]]
  Y_mat <- matrix(test[[outcome_name]])
  gAW_bounded <- as.matrix(pmax(gAW.pred, 0.01))
  Q_all_actions <- as.matrix(pmin(pmax(Q_all_actions, 0.001), 0.999))
  
  # 1. Uniform Policy Sampling
  random_policy <- matrix(NA_integer_, n, n_test)
  for (i in seq_len(n)) {
    allowed <- test_set[[i]]
    k <- length(allowed)
    if (k > 0) {
      random_policy[i, ] <- allowed[sample.int(k, n_test, replace = TRUE)]
    } else {
      random_policy[i, ] <- levels[sample.int(m, n_test, replace = TRUE)] %>% 
        as.numeric()
    }
  }
  
  # 2. TMLE Evaluator with Explicit Matrix Indexing [row, col]
  eval_tmle <- function(policy_mat, gAW_mat) {
    sapply(seq_len(n_test), function(p) {
      # Map chosen treatment actions to 1..m column indices
      d_raw <- policy_mat[, p]
      d_col <- level_map[as.character(d_raw)] %>% as.numeric()
      
      # Explicit 2D Matrix Indexing
      Q_d <- Q_all_actions[cbind(row_idx, d_col)]
      gAW_d <- gAW_mat[cbind(row_idx, d_col)]
      
      SL.ODTR::tmle.d.fun(
        A   = A_vec,
        Y   = Y_mat,
        d   = d_col,
        Qd  = Q_d,
        gAW = gAW_d,
        ab  = ab
      )$psi
    })
  }
  
  results <- eval_tmle(random_policy, gAW_bounded)
  
  # 3. Non-Random Policy
  if (!is.null(gAW.pred.spv)) {
    gAW_spv_bounded <- as.matrix(pmax(gAW.pred.spv, 0.01))
    
    non_random_policy <- matrix(NA_integer_, n, n_test)
    for (i in seq_len(n)) {
      allowed <- test_set[[i]]
      k <- length(allowed)
      if (k > 0) {
        col_idx <- level_map[as.character(allowed)] %>% as.numeric()
        probs_allowed <- gAW_spv_bounded[i, col_idx] / sum(gAW_spv_bounded[i, col_idx])
        non_random_policy[i, ] <- allowed[sample.int(k, n_test, replace = TRUE, prob = probs_allowed)]
      } else {
        probs_all <- gAW_spv_bounded[i, ] / sum(gAW_spv_bounded[i, ])
        non_random_policy[i, ] <- levels[sample.int(m, n_test, replace = TRUE, prob = probs_all)] %>% 
          as.numeric()
      }
    }
    
    results_non_random <- eval_tmle(non_random_policy, gAW_spv_bounded) 
    return(list(results, results_non_random))
  }
  
  return(results)
}

set_policy_value_plug_in <- function(test_set, test,
                             Q_all_actions,
                             gAW.pred.spv = NULL, 
                             ab, n_test = 1, levels) {
  n <- nrow(test)
  m <- length(levels)
  row_idx <- seq_len(n)
  
  # Guarantee levels match 1..m indexing
  level_map <- setNames(seq_along(levels), levels)
  
  if (!is.list(test_set)) {
    test_set <- as.list(test_set)
  }
  
  Q_all_actions <- as.matrix(pmin(pmax(Q_all_actions, 0.001), 0.999))
  
  # 1. Uniform Policy Sampling
  random_policy <- matrix(NA_integer_, n, n_test)
  for (i in seq_len(n)) {
    allowed <- test_set[[i]]
    k <- length(allowed)
    if (k > 0) {
      random_policy[i, ] <- allowed[sample.int(k, n_test, replace = TRUE)]
    } else {
      random_policy[i, ] <- levels[sample.int(m, n_test, replace = TRUE)] %>% 
        as.numeric()
    }
  }
  
  eval_plug.in <- function(policy_mat) {
    sapply(seq_len(n_test), function(p) {
      # Map chosen treatment actions to 1..m column indices
      d_raw <- policy_mat[, p]
      d_col <- level_map[as.character(d_raw)] %>% as.numeric()
      
      # Explicit 2D Matrix Indexing
      Q_all_actions[cbind(row_idx, d_col)] %>% mean()
      })
  }
  
  results <- eval_plug.in(random_policy)
  
  # 3. Non-Random Policy
  if (!is.null(gAW.pred.spv)) {
    gAW_spv_bounded <- as.matrix(pmax(gAW.pred.spv, 0.01))
    
    non_random_policy <- matrix(NA_integer_, n, n_test)
    for (i in seq_len(n)) {
      allowed <- test_set[[i]]
      k <- length(allowed)
      if (k > 0) {
        col_idx <- level_map[as.character(allowed)] %>% as.numeric()
        probs_allowed <- gAW_spv_bounded[i, col_idx] / sum(gAW_spv_bounded[i, col_idx])
        non_random_policy[i, ] <- allowed[sample.int(k, n_test, replace = TRUE, prob = probs_allowed)]
      } else {
        probs_all <- gAW_spv_bounded[i, ] / sum(gAW_spv_bounded[i, ])
        non_random_policy[i, ] <- levels[sample.int(m, n_test, replace = TRUE, prob = probs_all)] %>% 
          as.numeric()
      }
    }
    
    results_non_random <- eval_plug.in(non_random_policy) 
    return(list(results, results_non_random))
  }
  
  return(results)
}

set_policy_value_plug_in.bootstrap <- function(test_set, test,
                                     Q_all_actions,
                                     gAW.pred.spv = NULL, 
                                     ab, n_test = 1, levels) {
  n <- nrow(test)
  m <- length(levels)
  
  # Guarantee levels match 1..m indexing
  level_map <- setNames(seq_along(levels), levels)
  
  if (!is.list(test_set)) {
    test_set <- as.list(test_set)
  }
  
  Q_all_actions <- as.matrix(pmin(pmax(Q_all_actions, 0.001), 0.999))
  row_idx <- matrix(sample.int(n, size = n_test * (n %/% 2), replace = TRUE), 
                    ncol = n_test, nrow = n %/% 2)
  
  # 1. Uniform Policy Sampling
  random_policy <- matrix(NA_integer_, n)
  for (i in seq_len(n)) {
    allowed <- test_set[[i]]
    k <- length(allowed)
    if (k > 0) {
      random_policy[i] <- allowed[sample.int(k, 1, replace = TRUE)]
    } else {
      random_policy[i] <- levels[sample.int(m, 1, replace = TRUE)] %>% 
        as.numeric()
    }
  }
  
  eval_plug.in <- function(policy_mat) {
    sapply(seq_len(n_test), function(p) {
      # Map chosen treatment actions to 1..m column indices
      d_raw <- policy_mat[row_idx[,p],]
      d_col <- level_map[as.character(d_raw)] %>% as.numeric()
      
      # Explicit 2D Matrix Indexing
      Q_all_actions[cbind(row_idx[,p], d_col)] %>% mean()
    })
  }
  
  results <- eval_plug.in(random_policy)
  
  # 3. Non-Random Policy
  if (!is.null(gAW.pred.spv)) {
    gAW_spv_bounded <- as.matrix(pmax(gAW.pred.spv, 0.01))
    
    non_random_policy <- matrix(NA_integer_, n, n_test)
    for (i in seq_len(n)) {
      allowed <- test_set[[i]]
      k <- length(allowed)
      if (k > 0) {
        col_idx <- level_map[as.character(allowed)] %>% as.numeric()
        probs_allowed <- gAW_spv_bounded[i, col_idx] / sum(gAW_spv_bounded[i, col_idx])
        non_random_policy[i, ] <- allowed[sample.int(k, n_test, replace = TRUE, prob = probs_allowed)]
      } else {
        probs_all <- gAW_spv_bounded[i, ] / sum(gAW_spv_bounded[i, ])
        non_random_policy[i, ] <- levels[sample.int(m, n_test, replace = TRUE, prob = probs_all)] %>% 
          as.numeric()
      }
    }
    
    results_non_random <- eval_plug.in(non_random_policy) 
    return(list(results, results_non_random))
  }
  
  return(results)
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
#' @param n_test Integer indicating the number of policies to sample at random for the set-valued policy.
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
                                  ab, ab_xi, n_test=1,levels) {

  if(!is.list(test_set)){
    test_set <- as.list(test_set)
    # coords <- which(test_set == 1, arr.ind = TRUE)
    # test_set <- split(coords[, "col"], coords[, "row"])
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

  gAW.pred <- stats::predict(
    mod_ps,
    newdata = test[, covariates, drop = FALSE],
    type = "prob")
  gAW_bounded <- pmax(gAW.pred, 0.01)

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

table.evaluation <- function(test_set, optimal_policy_new,
                             prop_score_new, potential_outcomes, 
                             df_new_sample, levels_A, covariates_name, 
                             treatment_name = "A", outcome_name = "Y"){
  exact.matches <- sapply(1:length(test_set), function(i) {
    setequal(test_set[[i]], optimal_policy_new[[i]])%>% 
      as.numeric()
  })
  
  cardinality.mean <- sapply(1:length(test_set), function(i) {
    length(test_set[[i]])%>% 
      as.numeric()
  }) %>% mean()
  
  cov <- sapply(1:length(test_set), function(i) {
    coverage_strict_single(pred_set = test_set[[i]], true_set = optimal_policy_new[[i]])})
  
  spv.mean <- oracular_set_policy_value(test_set = test_set, 
                                        test = df_new_sample, 
                                        test_potential_outcome = potential_outcomes, 
                                        covariates = covariates_name, 
                                        treatment_name = treatment_name, 
                                        outcome_name = outcome_name, 
                                        gAW.pred = prop_score_new, 
                                        n_test = 1, levels = levels_A)
  spv.unif <- spv.mean[[1]]
  spv.propensity <- spv.mean[[2]]
  return(list(exact.matches, cardinality.mean, cov, spv.unif, spv.propensity))
}

