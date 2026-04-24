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
coverage_relaxed <- function(true_set, pred_set, levels=1:5){
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
coverage_strict <- function(true_set, pred_set, levels=1:5){
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
    coverage[i] <- coverage_strict_single(true_set = true_set[[i]], pred_set= pred_set[[i]])
  }
  mean(coverage)
}

#' Average width of prediction set 
#'
#' Computes the mean width the predicted set. 
#'
#' @param true_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param levels Vector of possible treatment levels. Defaults to `1:5`.
#'
#' @return A numeric value representing the mean of the prediction set across all observations.
#' @export
#'
#' @examples
#' pred <- list(c(1), c(3, 4))
#' width(pred)
width <- function(pred_set, levels=1:5){
  if(!(is.list(pred_set))){
    msg_list <- paste("Sets are not lists")
    warning(msg_list)
  }
  n <- length(pred_set)
  m_width <- matrix(0, nrow=n)
  for(i in 1:n){
    m_width[i] <- length(pred_set[[i]])
  }
  mean(m_width)
}

policy_value_tmle<- function(random_policy, test, 
                            mod_y, mod_ps, 
                            ab, covariates=c("x1","x2"), 
                            treatment_name = "A", outcome_name="Y",
                            family="gaussian", levels){
  ## ---- 2. Predict g once ----
  gAW.pred <- stats::predict(
    mod_ps,
    newdata = test[, covariates, drop = FALSE],
    type = "prob"
  )
  gAW_bounded <- pmax(gAW.pred, 0.01)
  
  ## ---- 3. Predict Q once for all actions ----
  base_newdata <- test[, covariates, drop = FALSE]
  Q_all_actions <- sapply(levels, function(a) {
    newdata_temp <- base_newdata
    newdata_temp[, treatment_name] <- factor(a, levels = levels)
    stats::predict(mod_y, 
                   newdata = newdata_temp, 
                   type = "response")$pred
  })
  
  ## ---- 4. Compute result ----
  Y_mat <- matrix(test[, outcome_name])
  tmle.obj.test = SL.ODTR::tmle.d.fun(A = test[,treatment_name], 
                                      Y = Y_mat,
                                      d = random_policy, 
                                      Qd =  Q_all_actions[cbind(1:nrow(test),random_policy)], 
                                      gAW = gAW_bounded[cbind(1:nrow(test),random_policy)], 
                                      ab = ab)
  return(tmle.obj.test$psi)
}

bounds_set_policy_value <- function(test_set, test, 
                                    mod_y,
                                    treatment_name = "A",
                                    outcome_name = "Y",
                                    mod_ps, ab, n_test = 1,
                                    covariates = c("x1","x2"),
                                    levels) {
  
  n <- nrow(test)
  m <-length(levels)
  row_idx <- seq_len(n)
  col_offset <- (0:(m - 1)) * n
  
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
  
  ## ---- 2. Predict g once ----
  gAW.pred <- stats::predict(
    mod_ps,
    newdata = test[, covariates, drop = FALSE],
    type = "prob"
  )
  gAW_bounded <- pmax(gAW.pred, 0.01)
  
  ## ---- 3. Predict Q once for all actions ----
  base_newdata <- test[, covariates, drop = FALSE]
  
  Q_all_actions <- sapply(levels, function(a) {
    newdata_temp <- base_newdata
    newdata_temp[, treatment_name] <- factor(a, levels = levels)
    stats::predict(mod_y, 
                   newdata = newdata_temp, 
                   type = "response")$pred
  })
  
  ## ---- 4. Compute result ----
  Y_mat <- matrix(test[, outcome_name])
  
  results <- unlist(
    parallel::mclapply(seq_len(n_test), function(p) {
      d <- random_policy[, p]
      lin_idx <- row_idx + col_offset[d]
      
      SL.ODTR::tmle.d.fun(
        A   = test[, treatment_name],
        Y   = Y_mat,
        d   = d,
        Qd  = Q_all_actions[lin_idx],
        gAW = gAW_bounded[lin_idx],
        ab  = ab)$psi
    }, mc.cores = parallel::detectCores())
  )
  
  results
}

oracular_set_policy_value <- function(test_set, test, test_potential_outcome,
                                    treatment_name = "A",
                                    outcome_name = "Y",
                                     n_test = 1,
                                    covariates = c("x1","x2"),
                                    levels) {
  
  n <- nrow(test)
  m<- length(levels)
  row_idx <- seq_len(n)
  col_offset <- (0:(m - 1)) * n
  
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
  
  results
}


ivf_set_policy_values <- function(test_set, test, 
                                    mod_y, mod_xi, mod_ps, 
                                    treatment_name = "A",
                                    outcome_name = "Y", second_outcome ="xi",
                                    ab, levels, n_test=1,
                                    covariates = c("x1","x2")) {
  
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
  
  compute_psi <- function(policy_mat, outcome_vec, Q_mat) {
    # Overhead check: Don't fork processes if n_test is 1
    if (n_test <= 1) {
      d <- policy_mat[, 1]
      lin_idx <- row_idx + col_offset[d]
      return(SL.ODTR::tmle.d.fun(A = test_A, Y = outcome_vec, d = d, 
                                 Qd = Q_mat[lin_idx], gAW = gAW_bounded[lin_idx], 
                                 ab = ab)$psi)
    }
  
    unlist(parallel::mclapply(seq_len(n_test), function(p) {
      d <- policy_mat[, p]
      lin_idx <- row_idx + col_offset[d]
      SL.ODTR::tmle.d.fun(A = test_A, Y = outcome_vec, d = d, 
                          Qd = Q_mat[lin_idx], gAW = gAW_bounded[lin_idx], 
                          ab = ab)$psi
    }, mc.cores = parallel::detectCores()))
  } 
    
  list(
    results_random_Y  = compute_psi(random_policy, Y_vec,  Q_all_Y),
    results_random_xi = compute_psi(random_policy, xi_vec, Q_all_xi),
    results_min_Y     = compute_psi(lowest_policy, Y_vec,  Q_all_Y),
    results_min_xi    = compute_psi(lowest_policy, xi_vec, Q_all_xi)
  )
}
