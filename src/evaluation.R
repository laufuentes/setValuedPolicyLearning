#' Average Jaccard distance for prediction sets
#'
#' Computes the average Jaccard Distance between a true set of labels and 
#' a predicted set of labels. The distance is defined as `1 - (intersection / union)`.
#'
#' @param true_set A `list` of numeric or character vectors representing the ground truth sets.
#' @param pred_set A `list` of numeric or character vectors representing the predicted sets.
#' @param levels Vector of possible treatment levels. Defaults to `1:5`.
#'
#' @return A numeric value representing the mean Jaccard Distance across all observations.
#' @export
#'
#' @examples
#' true <- list(c(1, 2), c(3))
#' pred <- list(c(1, 3), c(3))
#' jaccard_distance(true, pred)
jaccard_distance <- function(true_set, pred_set, levels=1:5){
  if(!(is.list(pred_set) & is.list(true_set))){
    msg_list <- paste("Sets are not lists")
    warning(msg_list)
  }
  if(!(length(pred_set)== length(true_set))){
    msg_length_list <- paste("Sets of different size")
    warning(msg_length_list)
  }
  n <- length(true_set)
  jac_index <- matrix(0, nrow=n)
  for(i in 1:n){
    intersection <- length(intersect(true_set[[i]], pred_set[[i]]))
    union <- length(union(true_set[[i]], pred_set[[i]]))
    jac_index[i] <- 1 - (intersection / union)
  }
  mean(jac_index)
}

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

#' Coverage
#'
#' Computes the average proportion of the true set that is contained within 
#' the predicted set. This represents the average recall across all observations.
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
#' coverage(true, pred)
coverage <- function(true_set, pred_set){
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
    size_true <- length(true_set[[i]])
    coverage[i] <- intersection / size_true
  }
  mean(coverage)
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
#' coverage_unif(true, pred)
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

set_policy_value<- function(test_set, test, 
                            QAW.reg.train, g.reg.train, 
                            ab, covariates=c("x1","x2"), 
                            family="gaussian", levels){
  random_policy <- apply(data.frame(1:nrow(test)), 1, 
                         function(i)if(length(test_set[[i]])>0){
                           sample(test_set[[i]],1)
                           }else{sample(levels,1)})
  gAW.pred = stats::predict(g.reg.train,
                            newdata = data.frame(test[,covariates]), 
                            type = "prob")
  gAW_bounded <- pmax(gAW.pred, 0.01)
  
  newdata_temp <- data.frame(test[,covariates[1:2]],
                             A = factor(random_policy, 
                                        levels=levels))
  Qd.test = stats::predict(QAW.reg.train, 
                           newdata = newdata_temp, 
                           type = "response")$pred
  
  tmle.obj.test = SL.ODTR::tmle.d.fun(A = test$A, 
                                      Y = matrix(test$Y),
                                      d = random_policy, 
                                      Qd = Qd.test, 
                                      gAW = gAW_bounded[cbind(1:nrow(test),random_policy)], 
                                      ab = ab)
  return(tmle.obj.test$psi)
}

bounds_set_policy_value <- function(test_set, test, QAW.reg.train, 
                                         g.reg.train, ab, n_test = 100,
                                         covariates = c("x1","x2"),
                                         levels, m=5) {
  
  n <- nrow(test)
  row_idx <- seq_len(n)
  col_offset <- (0:(m - 1)) * n
  
  ## ---- 1. FAST policy sampling ----
  random_policy <- matrix(NA_integer_, n, n_test)
  
  for (i in seq_len(n)) {
    allowed <- test_set[[i]]
    
    if (length(allowed) > 0) {
      random_policy[i, ] <- sample(allowed, n_test, replace = TRUE)
    } else {
      random_policy[i, ] <- sample.int(m, n_test, replace = TRUE)
    }
  }
  
  ## ---- 2. Predict g once ----
  gAW.pred <- stats::predict(
    g.reg.train,
    newdata = test[, covariates, drop = FALSE],
    type = "prob"
  )
  gAW_bounded <- pmax(gAW.pred, 0.01)
  
  ## ---- 3. Predict Q once for all actions ----
  base_newdata <- test[, covariates, drop = FALSE]
  
  Q_all_actions <- sapply(seq_len(m), function(a) {
    newdata_temp <- base_newdata
    newdata_temp$A <- factor(a, levels = seq_len(m))
    stats::predict(QAW.reg.train, 
                   newdata = newdata_temp, 
                   type = "response")$pred
  })
  
  ## ---- 4. Compute result ----
  Y_mat <- matrix(test$Y)
  
  results <- unlist(
    parallel::mclapply(seq_len(n_test), function(p) {
      d <- random_policy[, p]
      lin_idx <- row_idx + col_offset[d]
      
      SL.ODTR::tmle.d.fun(
        A   = test$A,
        Y   = Y_mat,
        d   = d,
        Qd  = Q_all_actions[lin_idx],
        gAW = gAW_bounded[lin_idx],
        ab  = ab)$psi
    }, mc.cores = parallel::detectCores())
  )
  
  results
}