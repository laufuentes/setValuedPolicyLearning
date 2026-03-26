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

#' Adaptive Prediction Sets (APS) Nonconformity Score 
#' From: Romano, Y., Sesia, M., & Candes, E. (2020). 
#' Classification with valid and adaptive coverage. 
#'
#' Generates the APS nonconformity scores for a matrix of probabilities. This 
#' score uses the cumulative sum of sorted probabilities and a uniform 
#' random variable for tie-breaking.
#'
#' @param prob_matrix Matrix of probabilities (observations by treatments/classes).
#'
#' @return A matrix of the same dimensions as `prob_matrix` containing 
#'   the APS nonconformity scores.
#' @export
#'
#' @examples
#' APS_score(matrix(runif(10 * 5), 10, 5))
APS_score <- function(prob_matrix) {
  n <- nrow(prob_matrix)
  K <- ncol(prob_matrix)
  
  scores <- matrix(0, n, K)
  
  for (i in 1:n) {
    p <- prob_matrix[i, ]
    U <- runif(K)   # independent U for each class
    
    for (a in 1:K) {
      p_a <- p[a]
      scores[i, a] <- sum(p[p > p_a]) + p_a * U[a]
    }
  }
  
  return(scores)
}

#' Weighted Score for a Specific Action
#'
#' Computes the negative weighted probability (loss score) for a specific 
#' treatment action $a$ across a set of experts.
#'
#' @param df_pred Data frame containing the covariates for prediction.
#' @param fitted_experts Matrix where rows are observations and columns are 
#'   individual expert policy predictions.
#' @param weights Numeric vector of weights assigned to each expert. 
#'   Length must match the number of columns in `fitted_experts`.
#' @param a The specific treatment action level to evaluate.
#'
#' @return A numeric vector of length `nrow(df_pred)` representing the 
#'   negative weighted score for action `a`.
#' @export
weighted_score_experts_a <- function(df_pred, fitted_experts, weights, a){
  n <- nrow(df_pred)
  pi_ax <- ifelse(fitted_experts==a,1,0)
  scores <- -as.vector(weights%*%t(pi_ax))
  return(scores)
}
