
#' Weighted Probability Distribution from Experts
#'
#' Computes the probability scores associated with covariates across a set of 
#' experts given a vector of expert weights.
#'
#'
#' @param df_pred Data frame containing the covariates for prediction.
#' @param fitted_experts Matrix where rows are observations and columns are 
#'   individual expert policy predictions.
#' @param weights Numeric vector of weights assigned to each expert. 
#'   Length must match the number of columns in `fitted_experts`.
#' @param levels Vector of possible treatment/action levels. Defaults to `1:5`.
#'
#' @return A numeric `matrix` of dimensions `nrow(df_pred)` by `length(levels)` 
#'   where each cell `(i, a)` represents the weighted probability of choosing 
#'   action `a` for observation `i`.
#' @export
weighted_probs_experts <- function(df_pred, fitted_experts, weights, levels = 1:5) {
  l <- length(levels)
  n <- nrow(df_pred)
  
  one_hot_experts <- outer(fitted_experts, levels, "==") + 0
  scores <- matrix(0, nrow = n, ncol = l)
  colnames(scores) <- levels  # Label columns for clarity
  for (i in seq_len(l)) {
    pi_ax <- one_hot_experts[, , i]
    scores[, i] <- as.vector(weights %*% t(pi_ax))
  }
  return(scores)
}

#' Convert a Binary Matrix to a Confidence Set List
#'
#' Converts a binary matrix (or any numeric matrix with non-zero elements) into a list 
#' where each element represents a row and contains the column indices of non-zero entries.
#'
#' @param binary_matrix A numeric or logical matrix where non-zero (or `TRUE`) entries 
#'   indicate membership or inclusion.
#'
#' @return A list of integer vectors of length `nrow(binary_matrix)`. The $i$-th element 
#'   contains the 1-based column indices where row $i$ has non-zero values. If a row has 
#'   no non-zero entries, an empty integer vector (`integer(0)`) is returned for that row.
#' @export
binary_to_confidence_set <- function(binary_matrix) {
  idx <- which(binary_matrix != 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(lapply(seq_len(nrow(binary_matrix)), function(i) integer(0)))
  }
  split(idx[, "col"], factor(idx[, "row"], levels = seq_len(nrow(binary_matrix))))
}

utils::globalVariables(c("val", "row_id", "exists"))
#' Create a standard data block
#'
#' Reshapes a 4D array into a long-format data frame.
#'
#' @param mech_index Integer index for the mechanism.
#' @param mech_name String label for the mechanism.
#' @param array_a 4D array with dimensions: (level, rep, mechanism, rate).
#' @param alpha_vec Vector of alpha levels.
#' @param rate_data Vector of rate labels.
#' @export
make_block <- function(mech_index, mech_name, array_a, alpha_vec, rate_data) {
  n_i <- dim(array_a)[4]   # random rate index
  n_a <- dim(array_a)[1]   # level index
  purrr::map_dfr(1:n_i, function(i) {
    purrr::map_dfr(1:n_a, function(a) {
      data.frame(
        value = array_a[a, , mech_index, i], mechanism = mech_name,
        level = paste0(alpha_vec[a]), type = paste0(rate_data[i])
      )
    })
  })
}

#' Create a simplified data block
#'
#' Reshapes a 3D array where the first dimension is already a composite or vector.
#'
#' @param mech_index Integer index for the mechanism.
#' @param mech_name String label for the mechanism.
#' @param array_a 3D array with dimensions: (value, mechanism, rate).
#' @param rate_data Vector of rate labels.
#' @export
make_smaller_block <- function(mech_index, mech_name, array_a, rate_data){
  purrr::map_dfr(1:length(array_a[1,mech_index,]), function(i) {
    data.frame(
      value = array_a[ , mech_index, i],
      mechanism = mech_name,
      type = paste0(rate_data[i])
    )
  })
}