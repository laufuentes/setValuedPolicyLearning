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

#' Create a high-dimensional data block
#'
#' Reshapes a 5D array into a long-format data frame including score and replication indices.
#'
#' @param mech_index Integer index for the mechanism.
#' @param mech_name String label for the mechanism.
#' @param array_a 5D array with dimensions: (level, rep, mechanism, score, rate).
#' @param alpha_vec Vector of alpha levels.
#' @param rate_data Vector of rate labels.
#' @param score_vec Vector of score labels.
make_bigger_block <- function(mech_index, mech_name, array_a, alpha_vec, rate_data, score_vec) {
  n_a      <- dim(array_a)[1]
  n_rep    <- dim(array_a)[2]
  n_score  <- dim(array_a)[4]
  n_rate   <- dim(array_a)[5]

  purrr::map_dfr(1:n_score, function(s) {
    purrr::map_dfr(1:n_rate, function(i) {
      purrr::map_dfr(1:n_a, function(a) {

        data.frame(
          value = array_a[a, , mech_index, s, i],   # length = n_rep
          rep = seq_len(n_rep),
          mechanism = mech_name,
          level = alpha_vec[a],                      # numeric
          type = factor(rate_data[i]),               # factor for discrete color
          score = score_vec[s]
        )

      })
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

#' Create a 4D data block with scores
#'
#' Reshapes a 4D array specifically ordered for level, mechanism, score, and rate.
#'
#' @param mech_index Integer index for the mechanism.
#' @param mech_name String label for the mechanism.
#' @param array_a 4D array with dimensions: (level, mechanism, score, rate).
#' @param alpha_vec Vector of alpha levels.
#' @param rate_data Vector of rate labels.
#' @param score_vec Vector of score labels.
#' @export
make_smaller_bigger_block <- function(mech_index, mech_name, array_a, alpha_vec, rate_data, score_vec) {
  n_rate  <- dim(array_a)[4]
  n_alpha <- dim(array_a)[1]
  n_score <- dim(array_a)[3]

  purrr::map_dfr(1:n_rate, function(i) {
    purrr::map_dfr(1:n_score, function(s) {
      data.frame(
        value = array_a[, mech_index, s, i],
        mechanism = mech_name,
        level = alpha_vec,
        type = factor(rate_data[i]),
        score = score_vec[s]
      )
    })
  })
}

#' Create a heatmap for treatments included in set-valued policies
#'
#' Creates a data frame one-hot encoding the inclusion of treatments
#' in set-valued policies. Rows represent observations and columns
#' the treatment levels.
#'
#' @param confidence_set Integer index for the mechanism.
#' @param levels_A String label for the mechanism.
#' @export
heatmap_treatments <- function(confidence_set, levels_A){
  `%>%`<- magrittr::`%>%`
  df <- tibble::tibble(
    row_id = 1:length(confidence_set),
    val = confidence_set
  ) %>%
    tidyr::unnest(val, keep_empty = TRUE) %>%
    dplyr::mutate(
      val = factor(val, levels = levels_A),
      exists = 1
    ) %>%
    tidyr::pivot_wider(
      names_from = val,
      values_from = exists,
      values_fill = 0,
      names_expand = TRUE
    ) %>% dplyr::select(dplyr::any_of(levels_A))
  return(df)
}

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

binary_to_confidence_set <- function(binary_matrix) {
  idx <- which(binary_matrix != 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(lapply(seq_len(nrow(binary_matrix)), function(i) integer(0)))
  }
  split(idx[, "col"], factor(idx[, "row"], levels = seq_len(nrow(binary_matrix))))
}
