#' Baseline Outcome Effect Function
#'
#' Calculates the non-treatment related baseline component of the outcome model 
#' based on a linear combination and an exponential transformation of covariates.
#'
#' @param X A numeric matrix of covariates of size n x 2.
#'
#' @return A numeric vector of length n containing the baseline effects.
#' @examples
#' X <- matrix(stats::runif(10 * 2), 10, 2)
#' baseline_effect(X)
#' @export
baseline_effect <- function(X){
  return(2*X[,1] - exp(X[,2] ))
}

#' Conditional Mean Outcome: Linear Frontier Scenario
#'
#' Simulates the conditional mean of the outcome under a "normal" scenario. 
#' It creates a linear boundary that partitions which 
#' treatment (1-2 vs. 3-4) receive an additive effect boost.
#'
#' @param X A numeric matrix of covariates of size n x 2.
#'
#' @return A numeric matrix of size n x 5, where each column represents the conditional mean for a specific treatment. 
#' @examples
#' X <- matrix(stats::runif(10 * 2), 10, 2)
#' mu_P0_normal(X)
#' @export
mu_P0_normal <- function(X){
  out <- matrix(1, nrow=nrow(X), ncol=5)*baseline_effect(X)
  cond <- (X[,1]+X[,2]<0.5)
  out[cond, 1:2] <- out[cond, 1:2] + exp(X[cond,1])
  out[!cond, 3:4] <- out[!cond, 3:4] + 2*(X[!cond,2])^2
  return(out)
}

#' Conditional Mean Outcome: Complex Simplex Scenario
#'
#' Defines optimal treatments across five distinct regions of the covariate space.
#' This scenario includes a circular boundary and quadrant-based logic for the 
#' remaining space, incorporating stochastic noise in the fifth treatment arm.
#'
#' @param X A numeric matrix of covariates of size n x 2.
#'
#' @return A numeric matrix of size n x 5 containing the transformed conditional means for each treatment. 
#' @examples
#' X <- matrix(stats::runif(10 * 2), 10, 2)
#' mu_P0_simplex_complicated(X)
#' @export
mu_P0_simplex_complicated <- function(X){
  out <- matrix(1, nrow=nrow(X), ncol=5)*baseline_effect(X)
  
  cond1 <- (X[,1]^2 + X[,2]^2) < 0.5
  out[cond1, ] <- out[cond1, ] + 5
  out[cond1, 5] <- out[cond1, 5] + 0.25*stats::rbinom(n=length(out[cond1, 5]), 
                                                      size=1, prob=0.5)
  
  cond2 <- (!cond1) & (X[,1] > 0) & (X[,2] > 0)
  out[cond2, 1] <- out[cond2, 1] + 4*abs(X[cond2,1])
  
  cond3 <- (!cond1) & (X[,1] <= 0) & (X[,2] > 0)
  out[cond3, 2] <- out[cond3, 2] + 4
  
  cond4 <- (!cond1) & (X[,1] > 0) & (X[,2] <= 0)
  out[cond4, 3] <- out[cond4, 3] + 5
  
  cond5 <- (!cond1) & (X[,1] <= 0) & (X[,2] <= 0)
  out[cond5, 4] <- out[cond5, 4] + exp(X[cond5, 2])
  return(out)
}


#' Synthetic data generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param seed Integer or NA (NA by default).
#' @param type String indicating the type of synthetic scenario ("normal", "complex")
#'
#' @return A list containing two data frames (\code{df_obs} with observed outcomes 
#' based on treatment and \code{df_complete} with all potential outcomes) and the 
#' oracular optimal treatment assignments. 
#' @examples
#' n <- 1e3 
#' generate_data(n)
#' @export
generate_data <- function(n, seed=NA, type=c("normal", "complex")){
  if(!is.na(seed)){
    set.seed(seed)
  }
  m <- 5
  A <- t(stats::rmultinom(n, 1, rep(1/m, m)))
  A_int <- max.col(A)        # integer 1..m
  A_factor <- factor(A_int)
  stopifnot(all(rowSums(A) == 1))
  X <- matrix(stats::rnorm(2*n), ncol=2)
  if(type=="complex"){
    potential_outcomes <- mu_P0_simplex_complicated(X)
  }else{
    potential_outcomes <- mu_P0_normal(X)
  }
   Y_obs <- rowSums(potential_outcomes*A) + rnorm(n, sd = 1)
   optimal_policy <- lapply(seq_len(nrow(potential_outcomes)),
                            function(i){
                              which(potential_outcomes[i, ] == max(potential_outcomes[i, ]))})
   df <- data.frame(
     x1 = X[,1],
     x2 = X[,2],
     A = A_factor,
     Y = Y_obs
   )
   
   df_complete <- data.frame(
     x1 = X[,1],
     x2 = X[,2],
     A = A_factor,
     Y = Y_obs, 
     Potential_outcomes = potential_outcomes)
   return(list(df, df_complete, optimal_policy))
}
