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
#' @param ncov Number of covariates (2 by default).
#' @param seed Integer or NA (NA by default).
#' @param type String indicating the type of synthetic scenario ("normal", "complex")
#' @param is_RCT Logical indicating if treatment allocated as in RCT (TRUE by default).
#'
#' @return A list containing two data frames (\code{df_obs} with observed outcomes 
#' based on treatment and \code{df_complete} with all potential outcomes) and the 
#' oracular optimal treatment assignments. 
#' @examples
#' n <- 1e3 
#' generate_data(n, type="normal")
#' @export
generate_data <- function(n, ncov=2, seed=NA, type=c("normal", "complex"), is_RCT= TRUE){
  ncov <- R.utils::Arguments$getIntegers(ncov, c(2, 15))
  type <- match.arg(type)
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  treatment_levels <- 5 
  X <- matrix(stats::rnorm(ncov*n), ncol=ncov)
  if(is_RCT){
    A <- t(stats::rmultinom(n, 1, rep(1/treatment_levels, treatment_levels)))
  }else{
    if(type=="normal"){
      beta_low_vec  <- c(10,10,1,3,2)  # 1,2 high
      beta_high_vec <- c(2,1,10,10,4)  # 3,4 low 
      w <- plogis(X[,1] + X[,2] - 0.5)
      beta <- (1 - w) * matrix(beta_low_vec,  nrow=nrow(X), ncol=5, byrow=TRUE) +
        w * matrix(beta_high_vec, nrow=nrow(X), ncol=5, byrow=TRUE)
      epsilon <- matrix(rnorm(nrow(X) * 5), nrow=nrow(X), ncol=5)
      
      treatment_assignment <- beta + epsilon
      
      probs <- exp(treatment_assignment - apply(treatment_assignment, 1, max))
      expit_treatment <- probs / rowSums(probs)
      
      A <- t(apply(expit_treatment, 1, function(p) rmultinom(1, 1, p)))
    }else{
      beta_high <- c(4, 10, 4, 2, 4)  # treatment 2 and then a bit of 1 and 4
      beta_medium <- c(10, 2, 2, 10, 15) # treatment 4, 1 and 5 and a bit of the others 
      beta_low <- c(4, 2, 10, 4, 2)  # treatment 3 and then a bit of 1 and 4
      z_axis <- X[,1] - X[,2]
      w_low  <- plogis(5*(z_axis - 0.25))
      w_high  <- plogis(5*(-z_axis - 0.5))
      w_medium <- pmax(0, 1 - (w_high + w_low))
      total_w  <- w_low + w_medium + w_high
      
      beta <- (w_low/total_w) * matrix(beta_low,  nrow=nrow(X), ncol=5, byrow=TRUE) +
        (w_high/total_w) * matrix(beta_high, nrow=nrow(X), ncol=5, byrow=TRUE) + 
        (w_medium/total_w) * matrix(beta_medium, nrow=nrow(X), ncol=5, byrow=TRUE)
      
      epsilon <- matrix(rnorm(nrow(X) * 5), nrow=nrow(X), ncol=5)
      treatment_assignment <- beta + epsilon
      probs <- exp(treatment_assignment - apply(treatment_assignment, 1, max))
      expit_treatment <- probs / rowSums(probs)
      A <- t(apply(expit_treatment, 1, function(p) rmultinom(1, 1, p)))
    }
  }
  A_int <- max.col(A)
  A_factor <- factor(A_int)
  stopifnot(all(rowSums(A) == 1))
 
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
