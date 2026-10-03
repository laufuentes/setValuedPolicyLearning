#' Baseline Outcome Effect Function
#'
#' Calculates the non-treatment related baseline component of the outcome model 
#' based on a linear combination of covariates.
#'
#' @param X A numeric matrix of covariates of size n x 2.
#'
#' @return A numeric vector of length n containing the baseline effects.
#' @examples
#' X <- matrix(stats::runif(10 * 2), 10, 2)
#' baseline_effect(X)
#' @export
baseline_effect <- function(X){
  beta <- c(1, -0.5, 1, -0.5, 1)
  return(X %*% beta)
}


#' Conditional Mean Outcome: Linear Scenario
#'
#' Simulates the conditional mean of the outcome under a "linear" scenario. 
#' It creates a linear boundary that partitions which 
#' treatment (1-2 vs. 3) receive an additive effect boost.
#'
#' @param X A numeric matrix of covariates of size n x 2.
#'
#' @return A numeric matrix of size n x 4, where each column represents the 
#' conditional mean for a specific treatment. 
#' 
#' @examples
#' X <- matrix(stats::runif(10 * 2), 10, 2)
#' mu_P_linear(X)
#' @export
mu_P_linear <- function(X){
  out <- matrix(baseline_effect(X), nrow = nrow(X), ncol = 4)
  cond <- (X[, 1] + X[, 2] < 0.5)
  out[cond, 1:2]  <- out[cond, 1] + 2
  out[!cond, 3] <- out[!cond, 3] + 3
  return(out)
}

#' Conditional Mean Outcome: Sinusoid Frontier Scenario
#'
#' Simulates the conditional mean of the outcome under a "complex" scenario. 
#' It creates a sinusoidal boundary that partitions which 
#' treatment (1-2 vs. 3) receive an additive effect boost.
#'
#' @param X A numeric matrix of covariates of size n x 5.
#'
#' @return A numeric matrix of size n x 4, where each column represents the conditional mean for a specific treatment. 
#' @examples
#' X <- matrix(stats::runif(10 * 2), 10, 2)
#' mu_P_complex(X)
#' @export
mu_P_complex <- function(X){
  out <- matrix(baseline_effect(X), nrow = nrow(X), ncol = 4)
  b <- 0.5 #2
  rho <- 1 #0.75
  xy <- cbind(X[,1],X[,2])
  cond <- as.matrix(xy) %*% c(1,1) <= 1 + b * sin(2*pi*rho*xy[,1])
  out[cond, 1:2]  <- out[cond, 1] + 2
  out[!cond, 3] <- out[!cond, 3] + 3
  return(out)
}

mu_P_tree <- function(X){
  out <- matrix(baseline_effect(X), nrow = nrow(X), ncol = 4)
  cond1 <- (X[,1] <= 0) & (X[,2] <= 0)
  out[cond1, 1:2] <- out[cond1, 1] + 3
  
  cond2 <- (X[,1] > 0) & (X[,2] > 0)
  out[cond2, 2:3] <- out[cond2, 2] +3
  
  cond3 <- (X[,1] <= 0) & (X[,2] > 0)
  out[cond3, 1] <- out[cond3, 1] + 2
  
  cond4 <- (X[,1] > 0) & (X[,2] <= 0)
  out[cond4, 3] <- out[cond4, 3] + 2
  return(out)
}



#' Synthetic data generator
#'
#' Generates a dataset simulating treatment assignment, covariates, and potential outcomes.
#'
#' @param n Number of observations to generate.
#' @param seed Integer or NA (NA by default).
#' @param is_RCT Logical indicating if treatment allocated as in RCT (TRUE by default).
#' @param type String indicating the type of synthetic scenario ("linear", "complex", "tree")
#'
#' @return A list containing two data frames (\code{df_obs} with observed outcomes 
#' based on treatment and \code{df_complete} with all potential outcomes) and the 
#' oracular optimal treatment assignments. 
#' @examples
#' n <- 1e3 
#' generate_data(n, type="simple")
#' @export
generate_data <- function(n, seed=NA, is_RCT= TRUE, type = c("linear", "complex", "tree")){
  type <- match.arg(type)
  if(!is.na(seed)){
    set.seed(seed)
  }
  ncov <- 5
  treatment_levels <- 4 
  X <- matrix(stats::rnorm(ncov*n), ncol=ncov)
  if(is_RCT){
    A <- t(stats::rmultinom(n, 1, rep(1/treatment_levels, treatment_levels)))
  }else{
    if(type=="linear"){
      w <- stats::plogis(X[,1] + X[,2] - 0.5)
      beta_low_vec  <- c(8,8,5,7)  # 1,2 high
      beta_high_vec <- c(4,4,8,7)  # 3 low 
      
      beta <- (1 - w) * matrix(beta_low_vec,  nrow=nrow(X), ncol=treatment_levels, byrow=TRUE) +
        w * matrix(beta_high_vec, nrow=nrow(X), ncol=treatment_levels, byrow=TRUE)
    }else if(type=="complex"){
      w <- stats::plogis(X[,1] + X[,2] - 1)
      beta_low_vec  <- c(10,10,5,7)  # 1,2 high
      beta_high_vec <- c(4,4,10,7)  # 3 low 
      
      beta <- (1 - w) * matrix(beta_low_vec,  nrow=nrow(X), ncol=treatment_levels, byrow=TRUE) +
        w * matrix(beta_high_vec, nrow=nrow(X), ncol=treatment_levels, byrow=TRUE)
    }else{
      s1_pos <- stats::plogis(X[, 1]) 
      s2_pos <- stats::plogis(X[, 2]) 
  
      w_Q1 <- (1 - s1_pos ) * (1 - s2_pos)  # X1 <= 0, X2 <= 0
      w_Q2 <- s1_pos * s2_pos  # X1 > 0,  X2 > 0
      w_Q3 <- (1 - s1_pos) * s2_pos  # X1 <= 0, X2 > 0
      w_Q4 <- s1_pos * (1 - s2_pos)  # X1 > 0,  X2 <= 0
      
      beta_Q1 <- c(10, 10,  5,  8)  # Q1 favors Treatments 1 & 2 
      beta_Q2 <- c( 5, 10, 10,  8)  # Q2 favors Treatments 2 & 3
      beta_Q3 <- c(10,  5,  5,  8)  # Q3 favors Treatment 1    
      beta_Q4 <- c( 5,  5, 10,  8)  # Q4 favors Treatment 3    
      
      beta <- w_Q1 %*% t(beta_Q1) + 
        w_Q2 %*% t(beta_Q2) + 
        w_Q3 %*% t(beta_Q3) + 
        w_Q4 %*% t(beta_Q4)
      
    }
    
      probs <- exp(beta - apply(beta, 1, max))
      expit_treatment <- probs / rowSums(probs)
      
      A <- t(apply(expit_treatment, 1, function(p) stats::rmultinom(1, 1, p)))
  }
  A_int <- max.col(A)
  A_factor <- factor(A_int)
  stopifnot(all(rowSums(A) == 1))
  
  if(type=="linear"){
    potential_outcomes <- mu_P_linear(X)
  }else if(type=="complex"){
    potential_outcomes <- mu_P_complex(X)
  }else{
    potential_outcomes <- mu_P_tree(X)
  }
  
  Y_obs <- rowSums(potential_outcomes*A) + 0.1*stats::rnorm(n, sd = 1)
  optimal_policy <- lapply(seq_len(nrow(potential_outcomes)),
                           function(i){
                             which(potential_outcomes[i, ] == 
                                     max(potential_outcomes[i, ]))})
  df <- data.frame(X,A = A_factor,Y = Y_obs)
  
  df_complete <- data.frame(
    X, A = A_factor, Y = Y_obs, 
    Potential_outcomes = potential_outcomes)
  
  if(is_RCT){
    prop_score <- matrix(0.5, nrow=n, ncol=treatment_levels)
  }else{
    prop_score <- expit_treatment
  }
  
  return(list(df, df_complete, optimal_policy, prop_score))
}
