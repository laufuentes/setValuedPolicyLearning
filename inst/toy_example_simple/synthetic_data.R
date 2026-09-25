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
  beta <- c(1, -0.5, 1, -0.5, 1)
  return(X %*% beta)
}

#TODO: essayer un effet de base plus compliqué
# a <- 0.1
# 
# rho <- 3
# 
# x <- y <- seq(0,1,length.out=1e3)
# xy <- expand.grid(x,y)
# z <- as.matrix(xy) %*% c(1,1) <= 1 + a * sin(2*pi*rho*xy[, 1])
# bind_cols(xy, z) |> as_tibble() |> rename(x=Var1, y=Var2, z=`...3`) |>
#   ggplot2::ggplot() + ggplot2::geom_raster(ggplot2::aes(x=x,y=y,fill=z))

#' Conditional Mean Outcome: Linear Frontier Scenario
#'
#' Simulates the conditional mean of the outcome under a "normal" scenario. 
#' It creates a linear boundary that partitions which 
#' treatment (1-2 vs. 3) receive an additive effect boost.
#'
#' @param X A numeric matrix of covariates of size n x 2.
#'
#' @return A numeric matrix of size n x 4, where each column represents the conditional mean for a specific treatment. 
#' @examples
#' X <- matrix(stats::runif(10 * 2), 10, 2)
#' mu_P0_normal(X)
#' @export
mu_P0_normal <- function(X){
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
#' mu_P0_complex(X)
#' @export
mu_P0_complex <- function(X){
  out <- matrix(baseline_effect(X), nrow = nrow(X), ncol = 4)
  b <- 0.5 #2
  rho <- 1 #0.75
  xy <- cbind(X[,1],X[,2])
  cond <- as.matrix(xy) %*% c(1,1) <= 1 + b * sin(2*pi*rho*xy[,1])
  out[cond, 1:2]  <- out[cond, 1] + 2
  out[!cond, 3] <- out[!cond, 3] + 3
  return(out)
}

mu_P0_tree <- function(X){
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
#' @param ncov Number of covariates (4 by default).
#' @param seed Integer or NA (NA by default).
#' @param is_RCT Logical indicating if treatment allocated as in RCT (TRUE by default).
#'
#' @return A list containing two data frames (\code{df_obs} with observed outcomes 
#' based on treatment and \code{df_complete} with all potential outcomes) and the 
#' oracular optimal treatment assignments. 
#' @examples
#' n <- 1e3 
#' generate_data(n, type="normal")
#' @export
generate_data <- function(n, seed=NA, is_RCT= TRUE, type = "normal"){
  #ncov <- R.utils::Arguments$getIntegers(ncov, c(2, 15))
  if(!is.na(seed)){
    set.seed(seed)
  }
  ncov <- 5
  treatment_levels <- 4 
  X <- matrix(stats::rnorm(ncov*n), ncol=ncov)
  if(is_RCT){
    A <- t(stats::rmultinom(n, 1, rep(1/treatment_levels, treatment_levels)))
  }else{
    if(type=="normal"){
      w <- stats::plogis(X[,1] + X[,2] - 0.5)
      beta_low_vec  <- c(10,10,5,7)  # 1,2 high
      beta_high_vec <- c(4,4,10,7)  # 3 low 
      
      beta <- (1 - w) * matrix(beta_low_vec,  nrow=nrow(X), ncol=treatment_levels, byrow=TRUE) +
        w * matrix(beta_high_vec, nrow=nrow(X), ncol=treatment_levels, byrow=TRUE)
    }else if(type=="complex"){
      w <- stats::plogis(X[,1] + X[,2] - 1)
      beta_low_vec  <- c(10,10,5,7)  # 1,2 high
      beta_high_vec <- c(4,4,10,7)  # 3 low 
      
      beta <- (1 - w) * matrix(beta_low_vec,  nrow=nrow(X), ncol=treatment_levels, byrow=TRUE) +
        w * matrix(beta_high_vec, nrow=nrow(X), ncol=treatment_levels, byrow=TRUE)
    }else{
      beta_high <- c(4, 10, 10, 7)  # treatment 2 and 3 
      beta_medium <- c(10, 10, 10, 7) # treatment 1, 2 and 3 
      beta_low <- c(10, 10, 4, 7)  # treatments 1,2
      z_axis <- X[,1] - X[,2]
      w_low  <- stats::plogis(5*(z_axis - 0.25))
      w_high  <- stats::plogis(5*(-z_axis - 0.5))
      w_medium <- pmax(0, 1 - (w_high + w_low))
      total_w  <- w_low + w_medium + w_high
      
      beta <- (w_low/total_w) * matrix(beta_low,  nrow=nrow(X), 
                                       ncol=treatment_levels, byrow=TRUE) +
        (w_high/total_w) * matrix(beta_high, nrow=nrow(X), 
                                  ncol=treatment_levels, byrow=TRUE) + 
        (w_medium/total_w) * matrix(beta_medium, nrow=nrow(X), 
                                    ncol=treatment_levels, byrow=TRUE)
      
    }
    
      probs <- exp(beta - apply(beta, 1, max))
      expit_treatment <- probs / rowSums(probs)
      
      A <- t(apply(expit_treatment, 1, function(p) stats::rmultinom(1, 1, p)))
  }
  A_int <- max.col(A)
  A_factor <- factor(A_int)
  stopifnot(all(rowSums(A) == 1))
  
  if(type=="normal"){
    potential_outcomes <- mu_P0_normal(X)
  }else if(type=="complex"){
    potential_outcomes <- mu_P0_complex(X)
  }else{
    potential_outcomes <- mu_P0_tree(X)
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
