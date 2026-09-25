root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning"
setwd(root.path)
source("inst/toy_example_simple/synthetic_data.R")
library(tidyr)
library(dplyr)

seed<- 2026
n <- 10000
exp <- generate_data(n=n, is_RCT = FALSE, seed = seed, type = "normal")
df_obs <- exp[[1]]
df_complete <- exp[[2]]

summary(df_complete)

df_complete_plot <- df_complete %>%
  pivot_longer(
    cols = starts_with("Potential_outcomes."),
    names_to = "potential.outcome",
    values_to = "value"
  )

ggplot2::ggplot(data = df_complete_plot, 
                ggplot2::aes(x=value,color=potential.outcome))+ 
  ggplot2::geom_density(alpha=0.5)

X <- df_obs %>% 
  select(starts_with("X")) %>% 
  as.matrix()

A <- df_obs$A
levels_A <- levels(A)
m <- length(levels_A)

Y <- df_obs$Y

covariates_name <- colnames(df_obs)[1:5]

mod_glm <- stats::lm(formula = Y ~(X1 + X2 + X3 + X4 + X5)* A,
                     data = df_obs[1:5000,])

potential_outcomes_test_glm <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- df_obs[5001:n,] %>% select(-Y)
  new_data$A <- factor(val, levels=levels_A)
  stats::predict(mod_glm, newdata = new_data)}))

mod_grf <- grf::regression_forest(
  X = cbind(X[1:5000,], A[1:5000] %>% as.factor()), 
  Y = Y[1:5000], seed = seed)

potential_outcomes_test_grf <- do.call(cbind,lapply(1:m, function(val) {
  new_data <- cbind(X[5001:n,], factor(val, levels=levels_A) %>% as.numeric())
  stats::predict(mod_grf, newdata = new_data)$predictions}))

potential_outcome_estimations <- array(0, dim=c(nrow(potential_outcomes_test_grf), 4, 3))
potential_outcome_estimations[,,1] <- df_complete[5001:n,] %>% select(starts_with("Pot")) %>% 
  as.matrix()
potential_outcome_estimations[,,2] <- potential_outcomes_test_grf
potential_outcome_estimations[,,3] <- potential_outcomes_test_glm

dimnames(potential_outcome_estimations)<- list(
  rows = 1:nrow(potential_outcomes_test_grf), 
  potential.outcomes = 1:4, 
  methods = c("ground.truth", "regression.forest", "lm"))

df_plot <- as.data.frame.table(potential_outcome_estimations, responseName = "value") %>% 
  pivot_wider(names_from = methods, values_from = value)

ggplot2::ggplot(df_plot, ggplot2::aes(x = ground.truth)) +
  ggplot2::geom_point(ggplot2::aes(y = regression.forest, color = "Regression Forest"), 
                      alpha = 0.4, size = 1.5) +
  ggplot2::geom_point(ggplot2::aes(y = lm, color = "Linear Model"), 
                      alpha = 0.4, size = 1.5) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  ggplot2::facet_grid(rows = ggplot2::vars(potential.outcomes)) +
  ggplot2::scale_color_manual(
    values = c("Regression Forest" = "red", "Linear Model" = "blue")
  ) +
  ggplot2::labs(
    title = "Ground Truth vs. P.O. Estimations",
    x = "Ground Truth",
    y = "Estimated P.O",
    color = "Estimation Method"
  ) +
  ggplot2::theme_minimal()

######################################################
#################  propensity score  #################
######################################################
beta_low_vec  <- c(10,10,5,7)  # 1,2 high
beta_high_vec <- c(4,4,10,7)  # 3 low 
w <- stats::plogis(X[,3] + X[,4] - 0.5)
beta <- (1 - w) * matrix(beta_low_vec,  nrow=nrow(X), ncol=treatment_levels, byrow=TRUE) +
  w * matrix(beta_high_vec, nrow=nrow(X), ncol=treatment_levels, byrow=TRUE)

probs <- exp(beta - apply(beta, 1, max))
expit_treatment <- probs / rowSums(probs)

plot_expit_treatment <- expit_treatment %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), 
               names_to = "treatment", 
               values_to = "value")

ggplot2::ggplot(data = plot_expit_treatment,
                ggplot2::aes(x = value, color=treatment))+ 
  ggplot2::geom_density(alpha=0.2)

