# ── Set working directory  ──────────────────────────────────────────────────
root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning"
setwd(root.path)

# ── Required packages  ────────────────────────────────────────────────────────
source("inst/libraries.R")

# ── Load functions from R folder  ────────────────────────────────────────────
source("inst/toy_example_simple/synthetic_data.R")
source("R/utils.R")
source("R/evaluation.R")
source("inst/toy_example_simple/train_policies.R")

# ── General parameters  ───────────────────────────────────────────────────────
seed <- 2026
set.seed(seed)
VFolds <- 4 # folds to split data

n <- 10000
type <- "normal"
#random_rates <- seq(0,1,0.1)
alpha <- 0.1
z <- qnorm(1 - alpha/2)
n_bootstrap <- 30

# ── Simulations for varying sample sizes  ─────────────────────────────────────
SL.out<- list() # list where results will be saved

# ── Synthetic data generation  ──────────────────────────────────────────────
## Training observations
exp <- generate_data(n, is_RCT = FALSE, seed = seed, type = type)
# extract observational data
SL.out$df_obs <- exp[[1]]
# extract complete data
df_complete <- exp[[2]]
summary(df_complete)

# extract optimal policy
SL.out$optimal_policy <- exp[[3]]
# extract potential outcomes
SL.out$potential_outcomes <- df_complete %>%
  select(starts_with("Potential_outcomes."))

### Test observations
exp_new_sample <- generate_data(n/2, is_RCT = FALSE, seed = seed+1, type = type)
# extract observational data
SL.out$df_new_sample <- exp_new_sample[[1]]
# extract optimal policy
SL.out$optimal_policy_new <- exp_new_sample[[3]]
# extract potential outcomes
SL.out$potential_outcomes <- exp_new_sample[[2]] %>%
  select(starts_with("Potential_outcomes."))
SL.out$prop_score_new <- exp_new_sample[[4]] 

# ── Define data parameters  ─────────────────────────────────────────────────
covariates_name <- c("X1","X2", "X3", "X4", "X5")
X <- SL.out$df_obs[,covariates_name] %>% as.matrix()
X_new <- SL.out$df_new_sample[,covariates_name] %>% as.matrix()

treatment_name <- "A" # name of treatment indicator in dataset
A <- SL.out$df_obs[,treatment_name]
A_new <- SL.out$df_new_sample[,treatment_name]

levels_A <- levels(A) # treatment levels
m <- length(levels(A)) # number of treatment levels

outcome_name <- "Y" # name of outcome in dataset
Y <- SL.out$df_obs[,outcome_name]
Y_new <- SL.out$df_new_sample[,outcome_name]

ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new)))

set.seed(seed)
bootstrap_indices <- lapply(1:n_bootstrap, function(i) {
  sample(1:n, size = as.integer(n*0.75), replace = TRUE)})

source("inst/toy_example_simple/set-valued-policy-training.R")
saveRDS(results_list, file = "inst/toy_example_simple/images/results_list.rds")
saveRDS(SL.out, file = "inst/toy_example_simple/images/SL.out.rds")


source("inst/toy_example_simple/table.R")

# plot ECDF of nonconformity scores
# data_toghether <- cbind(r0_scores_policy,r0_scores_aggregation) %>%  
#   as.data.frame() %>% 
#   pivot_longer(cols = everything(), 
#                names_to = "Method",
#                values_to = "Value")
# 
# ggplot(data_toghether, aes(x = Value, colour = Method)) +
#   stat_ecdf(geom = "step", linewidth = 1, alpha=0.75) +
#   geom_hline(yintercept = 1-alpha, colour = "red") +
#   stat_ecdf(
#     data = as.data.frame(true_score),
#     aes(x = true_score,colour = "Oracular labels"),
#     linetype = "dashed", colour="black",
#     linewidth = 1.2
#   ) +
#   stat_ecdf(
#     data = as.data.frame(r1_score),
#     aes(x = r1_score, colour = "Random labels"),
#     linetype = "dashed", colour="gray",
#     linewidth = 1.2
#   ) +
#   labs(y = "ECDF", x = "Value", colour = "Method")
# 
# ggplot2::ggsave(filename = paste0("inst/toy_example_simple/images/ecdf_", n,".pdf"), 
#                 width = 10, height = 8)

