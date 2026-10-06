# ── Set working directory  ──────────────────────────────────────────────────
root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning"
setwd(root.path)

# ── Required packages  ────────────────────────────────────────────────────────
source("inst/libraries.R")

# ── Load functions from R folder  ────────────────────────────────────────────
source("R/synthetic_data.R")
source("R/utils.R")
source("R/evaluation.R")
source("inst/toy_examples/1.2-train_policies.R")

# ── General parameters  ───────────────────────────────────────────────────────
seed <- 2026
set.seed(seed)
VFolds <- 4 # folds to split data

n <- 10000
type <- "tree"

subdir_path <- file.path(root.path, "inst", 
                         "toy_examples", paste0("images_", type))

if (!dir.exists(subdir_path)) {
  message(sprintf("Creating subdirectory: %s", "images"))
  dir.create(subdir_path, recursive = TRUE)
}

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
SL.out$potential_outcomes_train <- df_complete |>
  select(starts_with("Potential_outcomes."))

### Test observations
exp_new_sample <- generate_data(n/2, is_RCT = FALSE, seed = seed+1, type = type)
# extract observational data
SL.out$df_new_sample <- exp_new_sample[[1]]
# extract optimal policy
SL.out$optimal_policy_new <- exp_new_sample[[3]]
# extract potential outcomes
SL.out$potential_outcomes <- exp_new_sample[[2]] |>
  select(starts_with("Potential_outcomes."))
SL.out$prop_score_new <- exp_new_sample[[4]] 

# ── Define data parameters  ─────────────────────────────────────────────────
# Baseline covariates
covariates_name <- c("X1","X2", "X3", "X4", "X5")
X <- SL.out$df_obs[,covariates_name] |> as.matrix()
X_new <- SL.out$df_new_sample[,covariates_name] |> as.matrix()

# Treatment
treatment_name <- "A"
A <- SL.out$df_obs[,treatment_name]
A_new <- SL.out$df_new_sample[,treatment_name]
levels_A <- levels(A) # treatment levels
m <- length(levels(A)) # number of treatment levels

# Outcome
outcome_name <- "Y" 
Y <- SL.out$df_obs[,outcome_name]
Y_new <- SL.out$df_new_sample[,outcome_name]
ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new)))

set.seed(seed)
bootstrap_indices <- lapply(1:n_bootstrap, function(i) {
  sample(1:n, size = as.integer(n*0.75), replace = TRUE)})

# Implement set-valued policy learning across n_bootstrap samples
source("inst/toy_examples/1.1-set-valued-policy-training.R")

# Save results
saveRDS(results_list, file = file.path(subdir_path, "results_list.rds"))
saveRDS(SL.out, file = file.path(subdir_path, "SL.out.rds"))

source("inst/toy_examples/1.3-figures.R") # Create figures