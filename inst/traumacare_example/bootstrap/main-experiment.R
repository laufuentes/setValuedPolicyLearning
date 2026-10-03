root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning/"
setwd(root.path)
seed <- 2026
set.seed(seed)

names <- readRDS("inst/traumacare_example/intermediate/preprocessing.rds")


# ── Load functions from R folder  ────────────────────────────────────────────
source("inst/libraries.R")
source("R/utils.R")
source("R/evaluation.R")


source("inst/toy_examples/train_policies.R")

# ── General parameters  ───────────────────────────────────────────────────────
random_rate <- c(0, 0.1, 0.25, 0.5)
n_rate <- length(random_rate)
alpha <- 0.1
z <- qnorm(1 - alpha/2)
n_bootstrap <- 30

# ── Load data  ──────────────────────────────────────────────
### Train samples
train_imp_onehot <- read.csv("inst/traumacare_example/intermediate/train_imp_one_hot.csv")
df_obs <- train_imp_onehot # use train_imp (if want to use no-one-hot encoded version)
n<- nrow(df_obs)  # number of observations. 

### Test samples (for final predictions)
test_imp_onehot <- read.csv("inst/traumacare_example/intermediate/test_imp_one_hot.csv")
df_new_sample <- test_imp_onehot # use test_imp (if want to use no-one-hot encoded version)

# ── Define data parameters  ─────────────────────────────────────────────────
# Baseline covariates
covariates_name <- names$covariates_name # use covariate_name (if want to use no-one-hot encoded version ) 
X <- df_obs[, covariates_name] %>% 
  as.matrix() %>% 
  apply(2, as.numeric) # Covariates for training data 

X_new <- df_new_sample[,covariates_name] %>% 
  as.matrix() %>% 
  apply(2, as.numeric) # Covariates for test data 

# Treatment
treatment_name <- names$treatment_name
A <- df_obs[,treatment_name] # Treatment vector for training data 
A_new <- df_new_sample[,treatment_name] # Treatment for test data
levels_A <- levels(A) # treatment levels
m <- length(levels(A)) # number of treatment levels

# Outcome
outcome_name <- names$outcome_name
Y <- df_obs[,outcome_name] 
Y_new <- df_new_sample[,outcome_name]
ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new)))

set.seed(seed)
bootstrap_indices <- lapply(1:n_bootstrap, function(i) {
  sample(1:n, size = as.integer(n*0.75), replace = TRUE)})

source("inst/traumacare_example/bootstrap/set-valued-policy-training.R")
saveRDS(results_list, file = "inst/traumacare_example/images/results_list.rds")

#source("inst/traumacare_example/bootstrap/figures.R")