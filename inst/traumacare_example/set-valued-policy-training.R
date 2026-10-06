root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning/"
setwd(root.path)
seed <- 2026
set.seed(seed)

names <- readRDS("inst/traumacare_example/intermediate/preprocessing.rds")


# ── Load functions from R folder  ────────────────────────────────────────────
source("inst/libraries.R")
source("R/utils.R")
source("R/evaluation.R")
source("inst/traumacare_example/train_policies.R")

# ── General parameters  ───────────────────────────────────────────────────────
random_rate <- c(0, 0.1, 0.25, 0.5)
n_rate <- length(random_rate)
alpha <- 0.1
z <- qnorm(1 - alpha/2)

# ── Load data  ────────────────────────────────────────────────────────────────
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
X <- df_obs[, covariates_name] |> 
  as.matrix() |> 
  apply(2, as.numeric) # Covariates for training data 

X_new <- df_new_sample[,covariates_name] |> 
  as.matrix() |> 
  apply(2, as.numeric) # Covariates for test data 

# Treatment
treatment_name <- names$treatment_name
A <- df_obs[,treatment_name] |> as.factor() # Treatment vector for training data 
levels_A <- levels(A) # treatment levels
m <- length(levels(A)) # number of treatment levels

# Outcome
outcome_name <- names$outcome_name
Y <- df_obs[,outcome_name] 
ab <- c(min(Y),max(Y))

# ── Split data into 3 folds  ────────────────────────────────────────────────
# Fold 1: Set-valued policy learning
# Fold 2: Predict set-valued policies (r selection)
# Fold 3: Train nuisances to evaluate on Fold2
set.seed(seed)
shuffled_indices <- sample(1:nrow(df_obs))
cut1 <- round(0.5 * nrow(df_obs))
cut2 <- round(0.2 * nrow(df_obs))
custom_folds <- list(
  Fold1 = shuffled_indices[1:cut1],
  Fold2 = shuffled_indices[(cut1 + 1): (cut1 + cut2)],
  Fold3 = shuffled_indices[(cut1 + cut2 + 1):nrow(df_obs)])
  
folds <- SuperLearner::CVFolds(nrow(df_obs), id = NULL, Y = df_obs[,outcome_name] ,
                                        cvControl = SuperLearner::SuperLearner.CV.control(V = 3L,
                                                                                          validRows = custom_folds))
  
train <- df_obs[folds[[1]],] # train set-valued policies 
n_train <- nrow(train)
pseudo.test.predict <-  df_obs[folds[[2]],] # predict set-valued policies 
evaluation <-  df_obs[folds[[3]],] # train nuisances for set-valued policy evaluation
  
# Train nuisances for evaluation
 QAW.reg.train.r = grf::probability_forest(
 X = cbind(evaluation[,covariates_name],
               evaluation[,treatment_name]), 
 Y = evaluation[,outcome_name] |> as.factor())
   
 Q.all.pseudo.r <-  do.call(cbind,lapply(0:1, function(val) {
     new_data <- cbind(pseudo.test.predict[,covariates_name], val)
     stats::predict(QAW.reg.train.r, newdata = new_data)$predictions[,2]}))
  
g.reg.train.r <- grf::probability_forest(X = evaluation[,covariates_name],
                                           Y = evaluation[,treatment_name]|> 
                                             as.factor())

gAW.pred.pseudo.r <- stats::predict(g.reg.train.r, 
                                        newdata = pseudo.test.predict[, covariates_name])$predictions
  
# ── SET_VALUED POLICY LEARNING ────────────────────────────────────────────────

# ── 1. CONFORMAL SET_VALUED POLICY LEARNING ───────────────────────────────────
# ── 0) Divide data into three even sets ───────────────────────────────────────
# ── Noisy label generation, scoring model & calibration ───────────────────────
set.seed(seed)
shuffled_indices <- sample(1:n_train)
cut1 <- cut2 <- round(0.4 * n_train)
custom_folds_conformal <- list(
    Fold1 = shuffled_indices[1:cut1],
    Fold2 = shuffled_indices[(cut1 + 1): (cut1 + cut2)],
    Fold3 = shuffled_indices[(cut1 + cut2 + 1):n_train])
  
folds_conformal <- SuperLearner::CVFolds(n_train, id = NULL, 
                                           Y = train[,outcome_name],
                                           cvControl = SuperLearner::SuperLearner.CV.control(V = 3L, validRows = custom_folds_conformal))
  
train1 <- train[folds_conformal[[1]],] # generate noisy labels
train2 <-  train[folds_conformal[[2]],] # score model and nuisances
calibration <-  train[folds_conformal[[3]],] # calibration
  
# ── 1) Black-box label generation (i.e. estimates of (X,A*)) ──────────────────
## 1.1) Generate random labels (i.e. A_rd)
A_rd <- sample(as.numeric(levels_A), size = nrow(calibration), replace = TRUE)
  
## 1.2) Estimate A* (OTR) using experts
# Training performed on train1
SL.library_cond <- c("SL.randomForest", "SL.mean", "SL.gam", "SL.glm", 
                     "SL.xgboost")
trained_policies_results <- train_policies(train_b = df_obs, train1 = train1, 
                 calibration = calibration, pseudo_test = pseudo.test.predict,
                 seed = seed, SL.library_cond = SL.library_cond) # GLB trained inside
  
saveRDS(trained_policies_results,
          file="inst/traumacare_example/images/trained_policies_results.rds")

# Extract results from training
doptFactorPredict_cal <- trained_policies_results$doptFactorPredict_cal 
doptFactorPredict_pseudo <- trained_policies_results$doptFactorPredict_pseudo
doptFactorPredict_new_naive <- trained_policies_results$doptFactorPredict_new_naive 
doptFactorPredict_pseudo_naive <- trained_policies_results$doptFactorPredict_pseudo_naive
model.glb.pf <- trained_policies_results$model.glb.pf 
model.glb.glm <- trained_policies_results$model.glb.glm
unweighted.naive_new <- trained_policies_results$unweighted.naive
unweighted.pseudo.naive <- trained_policies_results$unweighted.pseudo.naive
numalgs <- ncol(doptFactorPredict_cal) # number of baselines
  
selected_method <- "ql.SL" 

# Unweighted aggregation of baselines
unweighted_probs <- weighted_probs_experts(fitted_experts = doptFactorPredict_cal,
                                             weights =rep(1/numalgs, numalgs),
                                             df_pred = calibration,
                                             levels = as.numeric(levels_A))
  
unweighted_aggregation <- apply(apply(unweighted_probs, 1, function(x){
    rmultinom(1,1,prob=x)}), 2, which.max)-1 
  

# Extract baseline policies for conformal procedure
policy_cal <- doptFactorPredict_cal[,selected_method]
  
# Randomness injection
policy_cal_r <- matrix(0, nrow = nrow(calibration), ncol = n_rate)
for (i in seq_len(n_rate)){
    r <- random_rate[i]
    mix_factor<- stats::rbinom(nrow(calibration), 1, prob=r) # R ~ Ber(r)
    policy_cal_r[,i] <- mix_factor * A_rd + (1 - mix_factor) * policy_cal
}
  
selected_methods <- c(lapply(random_rate, function(r) {
                          paste0(selected_method, " (r=", r, ")")
                        }) |> unlist())
# ── 2) Train nonconformity score model (i.e. s(X,A)) ──────────────────────────
# Training performed on train2
# Two predictions:
# (i) on calibration
# (ii) on df_new_sample
 QAW.reg.train_conformal = grf::probability_forest(
     X = cbind(train2[,c(covariates_name,treatment_name)]), 
     Y = train2[,outcome_name] |> as.factor())

# ── 3) Calibration step  ──────────────────────────────────────────────────────
 potential_outcomes_cal <- do.call(cbind,lapply(0:1, function(val) {
     new_data <- calibration[,c(covariates_name,treatment_name)]
     new_data[,treatment_name] <- val 
     stats::predict(QAW.reg.train_conformal, newdata = new_data)$predictions[,2]}))

# Compute margin score on calibration data 
margin_po <-  margin_score(potential_outcomes_cal)

# Extract scores for different label types
r0_scores_policy <- apply(policy_cal_r,2,function(x){
  margin_po[cbind(1:nrow(calibration), x+1)]})  # Baseline policies and perturbed labels 
r0_scores_aggregation <- margin_po[cbind(1:nrow(calibration), 
                                           unweighted_aggregation+1)] # aggregation of baselines
r1_score <- margin_po[cbind(seq_len(nrow(calibration)), A_rd+1)] # random labels 

ecdf_data <- apply(cbind(doptFactorPredict_cal, 
                         unweighted=unweighted_aggregation), 2, 
                   function(x)margin_po[cbind(1:nrow(calibration),x+1)])|> 
  as.data.frame() |>
  pivot_longer(cols=everything(), names_to = "methods", values_to = "value")|> 
  ggplot(aes(x = value, colour = methods)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  geom_hline(yintercept = 1-alpha, colour = "red") +
  stat_ecdf(data= data.frame(r1_score = r1_score),
            aes(x = r1_score, colour = "Random labels"),
            linetype = "dashed", colour="gray",
            linewidth = 1.2
  ) +
  labs(y = "ECDF", x = "Value", colour = "Method")
ggplot2::ggsave(ecdf_data, 
                filename = paste0("inst/traumacare_example/images/ecdf.pdf"), 
                width = 10, height = 8)

 potential_outcomes_pseudo <- do.call(cbind,lapply(levels_A, function(val) {
     new_data <- pseudo.test.predict[, c(covariates_name, treatment_name)]
     new_data[,treatment_name] <- as.numeric(val)
     stats::predict(QAW.reg.train_conformal, newdata = new_data)$pred[,2]}))

# Compute margin score on pseudo test data (for r selection)
margin_po_pseudo <-  margin_score(potential_outcomes_pseudo)
  
 potential_outcomes_new <-  do.call(cbind, lapply(levels_A, function(val) {
   new_data <- df_new_sample[, covariates_name]
   new_data[,treatment_name] <- as.numeric(val)
   stats::predict(QAW.reg.train_conformal, newdata = new_data)$pred[,2]}))
  
# Compute margin score on new data (for final prediction)
margin_po_new <-  margin_score(potential_outcomes_new)
  
# Build and evaluate conformal set-valued policies   ─────────────────────────
pred_pseudo <- list()

# Baseline policy-based noisy labels
i<- 1
results_policy <- apply(r0_scores_policy, 2, function(x){
  quant <- stats::quantile(x, (1-alpha))
  conf_12 <- binary_to_confidence_set(margin_po_pseudo < quant)
  conf <- lapply(conf_12,function(x)x-1)
  i<- i+1
  pred_pseudo[[selected_methods[i]]] <- conf
  table.evaluation.real(conf, prop_score_new = gAW.pred.pseudo.r, 
                        potential_outcomes = Q.all.pseudo.r,
                        df_new_sample = pseudo.test.predict,
                        levels_A = levels_A, zero_indexed = TRUE,
                        treatment_name = treatment_name, ab = ab, 
                        outcome_name = outcome_name)
  }) 
names(results_policy) <- selected_methods 


# Aggregation-based noisy labels
quant_agg <- quantile(r0_scores_aggregation, 1 - alpha)
conf_set_agg_12 <- binary_to_confidence_set(margin_po_new < quant_agg)
conf_set_agg <- lapply(conf_set_agg_12,function(x)x-1)
pred_pseudo[["conf.agg"]] <- conf_set_agg
results_agg <- table.evaluation.real(conf_set_agg, 
                                     prop_score_new = gAW.pred.pseudo.r, 
                                     potential_outcomes = Q.all.pseudo.r,
                                     df_new_sample = pseudo.test.predict,
                                     levels_A = levels_A, ab = ab, 
                                     zero_indexed = TRUE,
                                     treatment_name = treatment_name, 
                                     outcome_name = outcome_name)
  
# ── 2.GREATEST LOWER BOUND (GLB) ──────────────────────────────────────────────
# ── Using regression forest for estimation ────────────────────────────────────
lowers <- uppers <- matrix(0, nrow=nrow(pseudo.test.predict), ncol=m)
for (l in as.numeric(levels_A)){
    data_l <- data.frame(pseudo.test.predict[,c(covariates_name, 
                                                treatment_name)])
    data_l[,treatment_name] <- l
    pred <- stats::predict(model.glb.pf, 
                           newdata = data_l, 
                           estimate.variance = TRUE)
    se <- sqrt(pred$variance.estimates[,2])
    lowers[,l+1] <- pred$predictions[,2] - z * se
    uppers[,l+1] <- pred$predictions[,2] + z * se}
uppest_lrw_bound <- apply(lowers, 1, max)
conf_set_pf_12 <- binary_to_confidence_set(uppers >= uppest_lrw_bound)
conf_set_pf <-  lapply(conf_set_pf_12,function(x)x-1)
pred_pseudo[["glb.grf"]] <- conf_set_pf
results_glb_pf <- table.evaluation.real(conf_set_pf, 
                                         prop_score_new = gAW.pred.pseudo.r, 
                                         potential_outcomes = Q.all.pseudo.r,
                                         df_new_sample = pseudo.test.predict,
                                         levels_A = levels_A, ab = ab,
                                         treatment_name = treatment_name, 
                                         outcome_name = outcome_name,
                                         zero_indexed = TRUE)
  
# ── Using linear model with interactions for estimation ─────────────────────
lowers <- uppers <- matrix(0, nrow=nrow(pseudo.test.predict), ncol=m)
for (l in as.numeric(levels_A)){
    data_l <- pseudo.test.predict[,c(covariates_name,treatment_name)]
    data_l[,treatment_name] <- l
    pred <- stats::predict(model.glb.glm, newdata = data_l, se.fit = TRUE, type = "response")
    se <- pred$se.fit 
    lowers[,l+1] <- (pred$fit - z * se) |> as.numeric()
    uppers[,l+1] <- (pred$fit + z * se) |> as.numeric()}
uppest_lrw_bound <- apply(lowers, 1, max)
conf_set_glm_12 <- binary_to_confidence_set(uppers>=uppest_lrw_bound)
conf_set_glm <- lapply(conf_set_glm_12, function(x)x-1)
pred_pseudo[["glb.glm"]] <- conf_set_glm
results_glb_glm <- table.evaluation.real(conf_set_glm, 
                                        prop_score_new = gAW.pred.pseudo.r, 
                                        potential_outcomes = Q.all.pseudo.r,
                                        df_new_sample = pseudo.test.predict,
                                        levels_A = levels_A, ab = ab,
                                        treatment_name = treatment_name, 
                                        outcome_name = outcome_name, 
                                        zero_indexed = TRUE)
  
# Save confidence intervals 
saveRDS(pred_pseudo, file = "inst/traumacare_example/images/conf.rds")

# ── 3.Create figures for the evaluation fold ──────────────────────────────────
dynamic_methods <- unlist(lapply(names(results_policy), function(m) {
  paste0("Conformal ", m)})) 
names(results_policy) <- dynamic_methods
all_results <- c(results_policy, 
                 list("Conformal aggregation" = results_agg), 
                 list("GLB PF" = results_glb_pf),
                 list("GLB GLM" = results_glb_glm))

# Save results
saveRDS(all_results, 
        file = "inst/traumacare_example/images/all_results.rds")

order_elements <- c("GLB PF", "GLB GLM", names(results_policy),
                    "Conformal aggregation")

# Cardinality plot 
cardinality <- bind_rows(lapply(all_results, `[[`, 1), 
                         .id = "Set-valued policy") 

cardinality_plot <- cardinality|> 
  select(all_of(order_elements)) |>
  pivot_longer(cols = everything(), names_to = "Method", 
               values_to = "Values") |> 
  mutate(Method = factor(Method, levels = order_elements)) |> 
  ggplot(aes(x = Method, y = Values)) +
  geom_point(shape = 3) + 
  labs(x = "Set-valued policies", y = "Mean cardinality" )+
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(cardinality_plot, 
       file = "inst/traumacare_example/images/cardinality_pseudo.pdf", 
       width = 10, height = 6)

# SPV plot 
spv_data <- lapply(all_results, `[[`, 3)|>
  bind_rows(.id = "Set-valued policy")

doctors <- mean(pseudo.test.predict[,outcome_name])
naive_baseline <- Q.all.pseudo.r[cbind(1:nrow(pseudo.test.predict), 
                                       unweighted.pseudo.naive+1)] |> mean()
plot_spv <- spv_data |> 
  filter(policy=="Propensity")|>
  mutate(`Set-valued policy` = factor(`Set-valued policy`, 
                                      levels = order_elements)) |> 
  ggplot( aes(x = `Set-valued policy`, y = value, 
              ymin = lower, ymax = upper, 
              color = `Set-valued policy`)) +
  geom_errorbar() + 
  geom_hline(yintercept = doctors, 
             color="black", 
             linetype = "dashed")+
  geom_hline(yintercept = naive_baseline, color="red",
             linetype = "dashed")+
  facet_grid(~estimator)+ 
  ylim(c(0.90,1))+
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none")+
  labs(x = "Set-valued policies", y = "SPV")

ggsave(plot_spv,
       filename = "inst/traumacare_example/images/SPV_boxplots_pseudo.pdf", 
       width = 10, height = 5)

# Latex table 
table_results <- bind_rows(lapply(all_results, `[[`, 2), .id = "Set-valued policy") |>
  mutate(across(-1, ~ sprintf("%.3f", coalesce(as.numeric(.x), 0)))) |>
  rename_with(~ paste0("Prop. of \\{", .x, "\\}"), -1)

spv_clean <- spv_data |>
  mutate(
    Metric = sprintf("%s SPV %s", estimator, tolower(policy)),
    Val = sprintf("[%.3f, %.3f]", lower, upper)
  ) |>
  select(`Set-valued policy` = 1, Metric, Val) |>
  pivot_wider(names_from = Metric, values_from = Val)


cardinality_clean <- cardinality |>
  pivot_longer(everything(), 
               names_to = "Set-valued policy", 
               values_to = "Mean card.") |>
  mutate(`Mean card.` = sprintf("%.3f", 
                                      coalesce(as.numeric(`Mean card.`), 0)))


joined_df <- table_results |>
  left_join(cardinality_clean, by = "Set-valued policy") |>
  left_join(spv_clean, by = "Set-valued policy") |>
  t()

colnames(joined_df) <- joined_df[1, ]
final_data <- as.data.frame(joined_df[-1, , drop = FALSE]) |>
  select(all_of(order_elements)) |> 
  tibble::rownames_to_column(var = "Metric")

sub_headers <- c("\\textbf{Set Type}", "PF", "GLM", 
                 paste0("r=", random_rate), "r=0")

header_groups <- c(
  " " = 1,
  "\\\\textbf{GLB}" = 2,
  "\\\\textbf{Conformal ql.SL}" = 4,
  "\\\\textbf{Conformal Agg.}" = 1)


latex_tbl <- kable(
  final_data,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  col.names = sub_headers,
  align = c("l", rep("c", ncol(final_data) - 1)),
  label = "tab:traumacare") |>
  add_header_above(header_groups, escape = FALSE) |> 
  kable_styling(latex_options = c("HOLD_position"), font_size = 9)


writeLines(
  c(
    "% latex table generated in R",
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    "\\setlength{\\tabcolsep}{1pt}",
    as.character(latex_tbl),
    "\\end{table}"
  ),
  con = "inst/traumacare_example/images/table.txt")

# ── 4. Final prediction  ──────────────────────────────────────────────────────
# Train nuisances for SPV plug-in estimator
QAW.reg.train.all = grf::probability_forest(
  X = cbind(X,A), Y = Y |> as.factor())

Q.all <-  do.call(cbind,lapply(0:1, function(val) {
  new_data <- cbind(df_new_sample[,covariates_name], val)
  stats::predict(QAW.reg.train.all, newdata = new_data)$predictions[,2]}))

g.reg.train.all <- grf::probability_forest(X = X,Y = A|> as.factor())

gAW.pred.all <- stats::predict(g.reg.train.all, 
                               newdata = df_new_sample[, covariates_name])$predictions

# Final predictions
# GLB GLM 
lowers <- uppers <- matrix(0, nrow=nrow(df_new_sample), ncol=m)
for (l in as.numeric(levels_A)){
  data_l <- df_new_sample[,covariates_name]
  data_l[,treatment_name] <- l
  pred <- stats::predict(model.glb.glm, newdata = data_l, se.fit = TRUE, type = "response")
  se <- pred$se.fit 
  lowers[,l+1] <- (pred$fit - z * se) |> as.numeric()
  uppers[,l+1] <- (pred$fit + z * se) |> as.numeric()}
uppest_lrw_bound <- apply(lowers, 1, max)
conf_set_glm_12_final <- binary_to_confidence_set(uppers>=uppest_lrw_bound)
conf_set_glm_final <- lapply(conf_set_glm_12_final, function(x)x-1)

spv_final.glb <- set_policy_value_plug_in(conf_set_glm_final, 
                                      test = df_new_sample,
                                      Q.all.actions = Q.all,
                                      gAX.pred = gAW.pred.all,
                                      zero_indexed = TRUE,
                                      levels = levels_A)

tab.glb <- table(sapply(conf_set_glm_final, function(el){
  paste0("{", paste0(el, collapse = ", "), "}")}))/length(conf_set_glm_final)

cardinality.mean.glb <- sapply(1:length(conf_set_glm_final), 
                           function(i) {
                             length(conf_set_glm_final[[i]])|> as.numeric()
                             }) |> mean()

# Conformal r=0
quant <- stats::quantile(r0_scores_policy[,which(random_rate==0)], (1-alpha))
conf_12_final <- binary_to_confidence_set(margin_po_new < quant)
conf_final <- lapply(conf_12_final,function(x)x-1)

spv_final.conf <- set_policy_value_plug_in(conf_final, 
                                      test = df_new_sample,
                                      Q.all.actions = Q.all,
                                      gAX.pred = gAW.pred.all,
                                      zero_indexed = TRUE,
                                      levels = levels_A)

tab.conf <- table(sapply(conf_final, function(el){
  paste0("{", paste0(el, collapse = ", "), "}")}))/length(conf_set_glm_final)

cardinality.mean.conf <- sapply(1:length(conf_final), 
                           function(i) {
                             length(conf_set_glm_final[[i]])|> as.numeric()
                           }) |> mean()

# Create a small summary table
tab_main <- cbind(
  "GLB GLM"       = tab.glb, 
  "Conformal r=0" = tab.conf)

rownames(tab_main) <- paste0("Prop. of ", rownames(tab_main))

results_final <- rbind(
  tab_main,
  "Propensity SPV" = c(spv_final.glb[[2]], spv_final.conf[[2]]))

little_latex_table <- kable(results_final, format = "latex", booktabs = TRUE, digits = 3)
writeLines(
  c(
    "% latex table generated in R",
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    "\\setlength{\\tabcolsep}{1pt}",
    as.character(little_latex_table),
    "\\end{table}"
  ),
  con = "inst/traumacare_example/images/table_little.txt")

