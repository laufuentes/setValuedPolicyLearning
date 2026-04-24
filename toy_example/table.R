types_optimal_treatment <- SL.out$optimal_policy_new %>% unique()
obs_types <- sapply(SL.out$optimal_policy_new, function(x) {
  which(sapply(types_optimal_treatment, function(y) setequal(x, y)))
})


train1 <- SL.out$df_obs[SL.out$folds[[1]],] # generate noisy labels 
train2 <-  SL.out$df_obs[SL.out$folds[[2]],] # score model and nuisances
test <-  SL.out$df_obs[SL.out$folds[[3]],]


alpha <- 0.1
# levels of noise
r_levels <- c(1,3,6,11)
cov_unif<- mean_width <- matrix(0, ncol = length(r_levels)+2)
confidence_sets <- list()

# oracular 
true_quant <- stats::quantile(SL.out$true_score, (1-alpha))
true_binary_confidence_set <-  ifelse(SL.out$new_scores<true_quant, 1, 0)
true_idx <- which(true_binary_confidence_set  != 0, arr.ind = TRUE) 
confidence_sets[[length(confidence_sets)+1]] <- split(true_idx[, "col"], 
                                                      factor(true_idx[, "row"], 
                                                             levels = seq_len(nrow(true_binary_confidence_set))))
cov_unif[,length(confidence_sets)] <- coverage_unif(true_set = SL.out$optimal_policy_new, confidence_sets[[length(confidence_sets)]])
mean_width[,length(confidence_sets)] <- width(confidence_sets[[length(confidence_sets)]])

# noisy ICP
for (r in r_levels){
  quant <- stats::quantile(SL.out$rate_scores_unweighted_cal[, r], (1-alpha))
  binary_confidence_set <-  ifelse(SL.out$new_scores<quant, 1, 0)
  idx <- which(binary_confidence_set  != 0, arr.ind = TRUE) 
  confidence_sets[[length(confidence_sets)+1]] <- split(idx[, "col"], 
                           factor(idx[, "row"], 
                                  levels = seq_len(nrow(binary_confidence_set))))
  cov_unif[,length(confidence_sets)] <- coverage_unif(true_set = SL.out$optimal_policy_new, confidence_sets[[length(confidence_sets)]])
  mean_width[,length(confidence_sets)] <- width(confidence_sets[[length(confidence_sets)]])
  
}
# uppest lower bound set
treatment<- matrix(0, nrow=nrow(rbind(train1, train2)), ncol=m)
treatment[cbind(1:nrow(treatment), c(train1[, treatment_name],
                                     train2[, treatment_name]))] <- 1
model <- grf::regression_forest(
  X = cbind(rbind(train1, train2)[, covariates_name], treatment),
  Y = rbind(train1, train2)[, outcome_name])


lowers <- uppers <- lowers_test <- uppers_test <- matrix(0, nrow=nrow(SL.out$df_new), ncol=m)
z <- stats::qnorm(1 - alpha/2)
for (l in as.numeric(levels_A)){
    treatment_l <- matrix(0, nrow=nrow(SL.out$df_new), ncol=m)
    treatment_l[,l] <- 1
    data_l <- data.frame(SL.out$df_new[,covariates_name], Treatment=treatment_l)
    pred <- stats::predict(model, newdata = data_l, estimate.variance = TRUE)
    se <- sqrt(pred$variance.estimates)
    lowers[,l] <- pred$predictions - z * se
    uppers[,l] <- pred$predictions + z * se
  }
uppest_lrw_bound <- apply(lowers, 1, max)
C_set_binary_naive <- ifelse(uppers>=uppest_lrw_bound, 1, 0)
indices_naive <- which(C_set_binary_naive != 0, arr.ind = TRUE)
confidence_sets[[length(confidence_sets)+1]] <- split(indices_naive[, "col"], indices_naive[, "row"])
cov_unif[,length(confidence_sets)] <- coverage_unif(true_set = SL.out$optimal_policy_new, confidence_sets[[length(confidence_sets)]])
mean_width[,length(confidence_sets)] <- width(confidence_sets[[length(confidence_sets)]])

prop_inclusion_type_treatment <- matrix(0, 
                                        nrow=length(types_optimal_treatment), 
                                        ncol=length(confidence_sets))
for(c in 1:length(confidence_sets)) {
  
  all_coverage <- mapply(coverage_strict_single, 
                         SL.out$optimal_policy_new, 
                         confidence_sets[[c]])
  
  # Group by type and calculate the mean (proportion)
  # Ensure obs_types matches the rows of your matrix
  prop_inclusion_type_treatment[, c] <- tapply(all_coverage, obs_types, mean)
}

# Create data frame
treatment_labels <- sapply(types_optimal_treatment, function(x) {
  paste0("\\{", paste(x, collapse = ", "), "\\}")
})

library(xtable)
results_df <- as.data.frame(prop_inclusion_type_treatment)
colnames(results_df) <- c("Oracular CP", 
                          apply(data.frame(1:length(r_levels)), 
                                1, 
                                function(i)paste0("$r=$",random_rate[r_levels[i]])),
                          "GLB")
results_df <- cbind("Optimal Treatments" = treatment_labels, results_df)  
total_row <- c(paste0("Coverage ($\\alpha$ = ", alpha, ")"), cov_unif)
width_row <- c("Mean cardinality", mean_width)

results_df <- rbind(results_df, total_row, width_row)
results_df[,-1] <- lapply(results_df[,-1], function(x) as.numeric(as.character(x)))

# 2. Apply the 1.2f format to all numeric columns
results_df[,-1] <- lapply(results_df[,-1], function(x) sprintf("%1.2f", x))

addline <- list()
addline$pos <- list(nrow(results_df) - 2, nrow(results_df)-1)
addline$command <- rep("\\hline \n", 2)

# 3. Print the LaTeX Code
print(xtable(results_df,
             align =  paste0("ll", paste(rep("c", ncol(results_df) - 1), 
                                                     collapse = ""))), 
      include.rownames = FALSE, 
      hline.after = c(-1, 0, nrow(results_df)), 
      add.to.row = addline, 
      sanitize.text.function = function(x) x, 
      file = paste0("images/", n, "/", "Table_coverage_", type, ".txt")) 

