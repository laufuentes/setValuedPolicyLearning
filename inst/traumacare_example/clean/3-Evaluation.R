set.seed(seed)
#SL.out <- readRDS("inst/traumacare_example/intermediate/SL.out.rds")
cat("evaluating for varying alpha...")

# General parameters ───────────────────────────────────────────────────────────
alphas <- c(0.01, seq(0.05,0.5, 0.05))
n_test <- 1e2 # for evaluation

spv <- array(0, dim=c(length(alphas), n_test, 3, 2))
mean_cardinality<- array(0, dim=c(length(alphas), 3))

# GLB 
pred.glb.levels <-  se.glb.levels <- list()
for (l in as.numeric(levels_A)){
  new_data <- cbind(pseudo.test.predict[,covariates_name], l)
  pred <- stats::predict(SL.out$model.glb, newdata = new_data, estimate.variance = TRUE, 
                         type= "response")
  pred.glb.levels[[l+1]] <- pred$predictions[,2]
  print(head(pred$predictions[,2]))
  se.glb.levels[[l+1]] <- sqrt(pred$variance.estimates[,2])
  print(head())
}

# Test varying alpha 
for(i in 1:length(alphas)){
  alpha <- alphas[i]
  
  # Conformal policy learning (unweighted)
  quant <- stats::quantile(SL.out$scores.density.unweighted, (1-alpha))
  binary_confidence_set <-  ifelse(SL.out$pseudo_scores<quant, 1, 0)
  idx <- which(binary_confidence_set  != 0, arr.ind = TRUE)
  confidence_set.unweighted <- split(idx[, "col"]-1, 
                                     factor(idx[, "row"], levels = seq_len(nrow(binary_confidence_set))))
  
  spv_res <- set_policy_value(confidence_set.unweighted, 
                              test= pseudo.test.predict,
                              levels=levels_A, n_test = n_test,
                              treatment_name = treatment_name,
                              outcome_name = outcome_name,
                              covariates = covariates_name, 
                              Q_all_actions = Q.all.pseudo,
                              gAW.pred = gAW.pred.pseudo, 
                              gAW.pred.spv = gAW.pred.pseudo, ab = ab)
  
  spv[i,,1,1] <- spv_res[[1]] # random 
  spv[i,,1,2] <- spv_res[[2]] # clinicians
  mean_cardinality[i,1]<- width(pred_set = confidence_set.unweighted)
  
  # Conformal policy learning (single policy)
  quant <- stats::quantile(SL.out$scores.density.single, (1-alpha))
  binary_confidence_set <-  ifelse(SL.out$pseudo_scores<quant, 1, 0)
  idx <- which(binary_confidence_set  != 0, arr.ind = TRUE)
  confidence_set.single <- split(idx[, "col"]-1,
                          factor(idx[, "row"],
                                 levels = seq_len(nrow(binary_confidence_set))))
  
  spv_res.single <- set_policy_value(confidence_set.single, test= pseudo.test.predict,
                              levels=levels_A, n_test = n_test,
                              treatment_name = treatment_name,
                              outcome_name = outcome_name,
                              covariates = covariates_name, 
                              Q_all_actions = Q.all.pseudo,
                              gAW.pred = gAW.pred.pseudo, 
                              gAW.pred.spv = gAW.pred.pseudo, ab = ab)
  
  spv[i,,2,1] <- spv_res.single[[1]] # random 
  spv[i,,2,2] <- spv_res.single[[2]] # clinicians
  mean_cardinality[i,2]<- width(pred_set = confidence_set.single)
  # GLB 
  lowers <- uppers <- matrix(0,nrow=nrow(pseudo.test.predict), ncol=m)
  z <- stats::qnorm(1 - alpha/2)
  for (l in as.numeric(levels_A)){
    lowers[,l+1] <- pred.glb.levels[[l+1]] - z * se.glb.levels[[l+1]]
    uppers[,l+1] <- pred.glb.levels[[l+1]] + z * se.glb.levels[[l+1]]
  }
  uppest_lrw_bound <- apply(lowers, 1, max)
  C_set_binary_naive <- ifelse(uppers>=uppest_lrw_bound, 1, 0)
  indices_naive <- which(C_set_binary_naive != 0, arr.ind = TRUE)
  naive.confidence_set <- split(indices_naive[, "col"]-1,
                                factor(indices_naive[, "row"],
                                       levels = seq_len(nrow(C_set_binary_naive))))
  
  spv_res.glb <- set_policy_value(naive.confidence_set, test= pseudo.test.predict,
                                  levels=levels_A, n_test = n_test,
                                  treatment_name = treatment_name,
                                  outcome_name = outcome_name,
                                  covariates = covariates_name, 
                                  Q_all_actions = Q.all.pseudo,
                                  gAW.pred = gAW.pred.pseudo, 
                                  gAW.pred.spv = gAW.pred.pseudo, ab = ab)
  
  spv[i,,3,1] <- spv_res.glb[[1]] # random 
  spv[i,,3,2] <- spv_res.glb[[2]] # clinicians
  mean_cardinality[i,3]<- width(pred_set = naive.confidence_set)
}

results <- list(spv=spv,mean_cardinality=mean_cardinality)
saveRDS(results, file = "inst/traumacare_example/intermediate/results_pseudotest.rds")

dimnames(spv) <- list(
  alpha = alphas,
  element = 1:n_test,
  method = c("Unweighted", "Single.policy", "GLB"),
  metric = c("uniform", "propensity"))


data_SPV_plot <- as.data.frame.table(spv, responseName = "value") %>%
  mutate(
    value = as.numeric(as.character(value)),
    element = as.numeric(as.character(element))
  )

doctors <- mean(pseudo.test.predict[,outcome_name])

plot_spv <- ggplot2::ggplot(data = data_SPV_plot, 
                            ggplot2::aes(x=alpha, 
                                         y=value, color=method))+
  ggplot2::geom_boxplot()+
  ggplot2::geom_hline(yintercept = doctors, 
                      show.legend = TRUE, color="black") +
  ggplot2::labs(title = "Set-policy values",
                y="Set-policy value") +
  ggplot2::facet_grid(~metric)+
  ggplot2::theme()

ggplot2::ggsave(plot = plot_spv, 
                filename = paste0("inst/traumacare_example/images/SPV_",outcome_name, "_seed", seed,".pdf"), 
                width = 8, height = 5)


dimnames(mean_cardinality) <- list(
  alpha = alphas,
  method = c("Unweighted", "Single.policy", "GLB"))

data_cardinality_plot <- as.data.frame.table(mean_cardinality, responseName = "value")
plot_cardinality <- ggplot2::ggplot(data = data_cardinality_plot,
                                    ggplot2::aes(x = alpha, y=value,
                                                 group = method,
                                                 colour = method))+
  ggplot2::geom_point(shape=4, size=3.5)+
  ggplot2::labs(title = "Mean cardinalities",
                y="Mean cardinality") +
  ggplot2::ylim(c(0,2)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                     vjust = 1,
                                                     hjust = 1))

ggplot2::ggsave(plot = plot_cardinality, 
                filename = paste0("inst/traumacare_example/images/Cardinalities_",outcome_name, "_seed", seed,".pdf"), 
                width = 8, height = 5)


# 3) Train the nuisances for final evaluation ────────────────────────────────────────
# Final evaluation 
## Outcome model (Q-model)  
QAW.reg.train.final= grf::probability_forest(
  X = cbind(X,A), 
  Y = Y %>% as.factor())

## Propensity score (G-model)
g.reg.train.final <- grf::probability_forest(X = X,
                                             Y = A %>% as.factor())


