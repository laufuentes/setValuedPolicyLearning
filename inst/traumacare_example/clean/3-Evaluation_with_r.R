set.seed(seed)
#SL.out <- readRDS("inst/traumacare_example/intermediate/SL.out.rds")
cat("evaluating for varying r...")

# General parameters ───────────────────────────────────────────────────────────
n_methods <- length(SL.out$perturbed_noisy_labels)+2

z <- stats::qnorm(1 - alpha/2)


spv <- array(0, dim=c(n_methods, 2, 3,  n_rate)) 
mean_cardinality<- array(0, dim=c(n_methods, n_rate))  

i<- 0
for(method in names(SL.out$perturbed_noisy_labels)){
  cal_score <- SL.out$calibration_scores[[method]] 
  i <- i+1
  for(r in seq_len(length(random_rate))){
    # Conformal policy learning (unweighted)
    quant <- stats::quantile(cal_score[,r], (1-alpha))
    binary_confidence_set <-  ifelse(SL.out$pseudo_scores<quant, 1, 0)
    idx <- which(binary_confidence_set  != 0, arr.ind = TRUE)
    confidence_set <- split(idx[, "col"]-1, 
                                       factor(idx[, "row"], levels = seq_len(nrow(binary_confidence_set))))
    
    spv_res.grf <- set_policy_value_plug_in(confidence_set, 
                                            test= pseudo.test.predict,
                                            levels=levels_A,
                                            Q.all.actions = Q.all.pseudo.grf,
                                            gAW.pred = gAW.pred.pseudo.grf, 
                                            zero_indexed=TRUE)
    
    spv[i, 1, 1, r] <- spv_res.grf[[1]] # random 
    spv[i, 2, 1, r] <- spv_res.grf[[2]] # clinicians
    
    spv_res.tmle <- set_policy_value_tmle(confidence_set, 
                                          Q.all.actions = Q.all.pseudo.grf, 
                                          gAW.pred = gAW.pred.pseudo.grf, 
                                          Y=pseudo.test.predict[,outcome_name], 
                                          A = pseudo.test.predict[,treatment_name],
                                          ab= ab, levels = levels_A, 
                                          zero_indexed = TRUE)
     
     spv[i, 1, 2, r] <- spv_res.tmle[[1]] # random 
     spv[i, 2, 2, r] <- spv_res.tmle[[2]] # clinicians
     
     
     spv_res.aipw <-set_policy_value_aipw(confidence_set, 
                           Q.all.actions = Q.all.pseudo.grf, 
                           gAW.pred = gAW.pred.pseudo.grf, 
                           Y=pseudo.test.predict[,outcome_name], 
                           A = pseudo.test.predict[,treatment_name], levels = levels_A, 
                           zero_indexed = TRUE)
     
     spv[i, 1, 3, r] <- spv_res.aipw[[1]] # random 
     spv[i, 2, 3, r] <- spv_res.aipw[[2]] # clinicians
     
    mean_cardinality[i, r]<- width(pred_set = confidence_set)
  } 
}
# GLB PF
lowers <- uppers <- matrix(0,nrow=nrow(pseudo.test.predict), ncol=m)
for (l in as.numeric(levels_A)){
  new_data <- pseudo.test.predict[,covariates_name]
  new_data[,treatment_name] <- l
  pred <- stats::predict(SL.out$model.glb.pf, newdata = new_data, estimate.variance = TRUE, 
                         type= "response")
    lowers[,l+1] <- pmax(pred$predictions[,2] - z * sqrt(pred$variance.estimates[,2]),0)
    uppers[,l+1] <- pmin(pred$predictions[,2] + z * sqrt(pred$variance.estimates[,2]),1)
}
uppest_lrw_bound <- apply(lowers, 1, max)
C_set_binary_naive <- ifelse(uppers>=uppest_lrw_bound, 1, 0)
indices_naive <- which(C_set_binary_naive != 0, arr.ind = TRUE)
naive.confidence_set.grf <- split(indices_naive[, "col"]-1,
                                factor(indices_naive[, "row"],
                                       levels = seq_len(nrow(C_set_binary_naive))))
  
spv_res.glb.grf <- set_policy_value_plug_in(naive.confidence_set.grf, test= pseudo.test.predict,
                                    levels=levels_A, Q.all.actions = Q.all.pseudo.grf, 
                                    gAW.pred = gAW.pred.pseudo.grf)
  
spv[n_methods-1, 1, 1, ] <- spv_res.glb.grf[[1]] # random 
spv[n_methods-1, 2, 1, ] <- spv_res.glb.grf[[2]] # clinicians
  

spv_res.tmle.grf <- set_policy_value_tmle(naive.confidence_set.grf, 
                                      Q.all.actions = Q.all.pseudo.grf, 
                                      gAW.pred = gAW.pred.pseudo.grf, 
                                      Y=pseudo.test.predict[,outcome_name], 
                                      A = pseudo.test.predict[,treatment_name],
                                      ab= ab, levels = levels_A, 
                                      zero_indexed = TRUE)

spv[n_methods-1, 1, 2, ] <- spv_res.tmle.grf[[1]] # random 
spv[n_methods-1, 2, 2,] <- spv_res.tmle.grf[[2]] # clinicians

spv_res.aipw.grf <-set_policy_value_aipw(naive.confidence_set.grf, 
                                     Q.all.actions = Q.all.pseudo.grf, 
                                     gAW.pred = gAW.pred.pseudo.grf, 
                                     Y=pseudo.test.predict[,outcome_name], 
                                     A = pseudo.test.predict[,treatment_name], levels = levels_A, 
                                     zero_indexed = TRUE)

spv[n_methods-1, 1, 3, ] <- spv_res.aipw.grf[[1]] # random 
spv[n_methods-1, 2, 3,] <- spv_res.aipw.grf[[2]] # clinicians

mean_cardinality[n_methods-1,]<- width(pred_set = naive.confidence_set.grf)

# GLB GLM
lowers <- uppers <- matrix(0, nrow=nrow(pseudo.test.predict), ncol=m)
for (l in as.numeric(levels_A)){
  data_l <- data.frame(pseudo.test.predict[,c(covariates_name,treatment_name)])
  data_l[,treatment_name] <- l
  pred <- stats::predict(SL.out$model.glb.glm, newdata = data_l, se.fit = TRUE)
  se <- sqrt(pred$se.fit)
  lowers[,l] <- (pred$fit - z * se) %>% as.numeric()
  uppers[,l] <- (pred$fit + z * se) %>% as.numeric()
}  
uppest_lrw_bound <- apply(lowers, 1, max)
C_set_binary_naive <- ifelse(uppers>=uppest_lrw_bound, 1, 0)
indices_naive <- which(C_set_binary_naive != 0, arr.ind = TRUE)
naive.confidence_set.glm <- split(indices_naive[, "col"]-1,
                              factor(indices_naive[, "row"],
                                     levels = seq_len(nrow(C_set_binary_naive))))

spv_res.glb.glm <- set_policy_value_plug_in(naive.confidence_set.glm, 
                                            test= pseudo.test.predict,
                                            levels=levels_A, 
                                            Q.all.actions =  Q.all.pseudo.grf,
                                            gAW.pred = gAW.pred.pseudo.grf)

spv[n_methods, 1, 1,] <- spv_res.glb.glm[[1]] # random 
spv[n_methods, 2, 1, ] <- spv_res.glb.glm[[2]] # clinicians

spv_res.tmle.glm <- set_policy_value_tmle(naive.confidence_set.glm, 
                                          Q.all.actions = Q.all.pseudo.grf, 
                                          gAW.pred = gAW.pred.pseudo.grf, 
                                          Y=pseudo.test.predict[,outcome_name], 
                                          A = pseudo.test.predict[,treatment_name],
                                          ab= ab, levels = levels_A, 
                                          zero_indexed = TRUE)

spv[n_methods, 1, 2,] <- spv_res.tmle.glm[[1]] # random 
spv[n_methods, 2, 2,] <- spv_res.tmle.glm[[2]] # clinicians

spv_res.aipw.glm <-set_policy_value_aipw(naive.confidence_set.glm, 
                                         Q.all.actions = Q.all.pseudo.grf, 
                                         gAW.pred = gAW.pred.pseudo.grf, 
                                         Y=pseudo.test.predict[,outcome_name], 
                                         A = pseudo.test.predict[,treatment_name], levels = levels_A, 
                                         zero_indexed = TRUE)

spv[n_methods, 1, 3, ] <- spv_res.aipw.glm[[1]] # random 
spv[n_methods, 2, 3,] <- spv_res.aipw.glm[[2]] # clinicians


mean_cardinality[n_methods,]<- width(pred_set = naive.confidence_set.glm)

results <- list(spv=spv,mean_cardinality=mean_cardinality)
saveRDS(results, file = "inst/traumacare_example/intermediate/results_pseudotest_r.rds")

dimnames(spv) <- list(
  #alpha = alphas,
  #element = 1:n_test,
  method = c(names(SL.out$calibration_scores), "GLB GRF", "GLB GLM"),
  metric = c("Uniform", "Propensity"), 
  technique =  c("plug-in", "TMLE", "AIPW"), 
  r_rates = random_rate)

data_SPV_plot <- as.data.frame.table(spv, responseName = "value") %>%
  mutate(
    value   = as.numeric(value),
  )

doctors <- mean(pseudo.test.predict[,outcome_name])
naive_method <- mean(Q.all.pseudo.grf[cbind(1:nrow(Q.all.pseudo.grf),
                                            SL.out$unweighted.pseudo.naive+1)])

plot_spv <- ggplot2::ggplot(data = data_SPV_plot, 
                            ggplot2::aes(x=method, 
                                         y=value, color= r_rates))+
  ggplot2::geom_point(data = ~ dplyr::filter(.x, !grepl("^GLB", method)),
                      ggplot2::aes(color = r_rates), shape = 3) +
  ggplot2::geom_point(data = ~ dplyr::filter(.x, grepl("^GLB", method)),
                      color = "black", shape = 3)+
  ggplot2::geom_hline(yintercept = doctors, 
                      show.legend = TRUE, color="black") +
  ggplot2::geom_hline(yintercept = naive_method, 
                      show.legend = TRUE, color="grey", linetype = "dashed") +
  ggplot2::labs(title = "Set-policy values",
                y="Set-policy value") +
  ggplot2::facet_grid(technique~metric)+
  ggplot2::ylim(c(0.9,1))+
  ggplot2::theme()

ggplot2::ggsave(plot = plot_spv, 
                filename = paste0("inst/traumacare_example/images/SPV_",outcome_name, "_seed", seed,".pdf"), 
                width = 10, height = 7)


dimnames(mean_cardinality) <- list(
  method = c(names(SL.out$calibration_scores), "GLB GRF", "GLB GLM"), 
  r_rate = random_rate)

data_cardinality_plot <- as.data.frame.table(mean_cardinality, responseName = "value")
plot_cardinality <- ggplot2::ggplot(data = data_cardinality_plot,
                                    ggplot2::aes(x = method, y=value, color=r_rate))+
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

stop()
# 3) Train the nuisances for final evaluation ────────────────────────────────────────
# Final evaluation 
## Outcome model (Q-model)  
QAW.reg.train.final= grf::probability_forest(
  X = cbind(X,A), 
  Y = Y %>% as.factor())

Q.all.pseudo.grf_final <-  do.call(cbind,lapply(0:1, function(val) {
  new_data <- cbind(SL.out$df_new_sample[,covariates_name], val)
  colnames(new_data) <- c(covariates_name, treatment_name)
  stats::predict(QAW.reg.train.final, newdata = new_data)$predictions[,2]}))

## Propensity score (G-model)
g.reg.train.final <- grf::probability_forest(X = X,
                                             Y = A %>% as.factor())
gAW.pred.pseudo.final <- stats::predict(g.reg.train.final, newdata = SL.out$df_new_sample[, covariates_name])$predictions


r <- 2 # that is (r=0.1)
spv_final <-  array(0, dim=c(n_test, n_methods, 2)) 
mean_cardinality_final<- array(0, dim=c(n_methods)) 
i <- 0
for(method in names(SL.out$perturbed_noisy_labels)){
  cal_score <- SL.out$calibration_scores[[method]] 
  i <- i+1
  # Conformal policy learning
    quant <- stats::quantile(cal_score[,r], (1-alpha))
    binary_confidence_set <-  ifelse(SL.out$new_scores<quant, 1, 0)
    idx <- which(binary_confidence_set  != 0, arr.ind = TRUE)
    confidence_set <- split(idx[, "col"]-1, 
                            factor(idx[, "row"], levels = seq_len(nrow(binary_confidence_set))))
    
    spv_res.grf <- set_policy_value_plug_in(confidence_set, 
                                            test= pseudo.test.predict,
                                            levels=levels_A, n_test = n_test,
                                            Q_all_actions = Q.all.pseudo.grf_final,
                                            gAW.pred.spv = gAW.pred.pseudo.final, 
                                            ab = ab)
    
    spv_final[ , i, 1] <- spv_res.grf[[1]] # random 
    spv_final[ , i, 2] <- spv_res.grf[[2]] # clinicians
    
    # spv_res.glm <- set_policy_value_plug_in(confidence_set, 
    #                                         test= pseudo.test.predict,
    #                                         levels=levels_A, n_test = n_test,
    #                                         Q_all_actions = Q.all.pseudo.glm,
    #                                         gAW.pred.spv = gAW.pred.pseudo.glm, 
    #                                         ab = ab)
    # 
    # spv[ , i, 1, 2, r] <- spv_res.glm[[1]] # random 
    # spv[ , i, 2, 2, r] <- spv_res.glm[[2]] # clinicians
    
    mean_cardinality_final[i]<- width(pred_set = confidence_set)
} 

# GLB GRF 
lowers <- uppers <- matrix(0,nrow=nrow(SL.out$df_new_sample), ncol=m)
for (l in as.numeric(levels_A)){
  new_data <- data.frame(SL.out$df_new_sample[,c(covariates_name)], A = l)
  colnames(new_data) <- c(covariates_name, treatment_name)
  pred <- stats::predict(SL.out$model.glb.pf, newdata = new_data, estimate.variance = TRUE, 
                         type= "response")
  lowers[,l+1] <- pmax(pred$predictions[,2] - z * sqrt(pred$variance.estimates[,2]),0)
  uppers[,l+1] <- pmin(pred$predictions[,2] + z * sqrt(pred$variance.estimates[,2]),1)
}
uppest_lrw_bound <- apply(lowers, 1, max)
C_set_binary_naive <- ifelse(uppers>=uppest_lrw_bound, 1, 0)
indices_naive <- which(C_set_binary_naive != 0, arr.ind = TRUE)
naive.confidence_set.grf <- split(indices_naive[, "col"]-1,
                                  factor(indices_naive[, "row"],
                                         levels = seq_len(nrow(C_set_binary_naive))))


spv_res.glb.grf <- set_policy_value_plug_in(naive.confidence_set.grf, 
                                            test= SL.out$df_new_sample,
                                            levels=levels_A, n_test = n_test,
                                            Q_all_actions = Q.all.pseudo.grf_final,
                                            gAW.pred.spv = gAW.pred.pseudo.final, ab = ab)

spv_final[ , n_methods-1, 1] <- spv_res.glb.grf[[1]] # random 
spv_final[ , n_methods-1, 2] <- spv_res.glb.grf[[2]]
mean_cardinality_final[n_methods-1]<- width(pred_set = naive.confidence_set.grf)
# GLB GLM
lowers <- uppers <- matrix(0,nrow=nrow(SL.out$df_new_sample), ncol=m)
for (l in as.numeric(levels_A)){
  data_l <- data.frame(SL.out$df_new_sample[,c(covariates_name)], A = l)
  colnames(data_l) <- c(covariates_name, treatment_name)
  pred <- stats::predict(SL.out$model.glb.glm, newdata = data_l, se.fit = TRUE)
  se <- sqrt(pred$se.fit)
  lowers[,l] <- (pred$fit - z * se) %>% as.numeric()
  uppers[,l] <- (pred$fit + z * se) %>% as.numeric()
}
uppest_lrw_bound <- apply(lowers, 1, max)
C_set_binary_naive <- ifelse(uppers>=uppest_lrw_bound, 1, 0)
indices_naive <- which(C_set_binary_naive != 0, arr.ind = TRUE)
naive.confidence_set.glm <- split(indices_naive[, "col"]-1,
                                  factor(indices_naive[, "row"],
                                         levels = seq_len(nrow(C_set_binary_naive))))


spv_res.glb.glm <- set_policy_value_plug_in(naive.confidence_set.glm, 
                                            test= SL.out$df_new_sample,
                                            levels=levels_A, n_test = n_test,
                                            Q_all_actions = Q.all.pseudo.grf_final,
                                            gAW.pred.spv = gAW.pred.pseudo.final, ab = ab)

spv_final[ , n_methods, 1] <- spv_res.glb.glm[[1]] # random 
spv_final[ , n_methods, 2] <- spv_res.glb.glm[[2]]
mean_cardinality_final[n_methods]<- width(pred_set = naive.confidence_set.glm)


dimnames(spv_final) <- list(
  element = 1:n_test,
  method = c(names(SL.out$calibration_scores), "GLB GRF", "GLB GLM"),
  metric = c("Uniform", "Propensity"))

data_SPV_plot <- as.data.frame.table(spv_final, responseName = "value") %>%
  mutate(
    value = as.numeric(as.character(value)),
    element = as.numeric(as.character(element))
  )

doctors <- read.csv("inst/traumacare_example/intermediate/Clinicians_Y_after_24h_before_28d.csv")$x
naive_method <- mean(Q.all.pseudo.grf_final[cbind(1:nrow(Q.all.pseudo.grf_final),
                                                  SL.out$unweighted.naive)])

plot_spv <- ggplot2::ggplot(data = data_SPV_plot, 
                            ggplot2::aes(x=method, 
                                         y=value))+
  geom_point(shape = 3) +
  ggplot2::geom_hline(yintercept = doctors, 
                      show.legend = TRUE, color="black") +
  ggplot2::geom_hline(yintercept = naive_method, 
                      show.legend = TRUE, color="grey", linetype = "dashed") +
  ggplot2::labs(title = "Set-policy values",
                y="Set-policy value") +
  ggplot2::facet_grid(~metric)+
  ggplot2::theme()

ggplot2::ggsave(plot = plot_spv, 
                filename = paste0("inst/traumacare_example/images/SPV_",outcome_name, "_seed", seed,"_final.pdf"), 
                width = 10, height = 7)


dimnames(mean_cardinality_final) <- list(
  method = c(names(SL.out$calibration_scores), "GLB GRF", "GLB GLM"))

data_cardinality_plot <- as.data.frame.table(mean_cardinality_final, responseName = "value")
plot_cardinality <- ggplot2::ggplot(data = data_cardinality_plot,
                                    ggplot2::aes(x = method, y=value))+
  ggplot2::geom_point(shape=4, size=3.5)+
  ggplot2::labs(title = "Mean cardinalities",
                y="Mean cardinality") +
  ggplot2::ylim(c(0,2)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                     vjust = 1,
                                                     hjust = 1))

ggplot2::ggsave(plot = plot_cardinality, 
                filename = paste0("inst/traumacare_example/images/Cardinalities_",outcome_name, "_seed", seed,"_final.pdf"), 
                width = 8, height = 5)
