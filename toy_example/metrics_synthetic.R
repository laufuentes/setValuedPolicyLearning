type_vals <- sort(unique(SL.out$data_toghether$type))
val_names <- c(paste0("type_", type_vals), "Oracular CP", "Behavioral policy")
color_values <- c(
  viridisLite::viridis(length(type_vals), option = "magma"),
  "aquamarine1","green")
names(color_values) <- val_names

label_values <- c(
  paste0("r = ", type_vals), 
  "Oracular CP", 
  "Behavioral policy"
)
names(label_values) <- val_names


# Density plot 
density_plot <- ggplot2::ggplot(SL.out$data_toghether,
                                ggplot2::aes(x = value, 
                                             color = paste0("type_", type))) +
  ggplot2::geom_density(linewidth = 1) +
  ggplot2::geom_density(data = behavioral_data,
                        ggplot2::aes(x = value, color = "Oracular CP"),
                        linetype = "dashed",
                        inherit.aes = FALSE, linewidth=1.2)+
  ggplot2::geom_density(data = data_true_all,
                        ggplot2::aes(x = value, color = "Behavioral policy"),
                        linetype = "dashed",
                        inherit.aes = FALSE, linewidth=1.2)+
  ggplot2::scale_color_manual(
    name = "Technique",
    values = color_values,
    labels = label_values
  ) +
  ggplot2::facet_grid(cols=ggplot2::vars(model)) +
  #gganimate:: transition_states(type, transition_length = 0, 
  #                              state_length = 1) + # remove for GIF
  ggplot2::labs(y = "Density",
                color = "Type of labels",
                linetype = "Score Type")#+gganimate::ease_aes('linear')

ggplot2::ggsave(density_plot, 
                filename=paste0("images/", n, "/", "density_",type,".pdf"), 
                width = 10, height = 8)

ecdf_plot <- ggplot2::ggplot(SL.out$data_toghether, 
                             ggplot2::aes(x = value, color = paste0("type_", type))) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::stat_ecdf(data = data_true_all,
                     ggplot2::aes(x = value, color = "Oracular CP"),
                     linetype = "dashed", linewidth = 1.2, inherit.aes = FALSE) +
  ggplot2::stat_ecdf(data = behavioral_data, 
                     ggplot2::aes(x = value, color = "Behavioral policy"),
                     linetype = "dashed", linewidth = 1.2, geom = "step", inherit.aes = FALSE) +
  ggplot2::facet_grid(cols = ggplot2::vars(model), 
                      scales = "free") +
  ggplot2::scale_color_manual(
    name = "Technique",
    values = color_values,
    labels = label_values
  ) +
  ggplot2::labs(y = "ECDF", x = "Value")

ggplot2::ggsave(ecdf_plot, 
                  filename=paste0("images/", n, "/", "ecdf_",type,".pdf"), 
                  width = 10, height = 8)


# # remove for GIF
# gganimate::anim_save(
#   paste0("images/density_", type, ".gif"),
#   animate(density_plot, fps = 5, width = 800, height = 600, renderer = gifski_renderer())
# )



cov_unif <- mean_width <- cov_relaxed <- array(0, dim=c(length(alphas), 3,
                                                         ncol(SL.out$rate_cal_labels_unweighted)))

spv<- array(0, dim=c(length(alphas), 1, 3,
                                  ncol(SL.out$rate_cal_labels_unweighted)))
heatmaps_r <- array(0, dim=c(nrow(SL.out$df_new_sample), m, length(alphas),ncol(SL.out$rate_cal_labels_unweighted),3))

# For uppest-lower bound 
treatment<- matrix(0, nrow=nrow(rbind(train1, train2)), ncol=m)
treatment[cbind(1:nrow(treatment), c(train1[, treatment_name],
                                     train2[, treatment_name]))] <- 1
model <- grf::regression_forest(
  X = cbind(rbind(train1, train2)[, covariates_name], treatment),
  Y = rbind(train1, train2)[, outcome_name])

for(i in 1:length(alphas)){
  alpha <- alphas[i]
  for (r in 1:ncol(SL.out$rate_cal_labels_unweighted)){
    # unweighted
    quant <- stats::quantile(SL.out$rate_scores_unweighted_cal[, r], (1-alpha))
    binary_confidence_set <-  ifelse(SL.out$new_scores<quant, 1, 0)
    idx <- which(binary_confidence_set  != 0, arr.ind = TRUE) 
    confidence_set <- split(idx[, "col"], 
                            factor(idx[, "row"], 
                                   levels = seq_len(nrow(binary_confidence_set))))
    spv[i,,1,r]<- oracular_set_policy_value(confidence_set, test= SL.out$df_new, 
                                            levels=levels_A,
                                            treatment_name = treatment_name,
                                            outcome_name = outcome_name,
                                            covariates = covariates_name,
                                            test_potential_outcome= SL.out$potential_outcomes)
    
    cov_relaxed[i,1,r]<- coverage_relaxed(true_set = SL.out$optimal_policy_new, 
                                          pred_set = confidence_set)
    cov_unif[i,1,r]<- coverage_unif(SL.out$optimal_policy_new, pred_set=confidence_set)
    mean_width[i,1,r]<- width(pred_set = confidence_set)
    heatmaps_r[,,i,r,1] <- heatmap_treatments(confidence_set, levels_A) %>% as.matrix()
  }
  # uppest lower bound set 
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
  naive.confidence_set <- split(indices_naive[, "col"], indices_naive[, "row"])
  
  spv[i,,2,]<- oracular_set_policy_value(naive.confidence_set, 
                                         test= SL.out$df_new, levels=levels_A,
                                         treatment_name = treatment_name,
                                         outcome_name = outcome_name,
                                         covariates = covariates_name,
                                         test_potential_outcome = SL.out$potential_outcomes)
  cov_relaxed[i,2,]<- coverage_relaxed(true_set = SL.out$optimal_policy_new, 
                                       pred_set = naive.confidence_set)
  cov_unif[i,2,]<- coverage_unif(SL.out$optimal_policy_new,pred_set=naive.confidence_set)
  mean_width[i,2,]<- width(pred_set = naive.confidence_set)
  heatmaps_r[,,i,r,2] <- heatmap_treatments(naive.confidence_set, levels_A) %>% as.matrix()
  
  # true 
  true_quant <- stats::quantile(SL.out$true_score, (1-alpha))
  true_binary_confidence_set <-  ifelse(SL.out$new_scores<true_quant, 1, 0)
  true_idx <- which(true_binary_confidence_set  != 0, arr.ind = TRUE) 
  true_confidence_set <- split(true_idx[, "col"], 
                                 factor(true_idx[, "row"], 
                                        levels = seq_len(nrow(true_binary_confidence_set))))
    
    mean_width[i,3,]<- width(pred_set = true_confidence_set)
    cov_relaxed[i,3,]<- coverage_relaxed(true_set = SL.out$optimal_policy_new, 
                                         pred_set = true_confidence_set)
    cov_unif[i,3,]<- coverage_unif(SL.out$optimal_policy_new,pred_set = true_confidence_set)
    spv[i,,3,]<- oracular_set_policy_value(true_confidence_set, 
                                           test= SL.out$df_new, levels=levels_A,
                                           treatment_name = treatment_name,
                                           outcome_name = outcome_name,
                                           covariates = covariates_name,
                                           test_potential_outcome=SL.out$potential_outcomes)
    
    heatmaps_r[,,i,r,3] <- heatmap_treatments(true_confidence_set, levels_A) %>% as.matrix()
  print(i)
}

results <- list(mean_width= mean_width, cov_relaxed=cov_relaxed, 
                  cov_unif=cov_unif, spv=spv,heatmaps_r=heatmaps_r)
  
saveRDS(object = results, file = paste0("predictions/plot_results_", type, "_", n, ".rds"))