type_vals <- sort(unique(SL.out$data_toghether$type))
val_names <- c(paste0("type_", type_vals), "Behavioral policy")
color_values <- c(
  viridisLite::viridis(length(type_vals), option = "magma"),"green")
names(color_values) <- val_names
label_values <- c(paste0("r = ", type_vals), "Behavioral policy")
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
                filename=paste0("images/density_",type,".pdf"), 
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
                filename=paste0("images/ecdf_",type,".pdf"), 
                width = 10, height = 8)

# # remove for GIF
# gganimate::anim_save(
#   paste0("images/experts/", folder, "margin/, density_", type, ".gif"),
#   animate(density_plot, fps = 5, width = 800, height = 600, renderer = gifski_renderer())
# )

mean_width <- array(0, dim=c(length(alphas), 2,
                             ncol(SL.out$rate_cal_labels_unweighted)))
spv<- array(0, dim=c(length(alphas), n_test, 2,
                     ncol(SL.out$rate_cal_labels_unweighted)))

heatmaps_r <- array(0, dim=c(nrow(SL.out$df_new_sample), m, 
                             length(alphas),
                             ncol(SL.out$rate_cal_labels_unweighted),2))


# For greatest lower bound (GLB)
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
    
    mean_width[i,1,r]<- width(pred_set = confidence_set)
    heatmaps_r[,,i,r,1] <- heatmap_treatments(confidence_set, levels_A) %>% as.matrix()
    spv[i,,1,r]<- bounds_set_policy_value(confidence_set, ab = ab,
                                          test= SL.out$df_new, levels=levels_A,
                                          treatment_name = treatment_name,
                                          outcome_name = outcome_name,
                                          covariates = covariates_name,
                                          mod_y=SL.out$QAW.reg.train,
                                          mod_ps = SL.out$g.reg.train)
  }
  # greatest lower bound (GLB) 
  lowers <- uppers <- lowers_test <- uppers_test <- matrix(0, 
                                                           nrow=nrow(SL.out$df_new), 
                                                           ncol=m)
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
  
  mean_width[i,2,]<- width(pred_set = naive.confidence_set)
  heatmaps_r[,,i,r,2] <- heatmap_treatments(naive.confidence_set, levels_A) %>% as.matrix()
  spv[i,,2,]<- bounds_set_policy_value(naive.confidence_set, ab = ab,
                                       test= SL.out$df_new, levels=levels_A,
                                       treatment_name = treatment_name,
                                       outcome_name = outcome_name,
                                       covariates = covariates_name,
                                       mod_y=SL.out$QAW.reg.train,
                                       mod_ps = SL.out$g.reg.train)
  
  print(i)
}

results <- list(mean_width= mean_width, spv=spv, 
                heatmaps_r=heatmaps_r)

saveRDS(object = results, 
        file = paste0("predictions/plot_results_", type, ".rds"))

spv_data <- dplyr::bind_rows(
  make_block(1, "Unweighted", results[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results[["spv"]])[1], function(a) {
    data.frame(
      value = results[["spv"]][a, ,2,1],
      mechanism = "GLB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
)

spv_means <- spv_data %>%
  dplyr::group_by(mechanism, level, type) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
                   .groups = "drop")

# unweighted SPV
unweighted_probs <- weighted_probs_experts(df_pred=SL.out$df_new, 
                                           fitted_experts=SL.out$doptFactorPredict_new, 
                                           weights=rep(1/numalgs, numalgs), 
                                           levels=as.numeric(levels_A))
est_set_unweighted <-  replicate(n_test,
                                 apply(unweighted_probs, 1, 
                                       function(x) which(stats::rmultinom(1,1,prob=x)!=0)))
unweighted_set <- apply(est_set_unweighted, 1, unique)

unweighted_SPV <-bounds_set_policy_value(unweighted_set, ab = ab,
                                         test= SL.out$df_new, levels=levels_A,
                                         treatment_name = treatment_name,
                                         outcome_name = outcome_name,
                                         covariates = covariates_name,
                                         mod_y=SL.out$QAW.reg.train,
                                         mod_ps = SL.out$g.reg.train)

hline_labels <- data.frame(
  x = 1,  xend = max(spv_means$level),
  y = c(mean(unweighted_SPV)),
  type = c("Unweighted probs"), 
  fill = scales::hue_pal()(1)
)

# Plot
if(!synthetic_scenario){
  pv_macf <- policy_value_tmle(model_macf, SL.out$df_new_sample, 
                               mod_y=SL.out$QAW.reg.train, mod_ps, 
                               SL.out$g.reg.tr, 
                               ab=ab, covariates=covariates_name, 
                               treatment_name = treatment_name, 
                               outcome_name=outcome_name,
                               family="gaussian", levels=levels_A)
  hline_labels <- data.frame(
    x = 1, xend = max(spv_means$level), 
    y = c(mean(unweighted_SPV), pv_macf ),
    type = c("Unweighted probs", "Classic policy learning"),  
    fill = scales::hue_pal()(2))
}

spv_plot <- ggplot2::ggplot(spv_means, 
                             ggplot2::aes(x = factor(level), 
                                          y = mean_value, color = type)) +
  ggplot2::geom_line(ggplot2::aes(group = type, color = type), 
                     position = position_dodge(width = 0.75)) +
  ggplot2::geom_point(ggplot2::aes(group = type, color = type), 
                      position = position_dodge(width = 0.75)) +
  ggplot2::geom_segment(data = hline_labels,
                        ggplot2::aes(x = x, xend = xend, y = y, yend = y, 
                                     color = type),
                        linetype = "dashed", linewidth = 1) +
  ggplot2::scale_color_manual(
    values = c(
      stats::setNames(viridisLite::viridis(length(unique(spv_means$type)), 
                                           option = "magma"), 
                      unique(spv_means$type)),
      stats::setNames(hline_labels$fill, hline_labels$type))) +
  ggplot2::facet_grid(rows=ggplot2::vars(mechanism)) +
  ggplot2::labs(title = "SPV",
                x = "Level",
                y = "Set policy value",
                color = "Legend") +
  ggplot2::theme_minimal()

ggplot2::ggsave(spv_plot, 
                filename=paste0("images/spv_plot_",type,".pdf"), 
                width = 30, height = 15)

mean_width_data <- dplyr::bind_rows(
  make_smaller_block(1, "Unweighted", results[["mean_width"]], random_rate) %>% 
    group_by(mechanism,type)%>% 
    mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
   data.frame(
    value = results[["mean_width"]][, 2, 1],
    mechanism = "GLB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
)

mean_width_plot <- ggplot2::ggplot(data=mean_width_data, 
                                   ggplot2::aes(x=factor(levels), y=value, 
                                                color=type))+
  ggplot2::geom_line(aes(group=type))+
  ggplot2::geom_point(aes(group=type))+
  ggplot2::scale_color_viridis_d(option="magma")+
  ggplot2::facet_grid(rows=ggplot2::vars(mechanism)) +
  ggplot2::labs(title="Mean width", x = "Level", y = "Coverage") +
  ggplot2::theme_minimal()

ggplot2::ggsave(mean_width_plot, 
                  filename=paste0("images/width_boxplots_", type, ".pdf"), 
                  width = 30, height = 15)

names_experts <- c("Unweighted", "GLB")
for (t in 1:dim(heatmaps_r)[5]){
  plots_completed <- list()
  for (i in 1:dim(heatmaps_r)[3]){
    plots <- list()
    for (r in 1:dim(heatmaps_r)[4]){
      file <- as.data.frame(heatmaps_r[,,i,r,t]) %>%
        mutate(row_id = row_number()) %>%
        pivot_longer(cols = -row_id, names_to = "column_m", values_to = "value")
      plots[[r]] <- ggplot(file, aes(x = column_m, y = row_id, fill = factor(value))) +
        geom_tile() +
        theme_minimal() +
        labs(title = paste0("r: ", random_rate[r]),
             x = "Treatment (m)",
             y = "Observations",
             fill = "Present")
    }
    plots_completed[[i]] <- gridExtra::arrangeGrob(grobs = plots, nrow = 1, ncol = dim(heatmaps_r)[4], top = paste0("Alpha = ", alphas[i]))
  }
  multi_page <- marrangeGrob(grobs = plots_completed, nrow = 1, ncol = 1)
  ggplot2::ggsave(
      filename = paste0("images/", score_name, "Heatmap_", names_experts[t], "_", type, ".pdf"),
      multi_page, width = 30, height = 15)
}