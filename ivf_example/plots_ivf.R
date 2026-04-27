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
spv_y_random <- spv_y_min <- spv_xi_random <- spv_xi_min <- array(0, dim=c(length(alphas), n_test, 2,
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
    spv_results <- ivf_set_policy_values(confidence_set, ab = ab,
                                         test= SL.out$df_new, levels=levels_A,
                                         treatment_name = treatment_name,
                                         outcome_name = outcome_name,
                                         covariates = covariates_name, 
                                         second_outcome = second_outcome_name, 
                                         mod_y=SL.out$QAW.reg.train, 
                                         mod_xi= SL.out$QxiAW.reg.train, 
                                         mod_ps = SL.out$g.reg.train)
    spv_y_random[i,,1,r]<- spv_results[["results_random_Y"]]
    spv_xi_random[i,,1,r]<- spv_results[["results_random_xi"]]
    spv_y_min[i,,1,r]<- spv_results[["results_min_Y"]]
    spv_xi_min[i,,1,r]<- spv_results[["results_min_xi"]]
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
  spv_results <- ivf_set_policy_values(naive.confidence_set, ab = ab,
                                       test= SL.out$df_new, levels=levels_A,
                                       treatment_name = treatment_name,
                                       outcome_name = outcome_name,
                                       covariates = covariates_name, 
                                       second_outcome = second_outcome_name, 
                                       mod_y=SL.out$QAW.reg.train, mod_xi = SL.out$QxiAW.reg.train,
                                       mod_ps = SL.out$g.reg.train)
  spv_y_random[i,,2,]<- spv_results[["results_random_Y"]]
  spv_xi_random[i,,2,]<- spv_results[["results_random_xi"]]
  spv_y_min[i,,2,]<- spv_results[["results_min_Y"]]
  spv_xi_min[i,,2,]<- spv_results[["results_min_xi"]]
  print(i)
}

results <- list(mean_width= mean_width, spv_y_random=spv_y_random, 
                spv_xi_random = spv_xi_random, spv_y_min = spv_y_min, 
                spv_xi_min = spv_xi_min, heatmaps_r=heatmaps_r)

saveRDS(object = results, 
        file = paste0("predictions/plot_results_", type, ".rds"))

spv_data_Y <- dplyr::bind_rows(
  make_block(1, "Unweighted", results[["spv_y_random"]], alphas, random_rate) %>% 
    rename(value_Y=value) %>% 
    mutate(choice="Random"), 
  map_dfr(1:dim(results[["spv_y_random"]])[1], function(a) {
    data.frame(
      value_Y = results[["spv_y_random"]][a, ,2,1],
      mechanism = "GLB",
      choice = "Random", 
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}), 
  make_block(1, "Unweighted", results[["spv_y_min"]], alphas, random_rate) %>%
    rename(value_Y = value) %>%
    mutate(choice="Lowest"), 
  map_dfr(1:dim(results[["spv_y_min"]])[1], function(a) {
    data.frame(
      value_Y = results[["spv_y_min"]][a, ,2,1],
      mechanism = "GLB",
      choice = "Lowest",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
)

spv_data_xi <- bind_rows( make_block(1, "Unweighted", results[["spv_xi_random"]], alphas, random_rate) %>% 
  rename(value_xi = value) %>% 
  mutate(choice="Random"), 
  map_dfr(1:dim(results[["spv_xi_random"]])[1], function(a) {
  data.frame(
    value_xi = results[["spv_xi_random"]][a, ,2,1],
    mechanism = "GLB",
    choice = "Random", 
    level = paste0(alphas[a]),
    type = paste0(random_rate[1])
  )}),
  make_block(1, "Unweighted", results[["spv_xi_min"]], alphas, random_rate) %>% 
  rename(value_xi = value) %>% 
  mutate(choice="Lowest"), 
  map_dfr(1:dim(results[["spv_xi_min"]])[1], function(a) {
  data.frame(
    value_xi = results[["spv_xi_min"]][a, ,2,1],
    mechanism = "GLB",
    choice = "Lowest",
    level = paste0(alphas[a]),
    type = paste0(random_rate[1])
  )}))

spv_data <- list(spv_data_Y, spv_data_xi) %>% 
  reduce(full_join, by = c("mechanism","level", "type", "choice")) %>%
  mutate(color_group = case_when(
    mechanism == "Unweighted" ~ paste0("type_", type),
    mechanism == "GLB" ~ "GLB"
  ))

type_vals <- sort(unique(spv_data$type))

spv_classic<- ivf_set_policy_values(SL.out$doptFactorPredict_new, ab = ab,
                                     test= SL.out$df_new, levels=levels_A,
                                     treatment_name = treatment_name,
                                     outcome_name = outcome_name,
                                     covariates = covariates_name, 
                                     second_outcome = second_outcome_name, 
                                     mod_y=SL.out$QAW.reg.train, 
                                    mod_xi = SL.out$QxiAW.reg.train,
                                     mod_ps = SL.out$g.reg.train)

hline_labels <- data.frame(
    x = 1, xend = max(spv_data$level), 
    y = spv_classic[["results_random_Y"]],
    type = c("Classic policy"),  
    fill = scales::hue_pal()(1))
  

spv_plot <- ggplot2::ggplot(spv_data, 
                            ggplot2::aes(x = factor(level),
                                         y= value_Y,
                                         color = color_group)) +
  ggplot2::geom_line(ggplot2::aes(group = color_group), 
                     position = position_dodge(width = 0.75)) +
  ggplot2::geom_point(ggplot2::aes(group = color_group), 
                      position = position_dodge(width = 0.75)) +
  ggplot2::geom_segment(data = hline_labels,
                        ggplot2::aes(x = x, xend = xend, y = y, yend = y, color="red"),
                        linetype = "dashed", linewidth = 1) +
  scale_color_manual(
    name = "Technique",
    values = c(
      stats::setNames(
        viridisLite::viridis(length(type_vals), option = "magma"),
        paste0("type_", type_vals)
      ),
      "GLB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "GLB"),
    labels = c(paste0("r = ", type_vals), "GLB")
  ) +
  ggplot2::facet_grid(rows=ggplot2::vars(choice)) +
  ggplot2::labs(x = expression("Confidence level ("* alpha *")"),
                y = "Set policy value (SPV)",
                color = "Legend") +
  ggplot2::theme_minimal()

ggplot2::ggsave(spv_plot, 
                filename=paste0("images/spv_plot_",type,".pdf"), 
                width = 30, height = 15)

level_choice <- 0.1
spv_Y_xi_plot <- ggplot2::ggplot(spv_data %>% filter(level==level_choice), 
                            ggplot2::aes(x = value_Y,
                                         y= value_xi,
                                         color = color_group, 
                                         group=color_group)) +
  ggplot2::geom_point(aes(shape=choice)) +
  ggplot2::geom_line(aes(group=color_group))+
  scale_color_manual(
    name = "Technique",
    values = c(
      stats::setNames(
        viridisLite::viridis(length(type_vals), option = "magma"),
        paste0("type_", type_vals)
      ),
      "GLB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals),  "GLB"),
    labels = c(paste0("r = ", type_vals), "GLB")
  ) +
  ggplot2::labs(x = expression("Set policy value ("* Y *")"),
                y = expression("Set policy value ("* xi *")"),
                color = "Legend") +
  ggplot2::theme_minimal()

ggplot2::ggsave(spv_Y_xi_plot, 
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
    dplyr::ungroup() ) %>% 
  mutate(color_group = case_when(
  mechanism == "Unweighted" ~ paste0("type_", type),
  mechanism == "GLB" ~ "GLB"
))

mean_width_plot <- ggplot2::ggplot(data=mean_width_data, 
                                   ggplot2::aes(x=factor(levels), y=value, 
                                                color=color_group))+
  ggplot2::geom_line(aes(group=color_group), alpha=0.5)+
  ggplot2::geom_point(data=mean_width_data%>% filter(mechanism=="Unweighted"),
                      aes(x=factor(levels), y=value, 
                          color=color_group, group=color_group))+
  ggplot2::geom_point(data = mean_width_data%>% filter(mechanism=="GLB"),
                      aes(x=factor(levels), y=value, 
                          color=color_group, group=color_group), shape=4)+
  scale_color_manual(
    name = "Technique",
    values = c(
      stats::setNames(
        viridisLite::viridis(length(type_vals), option = "magma"),
        paste0("type_", type_vals)
      ),
      "GLB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "Oracular CP", "GLB"),
    labels = c(paste0("r = ", type_vals), "Oracular CP", "GLB")
  ) +
  ggplot2::labs(x = expression("Confidence level ("* alpha *")"), 
                y = "Mean cardinality") +
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
    filename = paste0("images/", "Heatmap_", names_experts[t], "_", type, ".pdf"),
    multi_page, width = 30, height = 15)
}
