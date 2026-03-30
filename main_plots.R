# Plots

# Density plot 
density_plot <- ggplot2::ggplot(data_toghether,
                                ggplot2::aes(x = value, 
                                             color = as.factor(type))) +
  ggplot2::geom_density(linewidth = 1) +
  ggplot2::scale_color_viridis_d(option="magma")+
  ggplot2::facet_grid(rows=vars(mechanism),
             cols=vars(model)) +
  #gganimate:: transition_states(type, transition_length = 0, 
  #                              state_length = 1) + # remove for GIF
  ggplot2::labs(y = "Density",
       color = "Type of labels",
       linetype = "Score Type") +
  ggplot2::theme(strip.placement = "outside",
        strip.background = element_blank())#+gganimate::ease_aes('linear')
ggplot2::ggsave(density_plot, 
                filename=paste0("images/", score_name, n, "/", "density_", type, ".pdf"), 
                width = 10, height = 8)

if(synthetic_scenario){
  density_plot <- density_plot + 
    ggplot2::geom_density(data = data_true_all,
                        ggplot2::aes(x = value),
                        color = "aquamarine1",
                        linetype = "dashed",
                        inherit.aes = FALSE, linewidth=1.2)
} 

# # remove for GIF
# gganimate::anim_save(
#   paste0("images/experts/", folder, "margin/, n, "/", density_", type, ".gif"),
#   animate(density_plot, fps = 5, width = 800, height = 600, renderer = gifski_renderer())
# )
ecdf_plot <- ggplot2::ggplot(data_toghether,
                             ggplot2::aes(x = value,
                                          color = as.factor(type))) +
  ggplot2::scale_color_viridis_d(option="magma")+
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::facet_grid(rows=vars(mechanism), cols=vars(model), scales = "free") +
  ggplot2::labs(y = "ECDF",
                color = "Score Type",
                linetype = "Score Type") +
  ggplot2::theme(strip.placement = "outside",
                 strip.background = element_blank())

if(synthetic_scenario){
  ecdf_plot<- ecdf_plot + 
    ggplot2::stat_ecdf(data = data_true_all, aes(x=value),
                     color = "aquamarine1",
                     linetype = "dashed", linewidth=1.2, geom = "step")
}

ggplot2::ggsave(ecdf_plot, 
                filename=paste0("images/",score_name, n, "/","ecdf_",type,".pdf"), 
                width = 10, height = 8)


if(synthetic_scenario){
  cov_vec <- mean_width <- cov_relaxed <- array(0, dim=c(length(alphas), 5,
                                                         ncol(SL.out$rate_cal_labels_exp)))
  cov_unif<- spv<- array(0, dim=c(length(alphas), n_test, 5,
                                  ncol(SL.out$rate_cal_labels_exp)))
} else {
  mean_width <- array(0, dim=c(length(alphas), 4,
                               ncol(SL.out$rate_cal_labels_exp)))
  spv<- array(0, dim=c(length(alphas), n_test, 4,
                       ncol(SL.out$rate_cal_labels_exp)))
}

# For uppest-lower bound 
treatment<- matrix(0, nrow=nrow(rbind(train1, train2)), ncol=m)
treatment[cbind(1:nrow(treatment), c(train1[, treatment_name],
                                    train2[, treatment_name]))] <- 1

model <- grf::regression_forest(
  X = cbind(rbind(train1, train2)[, covariates_name], treatment),
  Y = rbind(train1, train2)[, outcome_name])

for(i in 1:length(alphas)){
  alpha <- alphas[i]
  for (r in 1:ncol(SL.out$rate_cal_labels_exp)){
    # unweighted
    quant <- stats::quantile(SL.out$rate_scores_unweighted_cal[, r], (1-alpha))
    binary_confidence_set <-  ifelse(SL.out$new_scores<quant, 1, 0)
    idx <- which(binary_confidence_set  != 0, arr.ind = TRUE) 
    confidence_set <- split(idx[, "col"], 
                            factor(idx[, "row"], 
                                   levels = seq_len(nrow(binary_confidence_set))))
    
    mean_width[i,1,r]<- width(pred_set = confidence_set)
    spv[i,,1,r]<- bounds_set_policy_value(confidence_set, ab = ab,
                                          test= SL.out$df_new, levels=levels_A,
                                          QAW.reg.train=SL.out$QAW.reg.train,
                                          g.reg.train = SL.out$g.reg.train)
    
    # SL 
    w.quantile <- stats::quantile(SL.out$rate_scores_sl_cal[,r], (1-alpha))
    w.binary_confidence_set <- ifelse(SL.out$new_scores<w.quantile, 1, 0)
    w.idx <- which(w.binary_confidence_set  != 0, arr.ind = TRUE) 
    w.confidence_set <- split(w.idx[, "col"], 
                              factor(w.idx[, "row"], 
                                     levels = seq_len(nrow(w.binary_confidence_set))))
    
    mean_width[i,2,r]<- width(pred_set = w.confidence_set)
    spv[i,,2,r]<- bounds_set_policy_value(w.confidence_set, ab = ab,
                                          test= SL.out$df_new, levels=levels_A,
                                          QAW.reg.train=SL.out$QAW.reg.train,
                                          g.reg.train = SL.out$g.reg.train)
    
    # Exponential  
    exp.quantile <- stats::quantile(SL.out$rate_scores_exp_cal[,r], (1-alpha))
    exp.binary_confidence_set <- ifelse(SL.out$new_scores<exp.quantile, 1, 0)
    exp.idx <- which(exp.binary_confidence_set != 0, arr.ind = TRUE) 
    exp.confidence_set <- split(exp.idx[, "col"], 
                                factor(exp.idx[, "row"], 
                                       levels = seq_len(nrow(exp.binary_confidence_set))))
    
    mean_width[i,3,r]<- width(pred_set = exp.confidence_set)
    spv[i,,3,r]<- bounds_set_policy_value(exp.confidence_set, ab = ab,
                                          test= SL.out$df_new, levels=levels_A,
                                          QAW.reg.train=SL.out$QAW.reg.train,
                                          g.reg.train = SL.out$g.reg.train)
    if(synthetic_scenario){
      # Unweighted 
      cov_vec[i,1,r]<- coverage(true_set = SL.out$optimal_policy_new, 
                                pred_set = confidence_set)
      cov_relaxed[i,1,r]<- coverage_relaxed(true_set = SL.out$optimal_policy_new, 
                                            pred_set = confidence_set)
      cov_unif[i,,1,r]<- replicate(n_test,coverage_unif(SL.out$optimal_policy_new, 
                                                        pred_set=confidence_set))
      # SL 
      cov_vec[i,2,r]<- coverage(true_set = SL.out$optimal_policy_new, 
                                pred_set = w.confidence_set)
      cov_relaxed[i,2,r]<- coverage_relaxed(true_set = SL.out$optimal_policy_new, 
                                            pred_set = w.confidence_set)
      cov_unif[i,,2,r]<- replicate(n_test,coverage_unif(SL.out$optimal_policy_new,
                                                        pred_set = w.confidence_set))
      # Exponential 
      cov_vec[i,3,r]<- coverage(true_set = SL.out$optimal_policy_new, 
                                pred_set = exp.confidence_set)
      cov_relaxed[i,3,r]<- coverage_relaxed(true_set = SL.out$optimal_policy_new, 
                                            pred_set = exp.confidence_set)
      cov_unif[i,,3,r]<- replicate(n_test,coverage_unif(SL.out$optimal_policy_new,
                                                        pred_set = exp.confidence_set))
    }
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
  
  mean_width[i,4,]<- width(pred_set = naive.confidence_set)
  spv[i,,4,]<- bounds_set_policy_value(naive.confidence_set, ab = ab,
                                      test= SL.out$df_new, levels=levels_A,
                                      QAW.reg.train=SL.out$QAW.reg.train,
                                      g.reg.train = SL.out$g.reg.train)
  
  if(synthetic_scenario){
    # uppest-lower bound set  
    cov_vec[i,4,]<- coverage(true_set = SL.out$optimal_policy_new, 
                             pred_set = naive.confidence_set)
    cov_relaxed[i,4,]<- coverage_relaxed(true_set = SL.out$optimal_policy_new, 
                                         pred_set = naive.confidence_set)
    cov_unif[i,,4,]<- replicate(n_test,coverage_unif(SL.out$optimal_policy_new,
                                                     pred_set=naive.confidence_set))
   # true 
    true_quant <- stats::quantile(SL.out$true_score, (1-alpha))
    true_binary_confidence_set <-  ifelse(SL.out$new_scores<true_quant, 1, 0)
    true_idx <- which(true_binary_confidence_set  != 0, arr.ind = TRUE) 
    true_confidence_set <- split(true_idx[, "col"], 
                                 factor(true_idx[, "row"], 
                                        levels = seq_len(nrow(true_binary_confidence_set))))
    
    mean_width[i,5,]<- width(pred_set = true_confidence_set)
    cov_vec[i,5,]<- coverage(true_set = SL.out$optimal_policy_new, 
                             pred_set = true_confidence_set)
    cov_relaxed[i,5,]<- coverage_relaxed(true_set = SL.out$optimal_policy_new, 
                                         pred_set = true_confidence_set)
    cov_unif[i,,5,]<- replicate(n_test, coverage_unif(SL.out$optimal_policy_new,
                                                     pred_set = true_confidence_set))
    spv[i,,5,]<- bounds_set_policy_value(true_confidence_set, ab = ab, 
                                         test= SL.out$df_new, levels=levels_A,
                                         QAW.reg.train=SL.out$QAW.reg.train,
                                         g.reg.train = SL.out$g.reg.train) 
  }
  
  print(i)
}

if(synthetic_scenario){
  results <- list(mean_width= mean_width,
                  cov_vec=cov_vec, cov_relaxed=cov_relaxed, 
                  cov_unif=cov_unif, spv=spv)
} else{
  results <- list(mean_width= mean_width, spv=spv)
}

saveRDS(object = results, file = paste0("experts_pred/", score_name, 
                                        "plot_results_", type, "_", n, ".rds"))

#results <- readRDS(paste0("experts_pred/", score_name, 
#                          "plot_results_", type, ".rds"))

random_rate <- seq(0,1, 0.1)
spv_data <- dplyr::bind_rows(
  make_block(1, "Unweighted", results[["spv"]], alphas, random_rate), 
  make_block(2, "SL", results[["spv"]], alphas, random_rate),
  make_block(3, "Exp", results[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results[["spv"]])[1], function(a) {
    data.frame(
      value = results[["spv"]][a, ,4,1],
      mechanism = "Naive",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
)

if(synthetic_scenario){
  spv_data <- dplyr::bind_rows(spv_data, 
                        map_dfr(1:dim(results[["spv"]])[1], function(a) {
                          data.frame(
                            value = results[["spv"]][a, , 5, 1],
                            mechanism = "True",
                            level = paste0(alphas[a]),
                            type = paste0(random_rate[1]))})
                        )
}

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
unweighted_SPV <-bounds_set_policy_value(unweighted_set, test= SL.out$df_new,
                                         QAW.reg.train=SL.out$QAW.reg.train,
                                         g.reg.train = SL.out$g.reg.train, 
                                         ab = ab) 

# SL SPV
sl.probs<- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_new,
                                   weights =SL.out$w.coef, 
                                  df_pred = SL.out$df_new, 
                                  levels=as.numeric(levels_A))
est_set_sl <-  replicate(n_test,apply(sl.probs, 1,
                                      function(x)which(stats::rmultinom(1,1,prob=x)!=0)))
sl_set <- apply(est_set_sl, 1, unique)
sl_SPV <- bounds_set_policy_value(sl_set, test= SL.out$df_new,
                                  QAW.reg.train=SL.out$QAW.reg.train,
                                  g.reg.train = SL.out$g.reg.train, ab = ab)
# Exp SPV
exp_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_new,
                                    weights =SL.out$exp_weights, 
                                    df_pred = SL.out$df_new, 
                                    levels=as.numeric(levels_A))
est_set_exp <-  replicate(n_test,
                          apply(exp_probs, 1,
                                function(x)which(stats::rmultinom(1,1,prob=x)!=0)))
exp_set <- apply(est_set_exp, 1, unique)
exp_SPV <- bounds_set_policy_value(exp_set, test= SL.out$df_new,
                                   QAW.reg.train=SL.out$QAW.reg.train,
                                   g.reg.train = SL.out$g.reg.train, ab = ab)

hline_labels <- data.frame(
  x = 1,  # start of line (we'll extend it using xend)
  xend = max(spv_means$level),  # end of line
  y = c(mean(unweighted_SPV), mean(sl_SPV), mean(exp_SPV)),
  type = c("Unweighted probs", "SL probs", "Exp probs"),  # used for legend
  fill = scales::hue_pal()(3)
)
# True SPV
if(synthetic_scenario){
  true_SPV <-  bounds_set_policy_value(SL.out$optimal_policy_new,
                                       test= SL.out$df_new,
                                       QAW.reg.train=SL.out$QAW.reg.train,
                                       g.reg.train = SL.out$g.reg.train, ab = ab)
  hline_labels <- data.frame(
    x = 1,  # start of line (we'll extend it using xend)
    xend = max(spv_means$level),  # end of line
    y = c(mean(unweighted_SPV), mean(sl_SPV), mean(exp_SPV), mean(true_SPV)),
    type = c("Unweighted probs", "SL probs", "Exp probs", "True SPV"),  # used for legend
    fill = scales::hue_pal()(4)
  )
}


spv_boxplot <- ggplot2::ggplot(spv_data, aes(x = factor(level), y = value)) +
  ggplot2::geom_boxplot(aes(fill = type),position = position_dodge(width = 0.75)) +
  ggplot2::geom_segment(data = hline_labels,
                        ggplot2::aes(x = x, xend = xend, y = y, yend = y, color = type),
                        linetype = "dashed", linewidth = 1) +
  ggplot2::scale_fill_manual(
    values = c(
      stats::setNames(viridisLite::viridis(length(unique(spv_data$type)), 
                                           option = "magma"), 
                      unique(spv_data$type)),
      stats::setNames(hline_labels$fill, hline_labels$type)
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      stats::setNames(viridisLite::viridis(length(unique(spv_data$type)), 
                                           option = "magma"), 
                      unique(spv_data$type)),
      stats::setNames(hline_labels$fill, hline_labels$type)
    )
  ) +
  ggplot2::facet_grid(rows=vars(mechanism)) +
  ggplot2::labs(title = "SPV", x = "Level", y = "Set policy value") +
  ggplot2::theme_minimal()

# Plot
spv_lines <- ggplot2::ggplot(spv_means, 
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
      stats::setNames(hline_labels$fill, hline_labels$type)
    )
  ) +
  ggplot2::facet_grid(rows=vars(mechanism)) +
  ggplot2::labs(title = "SPV",
       x = "Level",
       y = "Set policy value",
       color = "Legend") +
  ggplot2::theme_minimal()

spv_plot <- gridExtra::grid.arrange(spv_boxplot, spv_lines, ncol=2)
ggplot2::ggsave(spv_plot, 
                filename=paste0("images/", score_name,  n, "/","spv_plot_",type,".pdf"), 
                width = 30, height = 15)

mean_width_data <- dplyr::bind_rows(
  make_smaller_block(1, "Unweighted", results[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  make_smaller_block(2, "SL", results[["mean_width"]], random_rate)%>%  group_by(mechanism,type)%>%mutate(levels=(row_number()-1)/(length(alphas)-1))%>% ungroup(),
  make_smaller_block(3, "Exp", results[["mean_width"]], random_rate)%>%  group_by(mechanism,type)%>%mutate(levels=(row_number()-1)/(length(alphas)-1))%>% ungroup(),
  data.frame(
    value = results[["mean_width"]][, 4, 1],
    mechanism = "True",
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
  ggplot2::facet_grid(rows=vars(mechanism)) +
  ggplot2::labs(title="Mean width", x = "Level", y = "Coverage") +
  ggplot2::theme_minimal()

if(synthetic_scenario){
  cov_data <- dplyr::bind_rows(
    make_block(1, "Unweighted", results[["cov_unif"]], alphas, random_rate), 
    make_block(2, "SL", results[["cov_unif"]], alphas, random_rate),
    make_block(3, "Exp", results[["cov_unif"]], alphas, random_rate), 
    map_dfr(1:dim(results[["cov_unif"]])[1], function(a) {
      data.frame(
        value = results[["cov_unif"]][a, , 4, 1],
        mechanism = "True",
        level = paste0(alphas[a]),
        type = paste0(random_rate[1])
      )})
  )
  
  cov_means <- cov_data %>%
    dplyr::group_by(mechanism, level, type) %>%
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  
  
  cov_plot <- ggplot2::ggplot(cov_data, 
                              ggplot2::aes(x = factor(level), y = value, 
                                           fill = type)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_line(data = cov_means,
                       ggplot2::aes(x = factor(level), y = mean_value, 
                                    group = type, colour = type)) +
    ggplot2::geom_abline(slope = -1/length(alphas), colour= "#00BFC4", 
                         intercept = 1, linetype = "dashed",
                         linewidth = 2, alpha=0.7) +
    ggplot2::scale_color_viridis_d(option = "magma") +
    ggplot2::scale_fill_viridis_d(option = "magma") +
    ggplot2::facet_grid(rows=vars(mechanism)) +
    ggplot2::labs(title = "Coverage", x = "Level", y = "Coverage") +
    ggplot2::theme_minimal()
  cov_width_plot <- gridExtra::grid.arrange(cov_plot, mean_width_plot, ncol=2)
  ggplot2::ggsave(cov_width_plot, 
                  filename=paste0("images/", score_name, n, "/", "/coverage-width_boxplots_", type, ".pdf"), 
                  width = 30, height = 15)
  
  
  cov_factor_data <- cov_means %>%
    dplyr::mutate(
      level = as.numeric(level),  # numeric for x-axis
      type = factor(type),         # factor for discrete color
      mechanism = factor(mechanism)) %>% 
    dplyr::mutate(cov_factor = mean_value - (1 - level))
  
  mean_lines <- cov_factor_data %>%
    dplyr::group_by(mechanism, type) %>%
    dplyr::summarise(mean_val = mean(cov_factor, na.rm = TRUE), .groups = "drop")
  
  cov_factor_plot <- ggplot2::ggplot(cov_factor_data, 
                                     ggplot2::aes(x = level, y = cov_factor, 
                                                  color = type, group = type)) +
    ggplot2::geom_line(linewidth = 1.2, alpha = 0.3) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(rows = vars(mechanism)) +
    # horizontal lines per mechanism and type
    ggplot2::geom_hline(data = mean_lines, 
                        ggplot2::aes(yintercept = mean_val, color = type),
                        linetype = "dashed", linewidth = 1) +
    # horizontal line at y = 0
    ggplot2::geom_hline(yintercept = 0, color = "aquamarine1", 
                        linetype = "solid", linewidth = 1) +
    ggplot2::scale_color_viridis_d(option = "magma") +
    ggplot2::labs(
      title = "Coverage factor",
      x = "Level",
      y = "Coverage Factor",
      color = "Random Rate") +
    ggplot2::theme_minimal(base_size = 14)
  
  ggplot2::ggsave(cov_factor_plot, filename=paste0("images/",score_name, n, "/","coverage_factor_plot_",type,".pdf"), width = 10, height = 8)
  
  
  all_data <- dplyr::left_join(spv_means, cov_means, 
                               by=c("level", "mechanism", "type"), 
                               suffix = c("_spv", "_cov"))
  
  
  plot_2D <- ggplot2::ggplot(all_data, 
                             ggplot2::aes(x=mean_value_spv,y=mean_value_cov,
                                  color = type)) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(rows = vars(mechanism), cols=vars(level), scales="fixed") +
    ggplot2::scale_color_viridis_d(option = "magma") +
    ggplot2::labs(
      x = "SPV",
      y = "Coverage",
      color = "Random Rate"
    ) +
    ggplot2::theme_minimal(base_size = 14)
  
  
}else{
  ggplot2::ggsave(mean_width_plot, 
                  filename=paste0("images/", score_name, n,"/", "width_boxplots_", type, ".pdf"), 
                  width = 30, height = 15)
}




