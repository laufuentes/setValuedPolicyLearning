results_6000 <- readRDS(paste0("experts_pred/", score_name, RCT_file, "plot_results_", type, "_6000", ".rds"))
results_12000 <- readRDS(paste0("experts_pred/", score_name, RCT_file, "plot_results_", type, "_12000", ".rds"))
results_18000 <- readRDS(paste0("experts_pred/", score_name, RCT_file, "plot_results_", type, "_18000", ".rds"))

SL.out_6000 <- readRDS(paste0("experts_pred/", score_name, RCT_file, "_", type, "_6000.rds"))
SL.out_12000 <- readRDS(paste0("experts_pred/", score_name, RCT_file, "_", type, "_12000.rds"))
SL.out_18000 <- readRDS(paste0("experts_pred/", score_name, RCT_file, "_", type, "_18000.rds"))
numalgs <- SL.out_6000$libraryNames %>% length()

# Set policy value plot
spv_data_6000 <- dplyr::bind_rows(
  make_block(1, "Unweighted", results_6000[["spv"]], alphas, random_rate), 
  make_block(2, "SL", results_6000[["spv"]], alphas, random_rate),
  make_block(3, "Exp", results_6000[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_6000[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_6000[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000[["spv"]][a, , 5, 1],
      mechanism = "True",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size =6000)

spv_data_12000 <- dplyr::bind_rows(
  make_block(1, "Unweighted", results_12000[["spv"]], alphas, random_rate), 
  make_block(2, "SL", results_12000[["spv"]], alphas, random_rate),
  make_block(3, "Exp", results_12000[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_12000[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_12000[["spv"]])[1], function(a) {
      data.frame(
        value = results_12000[["spv"]][a, , 5, 1],
        mechanism = "True",
        level = paste0(alphas[a]),
        type = paste0(random_rate[1]))}))%>% mutate(size = 12000) 

spv_data_18000 <- dplyr::bind_rows(
  make_block(1, "Unweighted", results_18000[["spv"]], alphas, random_rate), 
  make_block(2, "SL", results_18000[["spv"]], alphas, random_rate),
  make_block(3, "Exp", results_18000[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_18000[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_18000[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000[["spv"]][a, , 5, 1],
      mechanism = "True",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 18000)

spv_data <- bind_rows(spv_data_6000, spv_data_12000, spv_data_18000)

spv_means <- spv_data %>%
  dplyr::group_by(mechanism, level, type, size) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
                   .groups = "drop")

# unweighted SPV
hline_results_list <- list()
i<- 1
sizes <- c(6000, 12000, 18000)
for (SL.out in list(SL.out_6000,SL.out_12000,SL.out_18000)){
  Y <- SL.out$df_obs[,outcome_name] 
  Y_new <- SL.out$df_new_sample[,outcome_name] 
  ab <- c(min(c(Y,Y_new)),max(c(Y,Y_new))) 
  
  unweighted_probs <- weighted_probs_experts(df_pred=SL.out$df_new, 
                                             fitted_experts=SL.out$doptFactorPredict_new, 
                                             weights=rep(1/numalgs, numalgs), 
                                             levels=as.numeric(levels_A))
  est_set_unweighted <-  replicate(n_test,
                                   apply(unweighted_probs, 1, 
                                         function(x) which(stats::rmultinom(1,1,prob=x)!=0)))
  unweighted_set <- apply(est_set_unweighted, 1, unique)
  unweighted_SPV <-oracular_set_policy_value(unweighted_set, 
                                             test= SL.out$df_new, levels=levels_A,
                                             treatment_name = treatment_name,
                                             outcome_name = outcome_name,
                                             covariates = covariates_name,
                                             test_potential_outcome = SL.out$potential_outcomes)
  
  # SL SPV
  sl.probs<- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_new,
                                    weights =SL.out$w.coef, 
                                    df_pred = SL.out$df_new, 
                                    levels=as.numeric(levels_A))
  est_set_sl <-  replicate(n_test,apply(sl.probs, 1,
                                        function(x)which(stats::rmultinom(1,1,prob=x)!=0)))
  sl_set <- apply(est_set_sl, 1, unique)
  sl_SPV <- oracular_set_policy_value(sl_set, 
                                      test= SL.out$df_new, levels=levels_A,
                                      treatment_name = treatment_name,
                                      outcome_name = outcome_name,
                                      covariates = covariates_name,
                                      test_potential_outcome = SL.out$potential_outcomes)
  
  # Exp SPV
  exp_probs <- weighted_probs_experts(fitted_experts = SL.out$doptFactorPredict_new,
                                      weights =SL.out$exp_weights, 
                                      df_pred = SL.out$df_new, 
                                      levels=as.numeric(levels_A))
  est_set_exp <-  replicate(n_test,
                            apply(exp_probs, 1,
                                  function(x)which(stats::rmultinom(1,1,prob=x)!=0)))
  exp_set <- apply(est_set_exp, 1, unique)
  exp_SPV <- oracular_set_policy_value(exp_set,
                                       test= SL.out$df_new, levels=levels_A,
                                       treatment_name = treatment_name,
                                       outcome_name = outcome_name,
                                       covariates = covariates_name,
                                       test_potential_outcome = SL.out$potential_outcomes)
  
  true_SPV <-  oracular_set_policy_value(SL.out$optimal_policy_new, 
                                       test= SL.out$df_new, levels=levels_A,
                                       treatment_name = treatment_name,
                                       outcome_name = outcome_name,
                                       covariates = covariates_name,
                                       test_potential_outcome = SL.out$potential_outcomes)
  
  temp_data <- data.frame(
    x = 1,  # start of line (we'll extend it using xend)
    xend = max(spv_means$level),  # end of line
    y = c(mean(unweighted_SPV), mean(sl_SPV), mean(exp_SPV), mean(true_SPV)),
    type = c("Unweighted", "SL", "Exp", "True"),  # used for legend
    size=sizes[i],
    fill = scales::hue_pal()(4))
  
  hline_results_list[[i]] <- temp_data
  i<- i+1
}

hline_labels<- do.call(rbind, hline_results_list)

# Mean width plot
mean_width_data_6000 <- dplyr::bind_rows(
  make_smaller_block(1, "Unweighted", results_6000[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  make_smaller_block(2, "SL", results_6000[["mean_width"]], random_rate)%>%  group_by(mechanism,type)%>%mutate(levels=(row_number()-1)/(length(alphas)-1))%>% ungroup(),
  make_smaller_block(3, "Exp", results_6000[["mean_width"]], random_rate)%>%  group_by(mechanism,type)%>%mutate(levels=(row_number()-1)/(length(alphas)-1))%>% ungroup(),
  data.frame(
    value = results_6000[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() ,
  data.frame(
    value = results_6000[["mean_width"]][, 5, 1],
    mechanism = "True",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=6000)

mean_width_data_12000 <- dplyr::bind_rows(
  make_smaller_block(1, "Unweighted", results_12000[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  make_smaller_block(2, "SL", results_12000[["mean_width"]], random_rate)%>%  group_by(mechanism,type)%>%mutate(levels=(row_number()-1)/(length(alphas)-1))%>% ungroup(),
  make_smaller_block(3, "Exp", results_12000[["mean_width"]], random_rate)%>%  group_by(mechanism,type)%>%mutate(levels=(row_number()-1)/(length(alphas)-1))%>% ungroup(),
  data.frame(
    value = results_12000[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() , 
  data.frame(
    value = results_12000[["mean_width"]][, 5, 1],
    mechanism = "True",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=12000)

mean_width_data_18000 <- dplyr::bind_rows(
  make_smaller_block(1, "Unweighted", results_18000[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  make_smaller_block(2, "SL", results_18000[["mean_width"]], random_rate)%>%  group_by(mechanism,type)%>%mutate(levels=(row_number()-1)/(length(alphas)-1))%>% ungroup(),
  make_smaller_block(3, "Exp", results_18000[["mean_width"]], random_rate)%>%  group_by(mechanism,type)%>%mutate(levels=(row_number()-1)/(length(alphas)-1))%>% ungroup(),
  data.frame(
    value = results_18000[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup()
  ,data.frame(
    value = results_18000[["mean_width"]][, 5, 1],
    mechanism = "True",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=18000)

mean_width_data <- bind_rows(mean_width_data_6000, mean_width_data_12000, mean_width_data_18000)

# Coverage plot
cov_data_6000 <- dplyr::bind_rows(
  make_block(1, "Unweighted", results_6000[["cov_unif"]], alphas, random_rate), 
  make_block(2, "SL", results_6000[["cov_unif"]], alphas, random_rate),
  make_block(3, "Exp", results_6000[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_6000[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}), 
  map_dfr(1:dim(results_6000[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000[["cov_unif"]][a, , 5, 1],
      mechanism = "True",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=6000)


cov_data_12000 <- dplyr::bind_rows(
  make_block(1, "Unweighted", results_12000[["cov_unif"]], alphas, random_rate), 
  make_block(2, "SL", results_12000[["cov_unif"]], alphas, random_rate),
  make_block(3, "Exp", results_12000[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_12000[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}),
  map_dfr(1:dim(results_12000[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000[["cov_unif"]][a, , 5, 1],
      mechanism = "True",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=12000)

cov_data_18000 <- dplyr::bind_rows(
  make_block(1, "Unweighted", results_18000[["cov_unif"]], alphas, random_rate), 
  make_block(2, "SL", results_18000[["cov_unif"]], alphas, random_rate),
  make_block(3, "Exp", results_18000[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_18000[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}),
  map_dfr(1:dim(results_18000[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000[["cov_unif"]][a, , 5, 1],
      mechanism = "True",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=18000)

cov_data <- bind_rows(cov_data_6000, cov_data_12000, cov_data_18000)

cov_means <- cov_data %>%
  dplyr::group_by(mechanism, level, type, size) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Coverage factor plot 
cov_factor_data <- cov_means %>%
  dplyr::mutate(
    level = as.numeric(level),  # numeric for x-axis
    type = factor(type),         # factor for discrete color
    mechanism = factor(mechanism)) %>% 
  dplyr::mutate(cov_factor = mean_value - (1 - level))

mean_lines <- cov_factor_data %>%
  dplyr::group_by(mechanism, type, size) %>%
  dplyr::summarise(mean_val = mean(cov_factor, na.rm = TRUE), .groups = "drop")


for(mech in c("Unweighted", "SL", "Exp", "True")){
  # SPV plot
  spv_mean_mech <- spv_means %>% filter(mechanism==mech)
  spv_lines <- ggplot2::ggplot(spv_mean_mech, 
                               ggplot2::aes(x = factor(level), 
                                            y = mean_value, color = type)) +
    ggplot2::geom_line(ggplot2::aes(group = type, color = type), 
                       position = position_dodge(width = 0.75)) +
    ggplot2::geom_point(ggplot2::aes(group = type, color = type), 
                        position = position_dodge(width = 0.75)) +
    ggplot2::geom_segment(data = hline_labels %>% filter(type %in% c(mech, "True")),
                          ggplot2::aes(x = x, xend = xend, y = y, yend = y, 
                                       color = type),
                          linetype = "dashed", linewidth = 1) +
    ggplot2::scale_color_manual(
      values = c(
        stats::setNames(
          viridisLite::viridis(length(unique(spv_means$type)), option = "magma"), 
          unique(spv_means$type)
        ),
        stats::setNames(hline_labels$fill, hline_labels$type),
        "ULB" = "blue"  
      )
    )+
    ggplot2::geom_point(
      data = spv_means %>% filter(mechanism == "ULB"), 
      ggplot2::aes(
        x = factor(level),
        y = mean_value,
        color = "ULB",  
        group = "ULB"),
      shape = 4,
      position = position_dodge(width = 0.75))+
    ggplot2::facet_grid(cols=ggplot2::vars(size)) +
    ggplot2::labs(title = "SPV",
                  x = expression("Confidence level (" * alpha * ")"),
                  y = "Set policy value",
                  color = "Legend") +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(spv_lines, filename=paste0("images/",score_name, RCT_file, "spv_lines_size_",mech, "_",type,".pdf"), width = 15, height = 8)
  
  # Mean width plot
  mean_width_data_mech <- mean_width_data %>% filter(mechanism==mech)
  mean_width_plot <- ggplot2::ggplot(data=mean_width_data_mech, 
                                     ggplot2::aes(x=factor(levels), y=value, 
                                                  color=type))+
    ggplot2::geom_line(aes(group=type))+
    ggplot2::geom_point(aes(group=type))+
    ggplot2::geom_point(
      data = mean_width_data %>% filter(mechanism == "ULB"), 
      ggplot2::aes(
        x = factor(levels), y = value, color = "ULB", group = "ULB"), 
      shape = 4) +
    ggplot2::scale_color_manual(
      values = c(
        stats::setNames(
          viridisLite::viridis(length(unique(mean_width_data_mech$type)), 
                               option = "magma"), 
          unique(mean_width_data_mech$type)),"ULB" = "blue")) +
    ggplot2::facet_grid(cols=ggplot2::vars(size)) +
    ggplot2::labs(title="Mean width", 
                  x = expression("Confidence level (" * alpha * ")"), 
                  y = "Coverage") +
    ggplot2::theme_minimal()
  
  ggplot2::ggsave(mean_width_plot, filename=paste0("images/",score_name, RCT_file, "mean_width_size_",mech,"_",type,".pdf"), width = 15, height = 8)
  
  # Coverage plot
  cov_data_mech <- cov_means %>% filter(mechanism==mech) # cov_data %>% filter(mechanism==mech)
  cov_plot <- ggplot2::ggplot(cov_data_mech, 
                              ggplot2::aes(x = factor(level), y = value, 
                                           fill = type)) +
    #ggplot2::geom_boxplot() +
    ggplot2::geom_line(ggplot2::aes(x = factor(level), y = mean_value, 
                                    group = type, colour = type)) +
    ggplot2::geom_abline(slope = -1/length(alphas), colour= "red", 
                         intercept = 1, linetype = 2,
                         linewidth = 2, alpha=0.7) +
    ggplot2::geom_point(data = cov_means %>% filter(mechanism=="ULB"),
                       ggplot2::aes(x = factor(level), y = mean_value, 
                                    color= "ULB", fill = "ULB"), shape=4)+
    ggplot2::scale_color_manual(
      values = c(
        stats::setNames(
          viridisLite::viridis(length(unique(cov_data_mech$type)), 
                               option = "magma"), 
          unique(cov_data_mech$type)),"ULB" = "blue")) +
    ggplot2::facet_grid(cols=ggplot2::vars(size)) +
    ggplot2::labs(title = "Coverage", 
                  x = expression("Confidence level (" * alpha * ")"), 
                  y = "Coverage") +
    ggplot2::theme_minimal()

  
  ggplot2::ggsave(cov_plot, filename=paste0("images/",score_name, RCT_file, "cov_size_", mech, "_",type,".pdf"), width = 15, height = 8)
  
  # Marginal coverage factor plot 
  cov_factor_data_mech <- cov_factor_data %>% filter(mechanism==mech)
  cov_factor_plot <- ggplot2::ggplot(cov_factor_data_mech, 
                                     ggplot2::aes(x = level, y = cov_factor, 
                                                  color = type, group = type)) +
    ggplot2::geom_line(linewidth = 1.2, alpha = 0.3) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_point(data=cov_factor_data %>% filter(mechanism=="ULB"),
                        ggplot2::aes(x = level, y = cov_factor, 
                                     color = "ULB", group ="ULB"), 
                        shape=4, size = 2) +
    ggplot2::facet_grid(cols=ggplot2::vars(size)) +
    # horizontal lines per mechanism and type
    ggplot2::geom_hline(data = mean_lines %>% filter(mechanism==mech), 
                        ggplot2::aes(yintercept = mean_val, color = type),
                        linetype = "dashed", linewidth = 1) +
    # horizontal line at y = 0
    ggplot2::geom_hline(yintercept = 0, color = "aquamarine1", 
                        linetype = "solid", linewidth = 1) +
    ggplot2::scale_color_manual(
      values = c(
        stats::setNames(
          viridisLite::viridis(length(unique(cov_factor_data_mech$type)), 
                               option = "magma"), 
          unique(cov_factor_data_mech$type)),"ULB" = "blue")) +
    ggplot2::labs(
      title = "Coverage factor",
      x = expression("Confidence level (" * alpha * ")"),
      y = "Coverage Factor",
      color = "Random Rate") +
    ggplot2::theme_minimal(base_size = 14)
  
  ggplot2::ggsave(cov_factor_plot, filename=paste0("images/",score_name, RCT_file, "cov_factor_size_", mech, "_", type,".pdf"), width = 15, height = 8)
}

