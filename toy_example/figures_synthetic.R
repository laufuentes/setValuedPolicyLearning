# non-RCT data 
results_6000_nonrct <- readRDS(paste0("predictions/plot_results_", type, "_6000", ".rds"))
SL.out_6000_nonrct <- readRDS(paste0("predictions/", type, "_6000.rds"))
results_12000_nonrct <- readRDS(paste0("predictions/plot_results_", type, "_12000", ".rds"))
SL.out_12000_nonrct <- readRDS(paste0("predictions/", type, "_12000.rds"))
results_18000_nonrct <- readRDS(paste0("predictions/plot_results_", type, "_18000", ".rds"))
SL.out_18000_nonrct <- readRDS(paste0("predictions/", type, "_18000.rds"))

# spv
spv_data_6000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_6000_nonrct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_6000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["spv"]][a, ,2,1],
      mechanism = "GLB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_6000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["spv"]][a, , 3, 1],
      mechanism = "Oracular CP",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size =6000)

spv_data_12000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_12000_nonrct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_12000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["spv"]][a, ,2,1],
      mechanism = "GLB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_12000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["spv"]][a, , 3, 1],
      mechanism = "Oracular CP",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 12000) 

spv_data_18000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_18000_nonrct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_18000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["spv"]][a, ,2,1],
      mechanism = "GLB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_18000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["spv"]][a, , 3, 1],
      mechanism = "Oracular CP",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 18000)

# mean width
mean_width_data_6000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", results_6000_nonrct[["mean_width"]], 
                     random_rate) %>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_6000_nonrct[["mean_width"]][, 2, 1],
    mechanism = "GLB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() ,
  data.frame(
    value = results_6000_nonrct[["mean_width"]][, 3, 1],
    mechanism = "Oracular CP",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=6000)

mean_width_data_12000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", 
                     results_12000_nonrct[["mean_width"]], 
                     random_rate) %>% 
    group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_12000_nonrct[["mean_width"]][, 2, 1],
    mechanism = "GLB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() , 
  data.frame(
    value = results_12000_nonrct[["mean_width"]][, 3, 1],
    mechanism = "Oracular CP",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=12000)

mean_width_data_18000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", 
                     results_18000_nonrct[["mean_width"]], 
                     random_rate) %>% 
    group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_18000_nonrct[["mean_width"]][, 2, 1],
    mechanism = "GLB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup()
  ,data.frame(
    value = results_18000_nonrct[["mean_width"]][, 3, 1],
    mechanism = "Oracular CP",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=18000)

# Coverage plot
cov_data_6000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", 
                     results_6000_nonrct[["cov_unif"]], 
                     random_rate)%>% 
    group_by(mechanism,type)%>%
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  map_dfr(1:dim(results_6000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["cov_unif"]][a, 2, 1],
      mechanism = "GLB",
      type = paste0(random_rate[1])
    )})%>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup() , 
  map_dfr(1:dim(results_6000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["cov_unif"]][a, 3, 1],
      mechanism = "Oracular CP",
      type = paste0(random_rate[1])
    )})%>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup()
) %>% mutate(size=6000)

cov_data_12000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", 
                     results_12000_nonrct[["cov_unif"]], 
                     random_rate)%>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  map_dfr(1:dim(results_12000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["cov_unif"]][a, 2, 1],
      mechanism = "GLB",
      type = paste0(random_rate[1])
    )})%>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  map_dfr(1:dim(results_12000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["cov_unif"]][a, 3, 1],
      mechanism = "Oracular CP",
      type = paste0(random_rate[1])
    )})%>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup()
) %>% mutate(size=12000)

cov_data_18000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", 
                     results_18000_nonrct[["cov_unif"]], 
                     random_rate)%>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(),  
  map_dfr(1:dim(results_18000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["cov_unif"]][a, 2, 1],
      mechanism = "GLB",
      type = paste0(random_rate[1])
    )})%>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  map_dfr(1:dim(results_18000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["cov_unif"]][a, 3, 1],
      mechanism = "Oracular CP",
      type = paste0(random_rate[1])
    )})%>% group_by(mechanism,type)%>% 
    mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup()
) %>% mutate(size=18000)


cov_data <- bind_rows(cov_data_6000_nonrct,
                      cov_data_12000_nonrct, 
                      cov_data_18000_nonrct) 
# Merge data frames 
spv_data <- bind_rows(spv_data_6000_nonrct, 
                      spv_data_12000_nonrct, 
                      spv_data_18000_nonrct)%>%
  rename(spv=value)%>%
  mutate(across(-spv, as.factor))

spv_means <- spv_data %>%
  dplyr::group_by(mechanism,level, type, size) %>%
  dplyr::summarise(mean_spv = mean(spv, na.rm = TRUE),
                   .groups = "drop")

mean_width_data <- bind_rows(mean_width_data_6000_nonrct,
                             mean_width_data_12000_nonrct,
                             mean_width_data_18000_nonrct) %>% 
  rename("mean_width"=value) %>% 
  mutate(across(-mean_width, as.factor))
  #mutate(mechanism_level = as.factor(paste(mechanism, level, sep = "_"))) %>%
  #select(-mechanism, -level)


cov_factor_data <- cov_data %>%
  dplyr::mutate(
    type = factor(type),         # factor for discrete color
    mechanism = factor(mechanism)) %>% 
  dplyr::mutate(cov_factor = cov_data$value - (1 - level)) %>% 
  mutate(across(-c(value,cov_factor ), as.factor))%>% 
  rename("cov_mean"=value)


complete_data <- list(spv_means, mean_width_data,cov_factor_data) %>%
  reduce(full_join, by = c("mechanism","level", "type", "size"))

# Random noise 
plot_data <- complete_data %>%
  mutate(color_group = case_when(
    mechanism == "Estimated labels" ~ paste0("type_", type),
    mechanism == "Oracular CP" ~ "Oracular CP",
    mechanism == "GLB" ~ "GLB"
  ))

type_vals <- sort(unique(plot_data$type))
plot_combined_level <- ggplot(plot_data,
                              aes(x = level,
                                  y = mean_spv)) +
  
  geom_line(data = subset(plot_data, mechanism == "Estimated labels"),
            aes(group = type, color = color_group),
            alpha = 0.7) +
  geom_point(aes(size = mean_width, color = color_group),
             alpha = 0.5, show.legend = c(size=FALSE)) +
  scale_color_manual(
    name = "Technique",
    values = c(
      stats::setNames(
        viridisLite::viridis(length(type_vals), option = "magma"),
        paste0("type_", type_vals)
      ),
      "Oracular CP" = "blue",
      "GLB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "Oracular CP", "GLB"),
    labels = c(paste0("r = ", type_vals), "Oracular CP", "GLB")
  ) +
  scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  facet_grid(~size) +
  
  labs(x = expression("Confidence level (" * alpha * ")"),
       y = "Set Policy Value (SPV)")
ggplot2::ggsave(plot_combined_level, filename=paste0("images/Level_SPV_", type,".pdf"), width = 15, height = 8)

plot_cov_level <- ggplot(plot_data,
                              aes(x = level,
                                  y = cov_mean)) +
  
  geom_line(data = subset(plot_data, mechanism == "Estimated labels"),
            aes(group = type, color = color_group),
            alpha = 0.7) +
  geom_point(aes(size = mean_width, color = color_group),
             alpha = 0.5, show.legend = c(size=FALSE)) +
  
  scale_color_manual(
    name = "Technique",
    values = c(
      stats::setNames(
        viridisLite::viridis(length(type_vals), option = "magma"),
        paste0("type_", type_vals)
      ),
      "Oracular CP" = "blue",
      "GLB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "Oracular CP", "GLB"),
    labels = c(paste0("r = ", type_vals), "Oracular CP", "GLB")
  ) +
  
  scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  facet_grid( ~ size) +
  
  labs(x = expression("Confidence level (" * alpha * ")"),
       y = expression("Coverage attained"))
ggplot2::ggsave(plot_cov_level, filename=paste0("images/Level_Coverage_", type,".pdf"), width = 15, height = 8)
 
plot_combined_cov_factor <- ggplot(plot_data,
                              aes(x = level,
                                  y = cov_factor, 
                                  group = color_group ,
                                  color=color_group)) +
  geom_line(alpha = 0.7) +
  geom_point(aes(size = mean_spv, color = color_group),
             alpha = 0.5, show.legend = c(size=FALSE)) +
  geom_hline(yintercept = 0, color="red")+
  scale_color_manual(
    name = "Technique",
    values = c(
      stats::setNames(
        viridisLite::viridis(length(type_vals), option = "magma"),
        paste0("type_", type_vals)
      ),
      "Oracular CP" = "blue",
      "GLB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "Oracular CP", "GLB"),
    labels = c(paste0("r = ", type_vals), "Oracular CP", "GLB")
  ) +

  scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  facet_grid(~ size) +
  labs(x = expression("Confidence level ("* alpha *")"),
       y = "Marginal coverage factor")
ggplot2::ggsave(plot_combined_cov_factor, filename=paste0("images/Level_cov_factor_", type,".pdf"), width = 15, height = 8)

plot_combined_level <- ggplot(plot_data,
                              aes(x = mean_width,
                                  y = mean_spv)) +
  
  geom_line(data = subset(plot_data, mechanism == "Estimated labels"),
            aes(group = type, color = color_group),
            alpha = 0.7) +
  geom_point(aes(size = mean_width, color = color_group),
             alpha = 0.5, show.legend = c(size=FALSE)) +
  scale_color_manual(
    name = "Technique",
    values = c(
      stats::setNames(
        viridisLite::viridis(length(type_vals), option = "magma"),
        paste0("type_", type_vals)
      ),
      "Oracular CP" = "blue",
      "GLB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "Oracular CP", "GLB"),
    labels = c(paste0("r = ", type_vals), "Oracular CP", "GLB")
  ) +
  scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  facet_grid(~size) +
  
  labs(x = "Mean cardinality",
       y = "Set Policy Value (SPV)")
ggplot2::ggsave(plot_combined_level, filename=paste0("images/Mean_width_SPV_", type,".pdf"), width = 15, height = 8)

plot_mean_level<- ggplot(plot_data,
                              aes(x = level,
                                  y = mean_width, 
                                  color = color_group)) +
  
  geom_line(aes(group = color_group),
            alpha = 0.7) +
  geom_point(aes(color = color_group),alpha = 0.5) +
  scale_color_manual(
    name = "Technique",
    values = c(
      stats::setNames(
        viridisLite::viridis(length(type_vals), option = "magma"),
        paste0("type_", type_vals)
      ),
      "Oracular CP" = "blue",
      "GLB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "Oracular CP", "GLB"),
    labels = c(paste0("r = ", type_vals), "Oracular CP", "GLB")
  ) +
  facet_grid(~size) +
  
  labs(x = expression("Confidence level ("* alpha *")"),
       y = "Mean cardinality")
ggplot2::ggsave(plot_mean_level, filename=paste0("images/Mean_width_level_", type,".pdf"), width = 15, height = 8)


optimal_treatments <- function(df) {
    df <- as.matrix(df)
    mat <- matrix(0, nrow = nrow(df), ncol = ncov)
    mat[, 1:2] <- df
    if(type=="normal"){
      p_o <- mu_P0_normal(mat) 
    }else{
      p_o <- mu_P0_simplex_complicated(mat) 
    }
    apply(data.frame(1:nrow(mat)), 1, function(i){paste0("{", 
                                                         paste(which(p_o[i,]==max(p_o[i,])), collapse = ","), "}")})
}
  
df <- tidyr::expand_grid(
    x = seq(-2, 2, length.out = 500),
    y = seq(-2, 2, length.out = 500))
df$optimal_treatments <- optimal_treatments(df)

plot_sythetic_scenario <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, 
                                         fill = as.factor(optimal_treatments))) +
    ggplot2::geom_raster() +
    ggplot2::labs(
      x = "X1",
      y = "X2",
      fill = "Optimal treatments"
    ) +
    ggplot2::theme_minimal()

ggplot2::ggsave(plot_sythetic_scenario, filename=paste0("images/Synthetic_data_", type,".pdf"), width = 10, height = 6)
