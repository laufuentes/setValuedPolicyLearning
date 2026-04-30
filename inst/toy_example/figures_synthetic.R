# ── Load prediction metrics  ────────────────────────────────────────────
# n = 6000
results_6000_nonrct <- readRDS(paste0("inst/predictions/plot_results_", type, "_6000", ".rds"))
SL.out_6000_nonrct <- readRDS(paste0("inst/predictions/", type, "_6000.rds"))
# n = 12000
results_12000_nonrct <- readRDS(paste0("inst/predictions/plot_results_", type, "_12000", ".rds"))
SL.out_12000_nonrct <- readRDS(paste0("inst/predictions/", type, "_12000.rds"))
# n = 18000
results_18000_nonrct <- readRDS(paste0("inst/predictions/plot_results_", type, "_18000", ".rds"))
SL.out_18000_nonrct <- readRDS(paste0("inst/predictions/", type, "_18000.rds"))

# ── Set-policy value data  ────────────────────────────────────────────
# n = 6000
spv_data_6000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_6000_nonrct[["spv"]], alphas, random_rate),
  purrr::map_dfr(1:dim(results_6000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["spv"]][a, ,2,1],
      mechanism = "GLB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}),
  purrr::map_dfr(1:dim(results_6000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["spv"]][a, , 3, 1],
      mechanism = "Oracular CP",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size =6000)
# n = 12000
spv_data_12000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_12000_nonrct[["spv"]], alphas, random_rate),
  purrr::map_dfr(1:dim(results_12000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["spv"]][a, ,2,1],
      mechanism = "GLB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}),
  purrr::map_dfr(1:dim(results_12000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["spv"]][a, , 3, 1],
      mechanism = "Oracular CP",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 12000)
# n = 18000
spv_data_18000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_18000_nonrct[["spv"]], alphas, random_rate),
  purrr::map_dfr(1:dim(results_18000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["spv"]][a, ,2,1],
      mechanism = "GLB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}),
  purrr::map_dfr(1:dim(results_18000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["spv"]][a, , 3, 1],
      mechanism = "Oracular CP",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 18000)
# merge data frames
spv_data <- dplyr::bind_rows(spv_data_6000_nonrct,
                      spv_data_12000_nonrct,
                      spv_data_18000_nonrct)%>%
  dplyr::rename(spv=value)%>%
  dplyr::mutate(across(-spv, as.factor))

spv_means <- spv_data %>%
  dplyr::group_by(mechanism, level, type, size) %>%
  dplyr::summarise(mean_spv = mean(spv, na.rm = TRUE),
                   .groups = "drop")

# ── Mean cardinality data  ────────────────────────────────────────────
# n = 6000
mean_cardinality_data_6000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", results_6000_nonrct[["mean_cardinality"]],
                     random_rate) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(),
  data.frame(
    value = results_6000_nonrct[["mean_cardinality"]][, 2, 1],
    mechanism = "GLB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1))%>%
    dplyr::ungroup() ,
  data.frame(
    value = results_6000_nonrct[["mean_cardinality"]][, 3, 1],
    mechanism = "Oracular CP",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1))%>%
    dplyr::ungroup()
) %>% dplyr::mutate(size=6000)
# n = 12000
mean_cardinality_data_12000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels",
                     results_12000_nonrct[["mean_cardinality"]],
                     random_rate) %>%
    dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(),
  data.frame(
    value = results_12000_nonrct[["mean_cardinality"]][, 2, 1],
    mechanism = "GLB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1))%>%
    dplyr::ungroup() ,
  data.frame(
    value = results_12000_nonrct[["mean_cardinality"]][, 3, 1],
    mechanism = "Oracular CP",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1))%>%
    dplyr::ungroup()
) %>% mutate(size=12000)
# n = 18000
mean_cardinality_data_18000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels",
                     results_18000_nonrct[["mean_cardinality"]],
                     random_rate) %>%
    dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1)) %>% 
    dplyr::ungroup(),
  data.frame(
    value = results_18000_nonrct[["mean_cardinality"]][, 2, 1],
    mechanism = "GLB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1))%>%
    dplyr::ungroup(),
  data.frame(
    value = results_18000_nonrct[["mean_cardinality"]][, 3, 1],
    mechanism = "Oracular CP",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1))%>%
    dplyr::ungroup()
) %>% mutate(size=18000)
# merge data frames
mean_cardinality_data <- bind_rows(mean_cardinality_data_6000_nonrct,
                             mean_cardinality_data_12000_nonrct,
                             mean_cardinality_data_18000_nonrct) %>%
  dplyr::rename("mean_cardinality"=value) %>%
  dplyr::mutate(across(-mean_cardinality, as.factor))

# ── Coverage data  ────────────────────────────────────────────
# n = 6000
cov_data_6000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels",
                     results_6000_nonrct[["cov_unif"]],
                     random_rate)%>%
    dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(),
  map_dfr(1:dim(results_6000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["cov_unif"]][a, 2, 1],
      mechanism = "GLB",
      type = paste0(random_rate[1])
    )})%>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1)) %>% 
    dplyr::ungroup() ,
  purrr::map_dfr(1:dim(results_6000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["cov_unif"]][a, 3, 1],
      mechanism = "Oracular CP",
      type = paste0(random_rate[1])
    )})%>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1)) %>% 
    dplyr::ungroup()) %>% dplyr::mutate(size=6000)
# n = 12000
cov_data_12000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels",
                     results_12000_nonrct[["cov_unif"]],
                     random_rate)%>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1)) %>% 
    dplyr::ungroup(),
  purrr::map_dfr(1:dim(results_12000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["cov_unif"]][a, 2, 1],
      mechanism = "GLB",
      type = paste0(random_rate[1])
    )})%>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1)) %>% 
    dplyr::ungroup(),
  purrr::map_dfr(1:dim(results_12000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["cov_unif"]][a, 3, 1],
      mechanism = "Oracular CP",
      type = paste0(random_rate[1])
    )})%>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1)) %>% 
    dplyr::ungroup()) %>% 
  dplyr::mutate(size=12000)
# n = 18000
cov_data_18000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels",
                     results_18000_nonrct[["cov_unif"]],
                     random_rate)%>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1)) %>% 
    dplyr::ungroup(),
  purrr::map_dfr(1:dim(results_18000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["cov_unif"]][a, 2, 1],
      mechanism = "GLB",
      type = paste0(random_rate[1])
    )})%>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(row_number()-1)/(length(alphas)-1)) %>% ungroup(),
  purrr::map_dfr(1:dim(results_18000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["cov_unif"]][a, 3, 1],
      mechanism = "Oracular CP",
      type = paste0(random_rate[1])
    )})%>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(level=(dplyr::row_number()-1)/(length(alphas)-1)) %>% 
    dplyr::ungroup()) %>% 
  dplyr::mutate(size=18000)
# merge data frames
cov_data <- dplyr::bind_rows(cov_data_6000_nonrct,cov_data_12000_nonrct,
                             cov_data_18000_nonrct)

# ── Marginal coverage factor computation   ────────────────────────────────────
cov_factor_data <- cov_data %>%
  dplyr::mutate(
    type = factor(type),         # factor for discrete color
    mechanism = factor(mechanism)) %>%
  dplyr::mutate(cov_factor = cov_data$value - (1 - level)) %>%
  dplyr::mutate(across(-c(value,cov_factor ), as.factor))%>%
  dplyr::rename("cov_mean"=value)

# ── Create complete data frame  ───────────────────────────────────────────────
complete_data <- list(spv_means, mean_cardinality_data,cov_factor_data) %>%
  purrr::reduce(full_join, by = c("mechanism","level", "type", "size"))

plot_data <- complete_data %>%
  dplyr::mutate(color_group = case_when(
    mechanism == "Estimated labels" ~ paste0("type_", type),
    mechanism == "Oracular CP" ~ "Oracular CP",
    mechanism == "GLB" ~ "GLB"
  ))

type_vals <- sort(unique(plot_data$type))

# ── Create SPV-level plots ────────────────────────────────────────────────────
plot_spv_level <- ggplot2::ggplot(plot_data,
                                  ggplot2::aes(x = level, y = mean_spv)) +
  ggplot2::geom_line(data = subset(plot_data, mechanism == "Estimated labels"),
                     ggplot2::aes(group = type, color = color_group),
                     alpha = 0.7) +
  ggplot2::geom_point(aes(size = mean_cardinality, color = color_group),
                      alpha = 0.5, show.legend = c(size=FALSE)) +
  ggplot2::scale_color_manual(
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
    labels = c(paste0("r = ", type_vals), "Oracular CP", "GLB")) +
  ggplot2::scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  ggplot2::facet_grid(~size) +
  ggplot2::labs(x = expression("Confidence level (" * alpha * ")"),
       y = "Set Policy Value (SPV)")
ggplot2::ggsave(plot_spv_level, filename=paste0("inst/images/Level_SPV_", type,".pdf"), width = 15, height = 8)

plot_cov_level <- ggplot2::ggplot(plot_data,
                                  ggplot2::aes(x = level,
                                  y = cov_mean)) +
  ggplot2::geom_line(data = subset(plot_data, mechanism == "Estimated labels"),
                     ggplot2::aes(group = type, color = color_group),
                     alpha = 0.7) +
  ggplot2::geom_point(ggplot2::aes(size = mean_cardinality, color = color_group),
             alpha = 0.5, show.legend = c(size=FALSE)) +
  ggplot2::scale_color_manual(
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
  ggplot2::scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  ggplot2::facet_grid( ~ size) +
  ggplot2::labs(x = expression("Confidence level (" * alpha * ")"),
       y = expression("Coverage attained"))
ggplot2::ggsave(plot_cov_level, filename=paste0("inst/images/Level_Coverage_", type,".pdf"), width = 15, height = 8)

plot_cov_factor_level <- ggplot2::ggplot(plot_data,
                                         ggplot2::aes(x = level, y = cov_factor,
                                  group = color_group, color=color_group)) +
  ggplot2::geom_line(alpha = 0.7) +
  ggplot2::geom_point(ggplot2::aes(size = mean_spv, color = color_group),
             alpha = 0.5, show.legend = c(size=FALSE)) +
  ggplot2::geom_hline(yintercept = 0, color="red")+
  ggplot2::scale_color_manual(
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
  ggplot2::scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  ggplot2::facet_grid(~ size) +
  ggplot2::labs(x = expression("Confidence level ("* alpha *")"),
       y = "Marginal coverage factor")
ggplot2::ggsave(plot_cov_factor_level, filename=paste0("inst/images/Level_cov_factor_", type,".pdf"), width = 15, height = 8)

plot_width_spv <- ggplot2::ggplot(plot_data,
                                  ggplot2::aes(x = mean_cardinality,
                                  y = mean_spv)) +

  ggplot2::geom_line(data = subset(plot_data, mechanism == "Estimated labels"),
                     ggplot2::aes(group = type, color = color_group),
            alpha = 0.7) +
  ggplot2::geom_point(ggplot2::aes(size = mean_cardinality, color = color_group),
             alpha = 0.5, show.legend = c(size=FALSE)) +
  ggplot2::scale_color_manual(
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
  ggplot2::scale_size_continuous(name = "Mean cardinality",
                        range = c(0.1, 2)) +
  ggplot2::facet_grid(~size) +
  ggplot2::labs(x = "Mean cardinality",
       y = "Set Policy Value (SPV)")
ggplot2::ggsave(plot_width_spv, filename=paste0("inst/images/mean_cardinality_SPV_", type,".pdf"), width = 15, height = 8)

plot_mean_level<- ggplot2::ggplot(plot_data,
                                  ggplot2::aes(x = level,
                                  y = mean_cardinality,
                                  color = color_group)) +
  ggplot2::geom_line(ggplot2::aes(group = color_group),
            alpha = 0.7) +
  ggplot2::geom_point(ggplot2::aes(color = color_group),alpha = 0.5) +
  ggplot2::scale_color_manual(
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
  ggplot2::facet_grid(~size) +
  ggplot2::labs(x = expression("Confidence level ("* alpha *")"),
       y = "Mean cardinality")
ggplot2::ggsave(plot_mean_level, filename=paste0("inst/images/mean_cardinality_level_", type,".pdf"), width = 15, height = 8)


optimal_treatments <- function(df) {
    df <- as.matrix(df)
    mat <- matrix(0, nrow = nrow(df), ncol = ncov)
    mat[, 1:2] <- df
    if(type=="normal"){
      p_o <- mu_P0_normal(mat)
    }else{
      p_o <- mu_P0_simplex_complicated(mat)
    }
    apply(data.frame(1:nrow(mat)), 1, function(i){
      paste0("{",
             paste(which(p_o[i,]==max(p_o[i,])), collapse = ","), 
             "}")})
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

ggplot2::ggsave(plot_sythetic_scenario, 
                filename=paste0("inst/images/Synthetic_data_", type,".pdf"), 
                width = 10, height = 6)
