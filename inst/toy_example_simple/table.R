# ── Generate the rownames  ─────────────────────────────────────────────────────────
types_optimal_treatment <- SL.out$optimal_policy_new |> unique()
obs_types <- sapply(SL.out$optimal_policy_new, function(x) {
  which(sapply(types_optimal_treatment, function(y) setequal(x, y)))
})
treatment_labels <- sapply(types_optimal_treatment, function(x) {
  paste0("$\\Pi^{\\star}(X_{i})$ = \\{", paste(x, collapse = ", "), "\\}")
})

cov.treatment_labels <- sapply(types_optimal_treatment, function(x) {
  paste0("Conditional coverage \\{", paste(x, collapse = ", "), "\\}")
})

obs_type_names <- treatment_labels[obs_types]
obs_cov_names <- cov.treatment_labels[obs_types]

# Aggregate scores by prediction type
exact_match_all<- array(0, dim = c(nrow(results_list[[1]]$exact_match), 
                                    ncol(results_list[[1]]$exact_match), 
                                    length(results_list)))

coverage_all  <- array(0, dim = c(nrow(results_list[[1]]$exact_match), 
                                    ncol(results_list[[1]]$exact_match), 
                                    length(results_list)))
for (i in seq_len(dim(coverage_all)[3])) {
  exact_match_all[, , i] <- results_list[[i]]$exact_match
  coverage_all[, , i] <- results_list[[i]]$coverage
}

exact_match_mean <- apply(exact_match_all, MARGIN = c(1, 2), FUN = mean)
coverage_mean    <- apply(coverage_all, MARGIN = c(1, 2), FUN = mean)
total.coverage <- colMeans(coverage_mean)

cardinality_all    <- do.call(rbind, lapply(results_list, `[[`, "cardinality"))
spv_uniform_all    <- do.call(rbind, lapply(results_list, `[[`, "spv_uniform"))
spv_propensity_all <- do.call(rbind, lapply(results_list, `[[`, "spv_propensity"))

# Prepare tables for latex 
eval_by_type <- data.frame(
  Set_Type = obs_type_names,
  exact_match_mean= exact_match_mean)

df_exact_match <- eval_by_type |>
  group_by(Set_Type) |>
  summarise(across(starts_with("exact_match"), ~ mean(.x, na.rm = TRUE)))

methods_pl <- results_list[[1]]$selected_methods
colnames(df_exact_match)<- c("Set Type","GLB GRF", "GLB GLM", 
                             "Policy Aggregation", "Conformal Aggregation", 
                             methods_pl, paste0("Conformal ", methods_pl))


eval_by_type.cov <- data.frame(
  Set_Type = obs_cov_names,
  cov.strict = coverage_mean)

df_strict_coverage <- eval_by_type.cov |>
  group_by(Set_Type) |>
  summarise(across(starts_with("cov."), ~ mean(.x, na.rm = TRUE)))

colnames(df_strict_coverage)<- c("Set Type","GLB GRF", "GLB GLM", 
                                 "Policy Aggregation", "Conformal Aggregation", 
                                 methods_pl, paste0("Conformal ", methods_pl))

order_elements <- c(
  "Set Type", 
  "GLB GRF", 
  "GLB GLM", 
  c(rbind(methods_pl, paste0("Conformal ", methods_pl))), # Interleaves method[i] and Conformal method[i]
  "Policy Aggregation", 
  "Conformal Aggregation"
)

df_summary <- rbind(df_exact_match,
                    c("Mean cardinality", colMeans(cardinality_all)), 
                    df_strict_coverage, 
                    c("Total coverage",total.coverage),
                    c("Uniform SPV", colMeans(spv_uniform_all)), 
                    c("Propensity SPV", colMeans(spv_propensity_all)))|>
  mutate(across(-`Set Type`, as.numeric)) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))|> 
  select(all_of(order_elements))


# ── Create latex table ────────────────────────────────────────────────────────
library(xtable)

x_tab <- xtable(
  df_summary,
  label = "tab:exact_matches_mean",
  digits = c(0, 0, 3, 3, 3, 3, 3,3) # First digit is for row indices, rest for columns
)

print(
  x_tab,
  include.rownames = FALSE,
  booktabs = TRUE,          
  caption.placement = "top",
  sanitize.text.function = identity,
  file = "inst/toy_example_simple/images/results.txt"
)

# ── Create SPV boxplot ────────────────────────────────────────────────────────
colnames(spv_uniform_all) <- colnames(spv_propensity_all) <- colnames(cardinality_all) <- 
  c("GLB GRF", "GLB GLM",  "Policy Aggregation", 
    "Conformal Aggregation", methods_pl, paste0("Conformal ", methods_pl))

df_spv.unif <- spv_uniform_all |> 
  as.data.frame()|>
  select(-all_of(c(methods_pl, "Policy Aggregation"))) |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "SPV Value") |> 
  mutate(Metric = "Uniform SPV")


df_pv <- spv_propensity_all |> 
  as.data.frame() |>
  select(all_of(c(methods_pl, "Policy Aggregation"))) |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "SPV Value") |> 
  mutate(Metric = "Policy value")

df_spv.prop <- spv_propensity_all |> 
  as.data.frame() |>
  select(-all_of(c(methods_pl, "Policy Aggregation"))) |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "SPV Value") |> 
  mutate(Metric = "Propensity SPV")


potential.outcomes <- SL.out$potential_outcomes
opt.treatment <- do.call(rbind, SL.out$optimal_policy_new)
opt.value <- SL.out$potential_outcomes[cbind(1:nrow(potential.outcomes), 
                                             opt.treatment[,1])] %>% mean()

A_rd <- sample(as.numeric(levels_A), size = nrow(SL.out$df_new_sample), replace = TRUE)
random.value <- SL.out$potential_outcomes[cbind(1:nrow(potential.outcomes), A_rd)] %>% mean()

spv_boxplot <- bind_rows(df_spv.unif, df_spv.prop, df_pv) |> 
  dplyr::mutate(Method = factor(`Method`, levels = order_elements)) |> 
  ggplot2::ggplot(
    ggplot2::aes(
      x = `Method`, 
      y = `SPV Value`, 
      color = `Metric`)) +
  ggplot2::geom_boxplot(na.rm = TRUE)+
  ggplot2::geom_hline(ggplot2::aes(yintercept = opt.value), colour = "black")+
  ggplot2::geom_hline(ggplot2::aes(yintercept = random.value), colour = "red")

ggplot2::ggsave(spv_boxplot, filename = "inst/toy_example_simple/images/SPV_boxplot.pdf", width = 18, height = 6)

# ── Create cardinality boxplot ────────────────────────────────────────────────
cardinality_boxplot <- cardinality_all |> 
  as.data.frame() |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "Cardinality") |> 
  dplyr::mutate(Method = factor(`Method`, levels = order_elements)) |> 
  ggplot2::ggplot(ggplot2::aes(x=`Method`, y = `Cardinality`)) +
  ggplot2::geom_boxplot()
  
ggplot2::ggsave(cardinality_boxplot, filename = "inst/toy_example_simple/images/Cardinality_boxplot.pdf", width = 18, height = 6)

# ── Create naive policy values boxplot ────────────────────────────────────────
values<- apply(data.frame(1:30), 1, function(i){
  dopt <- results_list[[i]]$doptFactorPredict_new_naive
  apply(data.frame(1:ncol(dopt)), 1, function(j){
    potential.outcomes[cbind(1:nrow(potential.outcomes), dopt[,j])] |> mean()})})

numalgs <- ncol(results_list[[1]]$doptFactorPredict_new_naive)
agg_value<- apply(data.frame(1:30), 1, function(i){
  dopt <- results_list[[i]]$doptFactorPredict_new_naive
  unweighted_probs <- weighted_probs_experts(fitted_experts = dopt,weights =rep(1/numalgs, numalgs),
                                             df_pred = SL.out$df_new_sample, levels = as.numeric(levels_A))
  
  unweighted_aggregation <- apply(apply(unweighted_probs, 1, function(x){
    rmultinom(1,1,prob=x)}), 2, which.max)
  potential.outcomes[cbind(1:nrow(potential.outcomes), unweighted_aggregation)] |> mean()})


values <- as.data.frame(t(rbind(values,agg_value)))
colnames(values) <- c(colnames(results_list[[1]]$doptFactorPredict_new_naive),"Aggregation")

values |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "Policy value") |> 
  ggplot2::ggplot(ggplot2::aes(x = `Method`, y = `Policy value`)) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_hline(ggplot2::aes(yintercept = opt.value), colour = "black")+
  ggplot2::geom_hline(ggplot2::aes(yintercept = random.value), colour = "red")


ggplot2::ggsave( filename=paste0("inst/toy_example_simple/images/Naive_policy_values.pdf"), 
                 width = 10, height = 6)


# ── Create heatmap for motivation ─────────────────────────────────────────────
levels_A <- levels(SL.out$df_new_sample$A)
m <- length(levels_A)
true_df <- as.data.frame(do.call(rbind, SL.out$optimal_policy_new))|>
  mutate(Row = row_number()) |>
  pivot_longer(
    cols = -Row, 
    values_to = "Levels"
  ) |>
  mutate(Levels = factor(Levels, levels=levels_A))

true_heatmap <- ggplot(true_df
                       |> filter(`Row`<=10), 
                       aes(x = name, y = `Row`, fill = `Levels`)) +
  geom_tile(linewidth = 0.1)+
  scale_fill_viridis_d(option = "viridis", drop = FALSE) +
  scale_color_discrete()+
  scale_y_reverse() +
  xlab("Ground truth")+
  theme(
    axis.text.y = element_blank(),
    legend.position = "none", 
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = ggplot2::margin(b = 44, r=-5)
  )

heatmap_data <- results_list[[1]]$doptFactorPredict_new_naive |> 
  as.data.frame() |>
  mutate(Row = row_number()) |>
  pivot_longer(
    cols = -Row,
    names_to = "Policy learning method",
    values_to = "Levels"
  ) |>
  mutate(Levels = factor(Levels)) # , levels = levels_A

heatmap_pl <- ggplot(heatmap_data |> 
                       filter(`Row`<=10), aes(x = `Policy learning method`, 
                                              y = `Row`, fill = `Levels`)) +
  geom_tile(linewidth = 0.1)+
  scale_fill_manual(values = viridisLite::viridis(m), 
                    drop=TRUE, 
                    limits = levels_A) +
  scale_y_reverse() +
  ylab("")+
  theme(
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    plot.margin =  ggplot2::margin(l=-10)
  )

heatmap_both <- grid.arrange(true_heatmap, heatmap_pl, ncol=2, widths = c(1, 4))
ggplot2::ggsave(heatmap_both, filename = paste0("inst/toy_example_simple/images/heatmap_recommendations.pdf"), 
                width = 10, height = 8)

# ── Create synthetic data plot  ───────────────────────────────────────────────
optimal_treatments <- function(df) {
  df <- as.matrix(df)
  mat <- matrix(0, nrow = nrow(df), ncol = 5)
  mat[, 1:2] <- df
  if(type=="normal"){
    p_o <- mu_P0_normal(mat)
  }else{
    p_o <- mu_P0_complex(mat)
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
    fill = "Optimal treatments") + 
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 18),
    legend.title = ggplot2::element_text(size = 18), 
    legend.text = ggplot2::element_text(size = 16), 
    strip.text = ggplot2::element_text(size = 16), 
    legend.key.size = ggplot2::unit(1, "cm"))

ggplot2::ggsave(plot_sythetic_scenario, 
                filename=paste0("inst/toy_example_simple/images/Synthetic_data_", type,".pdf"), 
                width = 10, height = 6)


# ── Create colored naive method boxplot  ──────────────────────────────────────
library(scales)
long_df <- cbind(
  SL.out$df_new_sample[, c("X1", "X2")],
  results_list[[1]]$doptFactorPredict_new_naive) |>
  pivot_longer(
    cols = -c(X1, X2),
    names_to = "Policy",
    values_to = "Value"
  )

default_2_colors <- hue_pal()(2)

ggplot(long_df, aes(x = X1, y = X2, color = factor(Value))) +
  geom_point() +
  scale_color_manual(
    values = c(
      "1" = default_2_colors[1], 
      "2" = default_2_colors[1], 
      "3" = default_2_colors[2], 
      "4" = "green",
    name = "Category")) +
  facet_wrap(~ Policy)

ggsave(filename = "inst/toy_example_simple/images/naive_colors.pdf")
