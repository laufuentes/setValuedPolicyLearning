add_r <- TRUE
idx_r <- 1:2
r_value <- ifelse(type == "tree", 0.15, 0.1)

# ── Generate the rownames for latex table ─────────────────────────────────────
types_optimal_treatment <- SL.out$optimal_policy_new |> unique()
obs_types <- sapply(SL.out$optimal_policy_new, function(x) {
  which(sapply(types_optimal_treatment, function(y) setequal(x, y)))
})
treatment_labels <- sapply(types_optimal_treatment, function(x) {
  paste0("$\\Pi^{\\star}(X_{i})$ = \\{", paste(x, collapse = ", "), "\\}")
})

cov.treatment_labels <- sapply(types_optimal_treatment, function(x) {
  paste0("Strict coverage \\{", paste(x, collapse = ", "), "\\}")
})
obs_type_names <- treatment_labels[obs_types]
obs_cov_names <- cov.treatment_labels[obs_types]

# ── Extract results from results_list ─────────────────────────────────────────
### Exact match and coverage
exact_match_all <-  coverage_all <- array(0, dim = c(nrow(results_list[[1]]$exact_match), 
                                    ncol(results_list[[1]]$exact_match), 
                                    length(results_list)))

for (i in seq_len(dim(coverage_all)[3])) {
  exact_match_all[, , i] <- results_list[[i]]$exact_match
  coverage_all[, , i] <- results_list[[i]]$coverage
}

exact_match_mean <- apply(exact_match_all, MARGIN = c(1, 2), FUN = mean)

df_exact_match <- data.frame(
  Set_Type = obs_type_names, exact_match_mean)|>
  group_by(Set_Type) |>
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))

coverage_mean    <- apply(coverage_all, MARGIN = c(1, 2), FUN = mean)

df_strict_coverage <- data.frame(
  Set_Type = obs_cov_names, coverage_mean)|>
  group_by(Set_Type) |>
  summarise(across(starts_with("X"), ~ mean(.x, na.rm = TRUE)))

### Other features (cardinality, SPVs, marginal and total strict coverage)
cardinality_all <- do.call(rbind, 
                           lapply(results_list, `[[`, "cardinality"))
mean_card    <- colMeans(cardinality_all)
strict_cov   <- colMeans(coverage_mean)
relaxed_cov  <- colMeans(do.call(rbind, 
                                 lapply(results_list, 
                                        `[[`, "relaxed_coverage"))) # marginal coverage
spv_unif_all <- do.call(rbind, 
                        lapply(results_list, `[[`, "spv_uniform")) 
spv_unif     <- colMeans(spv_unif_all)

spv_prop_all <- do.call(rbind, 
                        lapply(results_list, 
                               `[[`, "spv_propensity"))
spv_prop     <- colMeans(spv_prop_all)

### Names of baseline policies and conformal set-valued policies
methods_pl <- results_list[[1]]$selected_methods
clean_methods <- grep(paste0("\\(r=0\\.",ifelse(type=="tree",15,1),"\\)$"), methods_pl, value = TRUE, invert = TRUE)

# Ordered names
dynamic_methods <- unlist(lapply(clean_methods, function(m) {
  if (add_r && m %in% clean_methods[idx_r]) {
    c(m, paste0("Conformal ", m), paste0("Conformal ", m, " (r=",r_value,")"))
  } else {
    c(m, paste0("Conformal ", m))
  }})) 
order_elements <- c("Set Type", "GLB GRF", "GLB GLM", dynamic_methods, 
                    "Policy Aggregation", "Conformal Aggregation")

# ── Create latex table ────────────────────────────────────────────────────────
# Group results
df_summary <- rbind(
  df_exact_match, c("Mean cardinality", mean_card),
  df_strict_coverage, c("Total strict coverage", strict_cov),
  c("Marginal coverage", relaxed_cov),
  c("Uniform SPV", spv_unif), c("Propensity SPV", spv_prop))

colnames(df_summary)<- c(
  "Set Type", "GLB GRF", "GLB GLM", 
  "Policy Aggregation", "Conformal Aggregation", 
  clean_methods, 
  paste0("Conformal ", methods_pl))

df_summary <- df_summary |>
  select(all_of(order_elements)) |>
  mutate(across(-`Set Type`, ~ round(as.numeric(.x), 2)))

# Construct sub-headers vector
dynamic_cols <- lapply(clean_methods, function(x) {
  if (add_r && x %in% clean_methods[idx_r]) {
    c("Policy", "Conf.", paste0("Conf.$_{r=", r_value, "}$"))
  } else {
    c("Policy", "Conf.")
  }}) |> unlist()

sub_headers <- c("\\textbf{Set Type}", "GRF", "GLM", dynamic_cols, "Policy", "Conf.")

# Construct group headers structure 
header_groups <- c(
  " " = 1,
  "\\\\textbf{GLB}" = 2,
  setNames(
    if (add_r) ifelse(clean_methods %in% clean_methods[idx_r], 3, 2) else 2,
    paste0("\\\\textbf{", clean_methods, "}")
  ),
  "\\\\textbf{Aggregation}" = 2
)

latex_tbl <- kable(
  df_summary,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  col.names = sub_headers,
  align = c("l", rep("c", ncol(df_summary) - 1)),
  label = "tab:exact_matches_mean"
) %>%
  add_header_above(header_groups, escape = FALSE) %>% 
  kable_styling(latex_options = c("HOLD_position"), font_size = 9)

writeLines(
  c(
    "% latex table generated in R",
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    "\\setlength{\\tabcolsep}{1pt}",
    as.character(latex_tbl),
    "\\end{table}"
  ),
  con = file.path(subdir_path,"results.txt")
)


if(type =="linear"){
  # ── Create latex table for reduced linear scenario ──────────────────────────
  target_methods <- c("MACF", "ql.lm")
  sub_clean_methods <- intersect(clean_methods, target_methods)
  
  normal_dynamic_methods <- unlist(lapply(sub_clean_methods, function(m) {
    if (add_r && m %in% clean_methods[idx_r]) {
      c(m, paste0("Conformal ", m), paste0("Conformal ", m, paste0(" (r=",r_value,")")))
    } else {
      c(m, paste0("Conformal ", m))
    }
  }))
  
  target_cols <- c("Set Type", "GLB GRF", "GLB GLM", normal_dynamic_methods)
  
  # Subset rows and columns
  df_summary_reduced <- df_summary %>%
    filter(!`Set Type` %in% c("Uniform SPV", "Propensity SPV")) %>%
    select(all_of(target_cols))
  
  # Construct sub-headers vector
  dynamic_cols_reduced <- lapply(sub_clean_methods, function(x) {
    if (add_r && x %in% clean_methods[idx_r]) {
      c("Policy", "Conf.", paste0(" Conf.$_{r=", r_value, "}$"))
    } else {
      c("Policy", "Conf.")}
    }) |> unlist()
  sub_headers_reduced <- c("\\textbf{Set Type}", "GRF", "GLM", dynamic_cols_reduced)
  
  # Construct group headers structure 
  header_groups_reduced <- c(
    " " = 1,
    "\\textbf{GLB}" = 2,
    setNames(
      sapply(sub_clean_methods, function(m) {
        if (add_r && m %in% clean_methods[idx_r]) 3 else 2
      }),
      paste0("\\textbf{", sub_clean_methods, "}")
    )
  )
  
  latex_tbl_reduced <- kable(
    df_summary_reduced,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    col.names = sub_headers_reduced,
    align = c("l", rep("c", ncol(df_summary_reduced) - 1)),
    label = "tab:exact_matches_mean") %>%
    add_header_above(header_groups_reduced, escape = FALSE) %>%
    kable_styling(latex_options = c("HOLD_position"), font_size = 9)
  
  writeLines(
    c(
      "% latex table generated in R",
      "\\begin{table}[H]",
      "\\centering",
      "\\small",
      "\\setlength{\\tabcolsep}{1pt}",
      as.character(latex_tbl_reduced),
      "\\end{table}"
    ),
    con = file.path(subdir_path,"results_reduced.txt")
  )
}

# ── Create SPV boxplot figure  ────────────────────────────────────────────────
colnames(spv_prop_all) <- colnames(spv_unif_all) <- colnames(cardinality_all) <- 
  c("GLB GRF", "GLB GLM", 
    "Policy Aggregation", "Conformal Aggregation", 
    clean_methods, 
    paste0("Conformal ", methods_pl))

# Uniform SPVs
df_spv.unif <- spv_unif_all |> 
  as.data.frame()|>
  select(-all_of(c(clean_methods, "Policy Aggregation"))) |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "SPV Value") |> 
  mutate(Metric = "Uniform SPV")

# Policy values 
df_pv <- spv_prop_all |> 
  as.data.frame() |>
  select(all_of(c(clean_methods, "Policy Aggregation"))) |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "SPV Value") |> 
  mutate(Metric = "Policy value")

# Propensity SPVs
df_spv.prop <- spv_prop_all |> 
  as.data.frame() |>
  select(-all_of(c(clean_methods, "Policy Aggregation"))) |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "SPV Value") |> 
  mutate(Metric = "Propensity SPV")

# Horizontal lines
potential.outcomes <- SL.out$potential_outcomes
opt.treatment <- do.call(rbind, SL.out$optimal_policy_new)
opt.value <- SL.out$potential_outcomes[cbind(1:nrow(potential.outcomes), 
                                             opt.treatment[,1])] %>% mean()
levels_A <- levels(SL.out$df_new_sample$A)
m <- length(levels_A)
A_rd <- sample(as.numeric(levels_A), size = nrow(SL.out$df_new_sample), replace = TRUE)
random.value <- SL.out$potential_outcomes[cbind(1:nrow(potential.outcomes), A_rd)] %>% mean()

# Create figure
spv_boxplot <- bind_rows(df_spv.unif, df_spv.prop, df_pv) |> 
  dplyr::mutate(Method = factor(`Method`, levels = order_elements)) |> 
  ggplot2::ggplot(
    ggplot2::aes(
      x = `Method`, 
      y = `SPV Value`, 
      color = `Metric`)) +
  ggplot2::geom_boxplot(na.rm = TRUE)+
  ggplot2::geom_hline(ggplot2::aes(yintercept = opt.value), colour = "black")+
  ggplot2::geom_hline(ggplot2::aes(yintercept = random.value), colour = "red")+ 
  ggplot2::theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggplot2::ggsave(spv_boxplot, filename = file.path(subdir_path,"SPV_boxplot.pdf"), width = 18, height = 6)

if(type=="linear"){
  # Reduced figure for linear scenario
  spv_boxplot <- bind_rows(df_spv.unif, df_spv.prop, df_pv) |> 
    filter(!grepl("drql.lm|aggregation|drql.ksvm", Method, ignore.case = TRUE))|>
    dplyr::mutate(Method = factor(`Method`, levels = order_elements)) |> 
    ggplot2::ggplot(
      ggplot2::aes(
        x = `Method`, 
        y = `SPV Value`, 
        color = `Metric`)) +
    ggplot2::geom_boxplot(na.rm = TRUE)+
    ggplot2::geom_hline(ggplot2::aes(yintercept = opt.value), colour = "black")+
    ggplot2::geom_hline(ggplot2::aes(yintercept = random.value), colour = "red")+ 
    ggplot2::theme()
  
  ggplot2::ggsave(spv_boxplot, 
                  filename = file.path(subdir_path,"SPV_boxplot_reduced.pdf"), width = 12, height = 6)
  
}
# ── Create cardinality boxplot figure ─────────────────────────────────────────
cardinality_boxplot <- cardinality_all |> 
  as.data.frame() |>
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "Cardinality") |> 
  dplyr::mutate(Method = factor(`Method`, levels = order_elements)) |> 
  ggplot2::ggplot(ggplot2::aes(x=`Method`, y = `Cardinality`)) +
  ggplot2::geom_boxplot()
  
ggplot2::ggsave(cardinality_boxplot, 
                filename = file.path(subdir_path,"Cardinality_boxplot.pdf"), 
                width = 18, height = 6)

# ── Create Baseline policy values figure ──────────────────────────────────────
# Compute oracular policy values for all baselines
values<- apply(data.frame(1:30), 1, function(i){
  dopt <- results_list[[i]]$doptFactorPredict_new_naive
  apply(data.frame(1:ncol(dopt)), 1, function(j){
    potential.outcomes[cbind(1:nrow(potential.outcomes), dopt[,j])] |> mean()})})

numalgs <- ncol(results_list[[1]]$doptFactorPredict_new_naive)
# (including the aggregation of all baselines)
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

ggplot2::ggsave( filename= file.path(subdir_path,"Naive_policy_values.pdf"),
                 width = 10, height = 6)

# ── Create heatmap for introduction ───────────────────────────────────────────
true_df <- as.data.frame(do.call(rbind, SL.out$optimal_policy_new))|>
  mutate(Row = row_number()) |>
  pivot_longer(
    cols = -Row, 
    values_to = "Levels"
  ) |>
  mutate(Levels = factor(Levels, levels=levels_A))

# Ground truth column
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

# Baseline policies column
heatmap_data <- results_list[[1]]$doptFactorPredict_new_naive |> 
  as.data.frame() |>
  mutate(Row = row_number()) |>
  pivot_longer(
    cols = -Row,
    names_to = "Policy learning method",
    values_to = "Treatments"
  ) |>
  mutate(Treatments = factor(Treatments)) # , levels = levels_A

heatmap_pl <- ggplot(heatmap_data |> 
                       filter(`Row`<=10), aes(x = `Policy learning method`, 
                                              y = `Row`, fill = `Treatments`)) +
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
ggplot2::ggsave(heatmap_both, 
                filename = file.path(subdir_path,"heatmap_recommendations.pdf"),
                width = 10, height = 8)

# ── Create synthetic data plot  ───────────────────────────────────────────────
optimal_treatments <- function(df) {
  df <- as.matrix(df)
  mat <- matrix(0, nrow = nrow(df), ncol = 5)
  mat[, 1:2] <- df
  if(type=="linear"){
    p_o <- mu_P_linear(mat)
  }else if(type=="complex"){
    p_o <- mu_P_complex(mat)
  }else{
    p_o <- mu_P_tree(mat)
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
                filename= file.path(subdir_path,
                                    paste0("Synthetic_data_", type, ".pdf")),
                width = 10, height = 6)


