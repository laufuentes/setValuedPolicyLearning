setwd("~/Documents/PhD/Project 2 - Conformal Policy Sets /CPolicySets")

library(SL.ODTR)
library(hitandrun)
library(tidyr)
library(dplyr)
library(lava)
library(purrr)
library(grf)
library(randomForest)
library(gridExtra)
library(SuperLearner)
library(policytree)
library(glmnet)
library(tmle)
library(parallel)
library(caret)
library(polle)
library(viridisLite)

source("src/synthetic_data.R")
source("src/non-conformity-scores.R")
source("src/label-estimation.R")
source("src/utils.R")
source("src/evaluation.R")

score_name <- "margin/" # type of non conformity score

synthetic_scenario <- TRUE # if synthetic data, FALSE otherwise  
type <- "normal" # additional name for images (here: type of synthetic scenario)

score_name <- "margin/"
random_rate <- seq(0,1,0.1) # random rates to test 
n_rate <- length(random_rate) # number of random rates to test 


n_test <- 100 # number of repetitions for boxplots
alphas <- seq(0,1,0.05) # number of confidence levels to test 
covariates_name <- c("x1","x2") # name for covariates in dataset
treatment_name <- "A" # name of treatment indicator in dataset
outcome_name <- "Y" # name of outcome in dataset

# RCT data 
results_6000_rct <- readRDS(paste0("experts_pred/", score_name, "RCT/", "plot_results_", type, "_6000", ".rds"))
SL.out_6000_rct <- readRDS(paste0("experts_pred/", score_name, "RCT/", "_", type, "_6000.rds"))
results_12000_rct <- readRDS(paste0("experts_pred/", score_name, "RCT/", "plot_results_", type, "_12000", ".rds"))
SL.out_12000_rct <- readRDS(paste0("experts_pred/", score_name, "RCT/", "_", type, "_12000.rds"))
results_18000_rct <- readRDS(paste0("experts_pred/", score_name, "RCT/", "plot_results_", type, "_18000", ".rds"))
SL.out_18000_rct <- readRDS(paste0("experts_pred/", score_name, "RCT/", "_", type, "_18000.rds"))

# non-RCT data 
results_6000_nonrct <- readRDS(paste0("experts_pred/", score_name, "non_RCT/", "plot_results_", type, "_6000", ".rds"))
SL.out_6000_nonrct <- readRDS(paste0("experts_pred/", score_name, "non_RCT/", "_", type, "_6000.rds"))
results_12000_nonrct <- readRDS(paste0("experts_pred/", score_name, "non_RCT/", "plot_results_", type, "_12000", ".rds"))
SL.out_12000_nonrct <- readRDS(paste0("experts_pred/", score_name, "non_RCT/", "_", type, "_12000.rds"))
results_18000_nonrct <- readRDS(paste0("experts_pred/", score_name, "non_RCT/", "plot_results_", type, "_18000", ".rds"))
SL.out_18000_nonrct <- readRDS(paste0("experts_pred/", score_name, "non_RCT/", "_", type, "_18000.rds"))

# spv 
spv_data_6000_rct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_6000_rct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_6000_rct[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000_rct[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_6000_rct[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000_rct[["spv"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size =6000, study="RCT")

spv_data_6000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_6000_nonrct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_6000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_6000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["spv"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size =6000, study="observational")

spv_data_12000_rct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_12000_rct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_12000_rct[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000_rct[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_12000_rct[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000_rct[["spv"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 12000, study="RCT") 

spv_data_12000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_12000_nonrct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_12000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_12000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["spv"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 12000, study="observational") 


spv_data_18000_rct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_18000_rct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_18000_rct[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000_rct[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_18000_rct[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000_rct[["spv"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 18000, study="RCT")

spv_data_18000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_18000_nonrct[["spv"]], alphas, random_rate), 
  map_dfr(1:dim(results_18000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["spv"]][a, ,4,1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}), 
  map_dfr(1:dim(results_18000_nonrct[["spv"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["spv"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1]))}))%>% mutate(size = 18000, study="observational")

# mean width
mean_width_data_6000_rct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", results_6000_rct[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_6000_rct[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() ,
  data.frame(
    value = results_6000_rct[["mean_width"]][, 5, 1],
    mechanism = "True labels",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=6000, study="RCT")

mean_width_data_6000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", results_6000_nonrct[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_6000_nonrct[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() ,
  data.frame(
    value = results_6000_nonrct[["mean_width"]][, 5, 1],
    mechanism = "True labels",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=6000, study="observational")

mean_width_data_12000_rct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", results_12000_rct[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_12000_rct[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() , 
  data.frame(
    value = results_12000_rct[["mean_width"]][, 5, 1],
    mechanism = "True labels",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=12000, study="RCT")

mean_width_data_12000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", results_12000_nonrct[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_12000_nonrct[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() , 
  data.frame(
    value = results_12000_nonrct[["mean_width"]][, 5, 1],
    mechanism = "True labels",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=12000, study="observational")

mean_width_data_18000_rct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", results_18000_rct[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_18000_rct[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup()
  ,data.frame(
    value = results_18000_rct[["mean_width"]][, 5, 1],
    mechanism = "True labels",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=18000, study="RCT")

mean_width_data_18000_nonrct <- dplyr::bind_rows(
  make_smaller_block(1, "Estimated labels", results_18000_nonrct[["mean_width"]], random_rate) %>% group_by(mechanism,type)%>% mutate(levels=(row_number()-1)/(length(alphas)-1)) %>% ungroup(), 
  data.frame(
    value = results_18000_nonrct[["mean_width"]][, 4, 1],
    mechanism = "ULB",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup()
  ,data.frame(
    value = results_18000_nonrct[["mean_width"]][, 5, 1],
    mechanism = "True labels",
    type = paste0(random_rate[1])
  ) %>% dplyr::group_by(mechanism,type)%>%
    dplyr::mutate(levels=(row_number()-1)/(length(alphas)-1))%>% 
    dplyr::ungroup() 
) %>% mutate(size=18000, study="observational")

# Coverage plot
cov_data_6000_rct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_6000_rct[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_6000_rct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000_rct[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}), 
  map_dfr(1:dim(results_6000_rct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000_rct[["cov_unif"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=6000, study="RCT")

cov_data_6000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_6000_rct[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_6000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}), 
  map_dfr(1:dim(results_6000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_6000_nonrct[["cov_unif"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=6000, study="observational")

cov_data_12000_rct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_12000_rct[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_12000_rct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000_rct[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}),
  map_dfr(1:dim(results_12000_rct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000_rct[["cov_unif"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=12000, study="RCT")

cov_data_12000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_12000_nonrct[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_12000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}),
  map_dfr(1:dim(results_12000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_12000_nonrct[["cov_unif"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=12000, study="observational")

cov_data_18000_rct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_18000_rct[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_18000_rct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000_rct[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}),
  map_dfr(1:dim(results_18000_rct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000_rct[["cov_unif"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=18000, study="RCT")

cov_data_18000_nonrct <- dplyr::bind_rows(
  make_block(1, "Estimated labels", results_18000_nonrct[["cov_unif"]], alphas, random_rate), 
  map_dfr(1:dim(results_18000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["cov_unif"]][a, , 4, 1],
      mechanism = "ULB",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )}),
  map_dfr(1:dim(results_18000_nonrct[["cov_unif"]])[1], function(a) {
    data.frame(
      value = results_18000_nonrct[["cov_unif"]][a, , 5, 1],
      mechanism = "True labels",
      level = paste0(alphas[a]),
      type = paste0(random_rate[1])
    )})
) %>% mutate(size=18000, study="observational")


cov_data <- bind_rows(cov_data_6000_rct, cov_data_6000_nonrct,
                      cov_data_12000_rct, cov_data_12000_nonrct, 
                      cov_data_18000_rct, cov_data_18000_nonrct)

cov_means <- cov_data %>%
  dplyr::group_by(mechanism, level, type, size, study) %>%
  dplyr::summarise(cov_mean = mean(value, na.rm = TRUE), .groups = "drop")

# Merge data frames 
spv_data <- bind_rows(spv_data_6000_rct, spv_data_6000_nonrct, 
                      spv_data_12000_rct, spv_data_12000_nonrct, 
                      spv_data_18000_rct, spv_data_18000_nonrct)%>%
  rename(spv=value)%>%
  mutate(across(-spv, as.factor))

spv_means <- spv_data %>%
  #mutate(mechanism_level = as.factor(paste(mechanism, level, sep = "_"))) %>%
  dplyr::group_by(mechanism,level, type, size, study) %>%
  dplyr::summarise(mean_spv = mean(spv, na.rm = TRUE),
                   .groups = "drop")
  # dplyr::group_by(mechanism, level, type, size, study) %>%
  # dplyr::summarise(mean_spv = mean(spv, na.rm = TRUE),
  #                  .groups = "drop")

mean_width_data <- bind_rows(mean_width_data_6000_rct, mean_width_data_6000_nonrct,
                             mean_width_data_12000_rct, mean_width_data_12000_nonrct,
                             mean_width_data_18000_rct, mean_width_data_18000_nonrct) %>% 
  rename(level=levels, mean_width=value) %>%
  mutate(across(-mean_width, as.factor))
  #mutate(mechanism_level = as.factor(paste(mechanism, level, sep = "_"))) %>%
  #select(-mechanism, -level)


cov_factor_data <- cov_means %>%
  dplyr::mutate(
    level = as.numeric(level),  # numeric for x-axis
    type = factor(type),         # factor for discrete color
    mechanism = factor(mechanism)) %>% 
  dplyr::mutate(cov_factor = cov_mean - (1 - level)) %>% 
  mutate(across(-c(cov_mean,cov_factor ), as.factor))


complete_data <- list(spv_means, mean_width_data,cov_factor_data) %>%
  reduce(full_join, by = c("mechanism","level", "type", "size", "study"))


model <- lm(
  mean_spv ~
    level * mechanism * study + type +
    size -1,
  data = complete_data
)

summary(model)

# With random noise 
combined <- ggplot(complete_data, 
                     aes(x = cov_mean, y = mean_spv)) +
  geom_line(aes(group = interaction(size, type), 
                color = factor(size), 
                linetype = type),
            alpha = 0.3)+
  geom_point(aes(size = mean_width, 
                 color = factor(size)), 
             alpha = 0.5, 
             show.legend = c(size = FALSE)) +
  scale_size_continuous(range = c(0.5, 2)) +
  facet_grid(study ~ mechanism) +
  # Use labs to rename the aesthetics mapped in your geoms
  labs(
    x = expression("Coverage attained ( 1 - " * alpha * ")"), 
    y = "SPV",
    color = "Size",           
    linetype = "Noise level (r)" )
ggplot2::ggsave(combined, filename=paste0("images/",score_name, "Combined_cov_spv_", type,".pdf"), width = 15, height = 8)

# Random noise 
plot_data <- complete_data %>%
  mutate(color_group = case_when(
    mechanism == "Estimated labels" ~ paste0("type_", type),
    mechanism == "True labels" ~ "True labels",
    mechanism == "ULB" ~ "ULB"
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
      "True labels" = "blue",
      "ULB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "True labels", "ULB"),
    labels = c(paste0("r = ", type_vals), "True labels", "ULB")
  ) +
  
  scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  facet_grid(study ~ size) +
  
  labs(x = expression("Confidence level (" * alpha * ")"),
       y = "SPV")
ggplot2::ggsave(plot_combined_level, filename=paste0("images/",score_name, "Combined_level_", type,".pdf"), width = 15, height = 8)


plot_combined_cov <- ggplot(plot_data,
                            aes(x = cov_mean,
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
      "True labels" = "blue",
      "ULB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "True labels", "ULB"),
    labels = c(paste0("r = ", type_vals), "True labels", "ULB")
  ) +
  
  scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  facet_grid(study ~ size) +
  
  labs(x = expression("Coverage level attained ( 1- " * alpha * ")"),
       y = "SPV")
ggplot2::ggsave(plot_combined_cov, filename=paste0("images/",score_name, "Combined_cov_", type,".pdf"), width = 15, height = 8)


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
      "True labels" = "blue",
      "ULB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "True labels", "ULB"),
    labels = c(paste0("r = ", type_vals), "True labels", "ULB")
  ) +
  
  scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  facet_grid(study ~ size) +
  
  labs(x = expression("Confidence level (" * alpha * ")"),
       y = expression("Coverage attained ( 1- " * alpha * ")"))
ggplot2::ggsave(plot_cov_level, filename=paste0("images/",score_name, "Cov_level_", type,".pdf"), width = 15, height = 8)
 
plot_combined_cov_factor <- ggplot(plot_data,
                              aes(x = cov_factor,
                                  y = mean_spv)) +
  geom_line(data = subset(plot_data, mechanism == "Estimated labels"),
            aes(group = level, color = color_group),
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
      "True labels" = "blue",
      "ULB"  = "green"
    ),
    breaks = c(paste0("type_", type_vals), "True labels", "ULB"),
    labels = c(paste0("r = ", type_vals), "True labels", "ULB")
  ) +
  
  scale_size_continuous(name = "Mean width",
                        range = c(0.1, 2)) +
  facet_grid(study ~ size) +
  labs(x = "Marginal coverage factor ",
       y = "SPV")
ggplot2::ggsave(plot_combined_cov_factor, filename=paste0("images/",score_name, "Combined_cov_factor_", type,".pdf"), width = 8, height = 8)

