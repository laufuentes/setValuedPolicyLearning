
SL.out$rate_cal_labels_unweighted <- SL.out$rate_cal_labels_sl <- SL.out$rate_cal_labels_exp <- 
  SL.out$rate_scores_unweighted_cal <- SL.out$rate_scores_sl_cal<- SL.out$rate_scores_exp_cal<- matrix(0,nrow=nrow(test), ncol=n_rate)
if(synthetic_scenario){
  SL.out$rate_scores_unweighted_true<- SL.out$rate_scores_sl_true <- SL.out$rate_scores_exp_true <- matrix(0,nrow=nrow(test), ncol=n_rate)
}  

for(i in 1:n_rate){
  rate <- random_rate[i]
  mix_factor<- stats::rbinom(nrow(test),1,prob=rate)
  # Apply random rate to all tree kind of label generation mechanisms 
  SL.out$rate_cal_labels_unweighted[,i]<- mix_factor*A_rd +  (1-mix_factor)*unweighted_cal
  SL.out$rate_cal_labels_sl[,i]<- mix_factor*A_rd +  (1-mix_factor)*sl_cal
  SL.out$rate_cal_labels_exp[,i]<- mix_factor*A_rd +  (1-mix_factor)*exp_cal
  
  # Compute associated scores 
  SL.out$rate_scores_unweighted_cal[,i] <- margin_po[cbind(1:nrow(test), SL.out$rate_cal_labels_unweighted[,i])]
  SL.out$rate_scores_sl_cal[,i] <- margin_po[cbind(1:nrow(test), SL.out$rate_cal_labels_sl[,i])]
  SL.out$rate_scores_exp_cal[,i] <- margin_po[cbind(1:nrow(test), SL.out$rate_cal_labels_exp[,i])]
  
  if(synthetic_scenario){
    SL.out$rate_scores_unweighted_true[,i] <- true_marginal_scores[cbind(1:nrow(test), SL.out$rate_cal_labels_unweighted[,i])]
    SL.out$rate_scores_sl_true[,i] <- true_marginal_scores[cbind(1:nrow(test), SL.out$rate_cal_labels_sl[,i])]
    SL.out$rate_scores_exp_true[,i] <- true_marginal_scores[cbind(1:nrow(test), SL.out$rate_cal_labels_exp[,i])] 
  }
}

SL.out$rate_cal_labels_behavioral <- margin_po[cbind(1:nrow(test), test[,treatment_name])]

data_toghether <- dplyr::bind_rows(
  lapply(1:ncol(SL.out$rate_scores_unweighted_cal), function(i) {
    data.frame(
      value = SL.out$rate_scores_unweighted_cal[,i],
      model = "Estimated score",
      mechanism = "Unweighted",
      type = random_rate[i]
    )
  }),
  lapply(1:ncol(SL.out$rate_scores_sl_cal), function(i) {
    data.frame(
      value = SL.out$rate_scores_sl_cal[,i],
      model = "Estimated score",
      mechanism = "SL",
      type = random_rate[i]
    )
  }),
  lapply(1:ncol(SL.out$rate_scores_exp_cal), function(i) {
    data.frame(
      value = SL.out$rate_scores_exp_cal[,i],
      model = "Estimated score",
      mechanism = "Exp",
      type = random_rate[i]
    )
  }))

behavioral_data <-  dplyr::bind_rows(
data.frame(
  value = SL.out$rate_cal_labels_behavioral,
  model = "Estimated score",
  mechanism = "Unweighted",
  type = 0
), 
data.frame(
  value = SL.out$rate_cal_labels_behavioral,
  model = "Estimated score",
  mechanism = "SL",
  type = 0
),
data.frame(
  value = SL.out$rate_cal_labels_behavioral,
  model = "Estimated score",
  mechanism = "Exp",
  type = 0
))

if(synthetic_scenario){
  data_toghether<- dplyr::bind_rows(data_toghether,
            lapply(1:ncol(SL.out$rate_scores_unweighted_true), function(i) {
              data.frame(
                value = SL.out$rate_scores_unweighted_true[,i],
                model = "True score",
                mechanism = "Unweighted",
                type = random_rate[i]
              )
            }),
            lapply(1:ncol(SL.out$rate_scores_sl_true), function(i) {
              data.frame(
                value = SL.out$rate_scores_sl_true[,i],
                model = "True score",
                mechanism = "SL",
                type =random_rate[i]
              )
            }),
            lapply(1:ncol(SL.out$rate_scores_exp_true), function(i) {
              data.frame(
                value = SL.out$rate_scores_exp_true[,i],
                model = "True score",
                mechanism = "Exp",
                type = random_rate[i]
              )
            }))

  data_true_all <- dplyr::bind_rows(
    data.frame(
      value = as.vector(SL.out$true_score),
      model = "Estimated score",
      mechanism = c("Unweighted", "SL", "Exp")[1]
    ),
    data.frame(
      value = as.vector(SL.out$true_score),
      model = "Estimated score",
      mechanism = c("Unweighted", "SL", "Exp")[2]
    ),
    data.frame(
      value = as.vector(SL.out$true_score),
      model = "Estimated score",
      mechanism = c("Unweighted", "SL", "Exp")[3]
    ))
  
  data_true_all <- dplyr::bind_rows(data_true_all, 
                             data.frame(
                               value = as.vector(SL.out$true_score_true),
                               model = "True score",
                               mechanism = "Unweighted"
                             ),
                             data.frame(
                               value = as.vector(SL.out$true_score_true),
                               model = "True score",
                               mechanism = "SL"
                             ),
                             data.frame(
                               value = as.vector(SL.out$true_score_true),
                               model = "True score",
                               mechanism = "Exp"
                             ))
}