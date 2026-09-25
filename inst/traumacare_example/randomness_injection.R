# ── 3) Noisy calibration  ───────────────────────────────────────────────────
A_rd <- apply(data.frame(1:nrow(calibration)),1,function(i)sample(as.numeric(levels_A),size=1))

# 3.1) Generate perturbed labels  ────────────────────────────────────────────
rate_cal_labels_unweighted <- rate_cal_labels_single <- rate_scores_unweighted_cal <- rate_scores_single_cal<- matrix(0,nrow=nrow(calibration), 
                                                                                                                      ncol=n_rate)

## Combine noisy labels with random  using different randomness levels (r)
for(i in 1:n_rate){
  rate <- random_rate[i] # randomness level r
  mix_factor<- stats::rbinom(nrow(calibration),1,prob=rate) # R ~ Ber(r)
  # Combine noisy labels with random
  rate_cal_labels_unweighted[,i]<- mix_factor*A_rd +  (1-mix_factor)*SL.out$unweighted_cal
  rate_cal_labels_single[,i]<- mix_factor*A_rd +  (1-mix_factor)*SL.out$single_policy_cal
  # Compute associated scores
  rate_scores_unweighted_cal[,i] <- margin_po[cbind(1:nrow(calibration), 
                                                    rate_cal_labels_unweighted[,i]+1)]
  rate_scores_single_cal[,i] <- margin_po[cbind(1:nrow(calibration), 
                                                rate_cal_labels_single[,i]+1)]
}
