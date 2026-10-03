root.path <- "~/Documents/PhD/Project 2 - Conformal Policy Sets /setValuedPolicyLearning/"
setwd(root.path)
seed <- 2026
set.seed(seed)

subdirs <- c("images", "intermediate")
for (subdir in subdirs) {
  subdir_path <- file.path(root.path, "inst", 
                           "traumacare_example", 
                           subdir)
  
  if (!dir.exists(subdir_path)) {
    message(sprintf("Creating subdirectory: %s", subdir))
    dir.create(subdir_path, recursive = TRUE)
  }
}
# ── Required packages  ────────────────────────────────────────────────────────
source("inst/libraries.R")

# ── Load data ─────────────────────────────────────────────────────────────────
train <- read.csv("inst/traumacare_example/cohort_train.csv")
test <- read.csv("inst/traumacare_example/cohort_test.csv")

# ── Change NA to "Non_applicable" ────────────────────────────────────────────
train$CAUSE_DECES[which(train$CAUSE_DECES %>% is.na())] <- "Non_applicable"
test$CAUSE_DECES[which(test$CAUSE_DECES %>% is.na())] <- "Non_applicable"

# ── Replace ND, NR and IMP ────────────────────────────────────────────────────
train[(train=="ND" | train=="NR" | train =="IMP")] <- NA
test[(test=="ND" | test=="NR" | test =="IMP")] <- NA

# People that died and we dont know if it they had NEUROCHIR  ──────────────────
remove_ids <- which((train$DECES=="Oui") &
                      (as.numeric(train$DELARHDECH)<24) &
                      (train$CAUSE_DECES!="Trauma cranien") & 
                      (train$NEUROCHIR!=1))
train <- train[if(length(remove_ids)>0){-remove_ids}else{1:nrow(train)},]

# ── Transform numeric variables as numeric  ───────────────────────────────────
vars_to_num <- c(#"AGE", 
                 "FC_ARRIVEE_MEDIC", 
                 "GLASGOW_INITIAL", 
                 "GLASGOW_MOTEUR_INIT", 
                 "HEMOCUE_INITIAL", 
                 #"PAD_ARRIVEE_MEDIC",
                 "PAS_ARRIVEE_MEDIC", 
                 "SHOCK_INDEX", 
                 "SHOCK_INDEX_DIASTOLIQUE", 
                 "SHOCK_INDEX_INVERSE", 
                 "DELARHDECH")

train <- train %>% mutate(across(all_of(vars_to_num), as.numeric))
test  <- test  %>% mutate(across(all_of(vars_to_num), as.numeric))

# ── Define outputs  ──────────────────────────────────────────────────────────
# Died after 24h and before 28 days 
train$Y_after_24h_before_28d <- ifelse((train$DECES=="Oui") & 
                                         as.numeric(train$DELARHDECH)>=24 & 
                                         as.numeric(train$DELARHDECH)<=672, 0, 1) 
test$Y_after_24h_before_28d <- ifelse((test$DECES=="Oui") & 
                                        as.numeric(test$DELARHDECH)>=24 & 
                                        as.numeric(test$DELARHDECH)<=672, 0, 1) 

outcome_name <- "Y_after_24h_before_28d"

# Delete NAs in ouputs for train set (we don't know if they died)
to_delete_train <- which(train$Y_after_24h_before_28d %>% is.na()) # TODO : test that train$Y_less_28d has same NAs
train<- train[-to_delete_train,] 

# Transform categorical variables in factors
train <- train %>% mutate(across(where(is.character), 
                                 as.factor))
test <- test %>% mutate(across(where(is.character), 
                               as.factor))

# Change treatment (NEUROCHIR) to factor
train$HEMO_SHOCK <- as.factor(train$HEMO_SHOCK)
test$HEMO_SHOCK <- as.factor(test$HEMO_SHOCK)
#train$NEUROCHIR <- as.factor(train$NEUROCHIR)
#test$NEUROCHIR <- as.factor(test$NEUROCHIR)
# Note that 0 means died and 1 means survived 

# ── Data imputation  ──────────────────────────────────────────────────────────
# Variables post-NEUROCHIR (not to be imputed)
remove_variables <- c("HEMO_SHOCK", "DECES", "CAUSE_DECES", "DELARHDECH", #"NEUROCHIR"
                      "SUBJECT_REF","Y_after_24h_before_28d") #"Y_less_24h", "Y_less_28d",

treatment_name <- "HEMO_SHOCK" #"NEUROCHIR" # treatment indicator string

# Imputation
Subject_ref_train_test <- rbind(train,test)[,1]
imp.train <- mice(rbind(train,test)[,-1], ignore=c(rep(FALSE, nrow(train)),
                                                   rep(TRUE, nrow(test)))) # train
imp_df <- complete(imp.train) # complete data imputation

# train data complete
train_imp <- imp_df[1:nrow(train),] %>% 
  select(-c("DECES", "CAUSE_DECES", "DELARHDECH")) 

train_imp <- train_imp %>% 
  mutate("SUBJECT_REF"=Subject_ref_train_test[1:nrow(train_imp)])

# test data complete
test_imp <- imp_df[(nrow(train)+1):nrow(imp_df),] 

write.csv(mean(test_imp[,outcome_name]),
          file=paste0("inst/traumacare_example/intermediate/Clinicians", "_", outcome_name, ".csv"))

test_imp <- test_imp %>% 
  mutate("SUBJECT_REF"=Subject_ref_train_test[
    (nrow(train_imp)+1):length(Subject_ref_train_test)]) %>% 
  select(-all_of(remove_variables[remove_variables!="SUBJECT_REF"])) # all_of(remove_variables)

# re-convert categories to factors 
train_imp <- train_imp %>% 
  mutate(across(where(is.character), as.factor))

test_imp <- test_imp %>% 
  mutate(across(where(is.character), as.factor))

covariate_name <- setdiff(colnames(train_imp), 
                          c(remove_variables)) # covariates indicator string

true_outputs_test <- imp_df[(nrow(train_imp)+1):nrow(imp_df), c(outcome_name,treatment_name)] %>% 
  mutate("SUBJECT_REF" = Subject_ref_train_test[
    (nrow(train_imp)+1):length(Subject_ref_train_test)])

# ── One hot encoding ──────────────────────────────────────────────────────────
# variables to one-hot encode
categ_var <- c("AMPUTATION",
              "ANOMAL_PUPIL_PREHOSP", 
              "FRACAS_BASSIN", 
              "HEMORRAGIE_EXTERNE", 
              #"IOT_PREHOSP", 
              "ISCHEMIE_MEMBRE", 
              "MECANISME_CAUSE", 
              #"OSMOTHERAPIE", 
              #"PERTCONINT",
              "REA_CATECHO")  


dummies <- dummyVars(as.formula(paste("~", 
                                      paste(categ_var, 
                                            collapse = " + "))), 
                     data= train_imp)

# One-hot encoded training set
train_one_hot <- predict(dummies, 
                         newdata = train_imp)

# One-hot encoded test set
test_one_hot <- predict(dummies, 
                        newdata = test_imp)

# Adapt string format for some categories
colnames(train_one_hot) <- colnames(train_one_hot) %>%
  gsub(pattern = " ", replacement = "_") %>%
  gsub(pattern = "-", replacement = "_") %>%
  gsub(pattern = "/", replacement = "ou") %>%
  gsub(pattern = ",", replacement = "") %>%
  gsub(pattern = "'", replacement = "") %>% 
  gsub(pattern = "\\.", replacement = "_") %>%
  gsub(pattern = "\\(", replacement = "") %>%
  gsub(pattern = "\\)", replacement = "") # training one-hot encoded version


colnames(test_one_hot) <-  colnames(test_one_hot) %>%
  gsub(pattern = " ", replacement = "_") %>%
  gsub(pattern = "-", replacement = "_") %>%
  gsub(pattern = "/", replacement = "ou") %>%
  gsub(pattern = ",", replacement = "") %>%
  gsub(pattern = "'", replacement = "") %>% 
  gsub(pattern = "\\.", replacement = "_")%>% 
  gsub(pattern = "\\(", replacement = "") %>%
  gsub(pattern = "\\)", replacement = "") # testing one-hot encoded version

# Select only .Oui variables and exclude one category for MECANISME_CAUSE 
# (for identifiability)
categ_var_selected <- c(
  "AMPUTATION_Oui", 
  "ANOMAL_PUPIL_PREHOSP_Oui", 
  "FRACAS_BASSIN_Oui", 
  "HEMORRAGIE_EXTERNE_Oui", 
  #"IOT_PREHOSP_Oui", 
  "ISCHEMIE_MEMBRE_Oui",
  "MECANISME_CAUSE_Accident_en_montagne_ou_activité_en_plein_air",
  "MECANISME_CAUSE_Arme_à_feu",
  "MECANISME_CAUSE_Arme_blanche",
  "MECANISME_CAUSE_Autre_traumatisme_fermé",
  "MECANISME_CAUSE_Autre_traumatisme_pénétrant",
  "MECANISME_CAUSE_AVP_voiture_camion_bus",
  "MECANISME_CAUSE_AVP_autre_avion_bateau_train_autres",
  "MECANISME_CAUSE_AVP_bicyclette",
  "MECANISME_CAUSE_AVP_deux_roues_motorisé", 
  "MECANISME_CAUSE_AVP_engin_de_déplacement_personnel_motorisé_EDPM",
  "MECANISME_CAUSE_AVP_piéton",
  "MECANISME_CAUSE_Chute_dune_hauteur",
  "MECANISME_CAUSE_Chute_de_sa_hauteur",
  "MECANISME_CAUSE_Traumatisme_par_objet_contondant_non_pénétrant", 
  #"OSMOTHERAPIE_Oui",
  #"PERTCONINT_Oui", 
  "REA_CATECHO_Oui")

# NOTES: 
# for bool we kept ".Oui"
# For MECANISME_CAUSE: we excluded ".Inconnu" 

# Create the one-hot encoded version of training data
train_imp_onehot <- cbind(train_imp %>% select(-all_of(categ_var)),
                          train_one_hot[,categ_var_selected]) 
train_imp_onehot <- train_imp_onehot %>% 
  mutate(across(all_of(categ_var_selected), as.factor))

write.csv(train_imp_onehot, row.names = FALSE,
          file= file.path("inst/traumacare_example/intermediate/train_imp_one_hot.csv"))

# Create the one-hot encoded version of testing data
test_imp_onehot <- cbind(test_imp %>% select(-all_of(categ_var)),
                          test_one_hot[,categ_var_selected])
test_imp_onehot <- test_imp_onehot %>% 
  mutate(across(all_of(categ_var_selected), as.factor))

write.csv(test_imp_onehot, row.names = FALSE,
          file= file.path("inst/traumacare_example/intermediate/test_imp_one_hot.csv"))

write.csv(merge(test_imp_onehot, true_outputs_test, 
                by=intersect(names(test_imp_onehot), 
                             names(true_outputs_test))), 
          file= file.path("inst", "traumacare_example",
                          "intermediate",
                          "test-with-unavailable_at_ambulance.csv"), 
          row.names = FALSE)


covariate_name_one_hot <- setdiff(colnames(train_imp_onehot), 
                                  c(remove_variables)) # covariates indicator string (for one-hot encoded data)

config_conformal_preprocessing <- list("covariates_name" = covariate_name_one_hot, 
                                       "outcome_name" = outcome_name,
                                       "treatment_name" = treatment_name)

saveRDS(config_conformal_preprocessing,
           file = "inst/traumacare_example/intermediate/preprocessing.rds")
