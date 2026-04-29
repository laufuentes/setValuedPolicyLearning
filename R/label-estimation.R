#' Q-model builder for Policy Learning
#'
#' A wrapper function to build Q-models required for training policies using the 
#' polle package. It simplifies the creation of formula objects based on 
#' treatment restrictions.
#'
#' @param q_func String indicating the type of function for Q-learning. 
#'  c("q_glm", "q_rf", "q_xgboost", and "q_sl").
#' @param covariates Character vector of covariate names. Defaults to `c("x1", "x2")`.
#' @param sl_library String vector of libraries for learning policies 
#'   using SuperLearner. Required if q_func = "q_sl".
#' @param action_name String indicating the treatment/action variable 
#'   name in the data frame. Defaults to "A".
#'   
#' @return An object of class q_model as defined by the polle package.
#' @export
#' 
#' @examples
#' \dontrun{
#' make_q_model(q_func = "q_glm", action_name = "treatment")
#' }
make_q_model <- function(q_func, covariates=c("x1","x2"), sl_library = NULL, action_name = "A") {
  formula_obj <- stats::reformulate(covariates)
  switch(q_func,
         "q_glm" = polle::q_glm(formula = formula_obj, family = gaussian()),
         "q_rf" = polle::q_rf(formula = formula_obj),
         "q_xgboost" = polle::q_xgboost(formula = formula_obj),
         "q_sl" = {
           if (is.null(sl_library)) stop("'q_sl' requires sl_library.")
           polle::q_sl(formula = formula_obj, SL.library = sl_library)},
         stop("Unknown q_func '", q_func, "'")
  )
}

#' Add a Q-learner to a list of learners
#'
#' Appends a learner configuration to a list and automatically generates the 
#' associated `q_model` and `qv_model` using `polle` specifications.
#'
#' @param learners List to which the learner configuration will be added.
#' @param name String indicating the unique name of the expert.
#' @param type String indicating the type of learner. 
#'   ("ql", "drql", "ptl", "policytree").
#' @param action_name String indicating the treatment variable name. Defaults to "A".
#' @param covariates Character vector of covariate names. Defaults to `c("x1", "x2")`.
#' @param q_func String indicating the type of function for Q-learning. 
#'   ("q_glm", "q_rf", "q_xgboost", "q_sl"). Defaults to NULL.
#' @param sl_library Vector of libraries for SuperLearner. 
#'   Required if `q_func` is "q_sl". Defaults to NULL. 
#' @param depth Numeric value for the depth of growing trees. Defaults to 2. 
#' @param hybrid Logical indicating whether to use a hybrid approach. Defaults to FALSE.
#'
#' @return The updated `learners` list with the new learner and its characteristics.
#' @export
#'
#' @examples
#' \dontrun{
#' learners <- list()
#' learners <- add_qlearner(
#'   learners, 
#'   name = "ql_lm", 
#'   type = "ql", 
#'   q_func = "q_glm", 
#'   action_name = "A"
#' )
#' }
add_qlearner <- function(learners, name, 
                         type=c("ql", "drql", "ptl", "policytree"), 
                         action_name = "A", covariates=c("x1","x2"), 
                         q_func = NULL, sl_library = NULL,
                         depth = 2, hybrid = FALSE) {
  
  if (type != "policytree") {
    valid_q_funcs <- c("q_glm", "q_rf", "q_xgboost", "q_sl")
    if (is.null(q_func) || !q_func %in% valid_q_funcs)
      stop("Unknown q_func '", q_func, "'. Choose from: ", paste(valid_q_funcs, collapse = ", "))
    qv_model <- make_q_model(q_func,
                             covariates   = covariates,
                             sl_library   = sl_library,
                             action_name  = action_name)
  } else {
    qv_model <- NULL
  }
  
  learners[[name]] <- list(name     = name,
                           type     = type,
                           q_func   = q_func,
                           qv_model = qv_model, 
                           depth    = depth,
                           hybrid   = hybrid)
  message("+ '", name, "' added.")
  learners
}

#' Builder of G-model for Policy Learning
#'
#' A wrapper function to build G-models required for training policies using the 
#' `polle` package.
#'
#' @param g_func String indicating the type of function for Q-learning. 
#'  ("g_glm", "g_glmnet", "g_rf", "g_xgboost", and "q_sl").
#' @param sl_library String vector of libraries for learning policies 
#'   using SuperLearner. Required if g_func = "g_sl". Defaults to NULL.
#' @param num.trees Number of trees for growing trees. 
#' Required if g_func = "g_rf".  Defaults to NULL.
#'   
#' @return An object of class q_model as defined by the polle package.
#' @export
#' 
#' @examples
#' \dontrun{
#' make_g_model(g_func = "g_glm")
#' }
make_g_model <- function(g_func, sl_library = NULL, num.trees = NULL) {
  switch(g_func,
         "g_glm"     = polle::g_glm(family = "binomial"),
         "g_glmnet"  = polle::g_glmnet(family = "binomial"),
         "g_rf"      = {
           if (is.null(num.trees)) stop("'g_rf' requires num.trees.")
           polle::g_rf(num.trees = num.trees)
         },
         "g_xgboost" = polle::g_xgboost(objective = "binary:logistic"),
         "g_sl"      = {
           if (is.null(sl_library)) stop("'g_sl' requires sl_library.")
           polle::g_sl(SL.library = sl_library)
         },
         stop("Unknown g_func '", g_func, "'. Choose from: g_glm, 
              g_glmnet, g_rf, g_xgboost, g_sl")
  )
}

#' G-model creator
#'
#' Generates the associated `g_model` using `polle` specifications.
#'
#' @param m Numeric value for the number of treatment levels. 
#' @param g_func String indicating the type of function for G-learning. 
#'  ("g_glm", "g_glmnet", "g_rf", "g_xgboost", and "q_sl").
#' @param sl_library String vector of libraries for learning policies 
#'   using SuperLearner. Required if g_func = "g_sl".
#' @param num.trees Number of trees for growing trees.  Defaults to NULL.
#' Required if g_func = "g_rf". Defaults to NULL.
#'
#' @return Creates a learner 
#' @export
#'
#' @examples
#' \dontrun{
#' g_model <- glearner(
#'   m = 5, 
#'   g_func = "g_glm")
#' }
glearner <- function(m, g_func = NULL, sl_library = NULL, num.trees = NULL) {
  valid_g_funcs <- if (m > 2) {c("g_rf")}
  else {c("g_rf", "g_xgboost", "g_glm", "g_glmnet", "g_sl")}
  
  if (is.null(g_func) || !g_func %in% valid_g_funcs)
    stop("Unknown g_func '", g_func, "' for m=", m, ". Choose from: ",
         paste(valid_g_funcs, collapse = ", "))
  make_g_model(g_func = g_func, sl_library = sl_library, num.trees=num.trees)
}


#' Fit and Predict Expert Policies
#'
#' Fits multiple policy learners (Q-learning, DR Q-learning, PTL, or PolicyTree) 
#' using the polle and policytree packages, and generates predictions 
#' for test and optional new data.
#'
#' @param train Data frame used for training the models.
#' @param test Data frame used for prediction (calibration/test set).
#' @param new Optional data frame for prediction on new observations. Defaults to NULL.
#' @param covariates Character vector of covariate names. Defaults to `c("x1", "x2")`.
#' @param treatment_name String indicating the treatment variable. Defaults to "A".
#' @param outcome_name String indicating the outcome variable. Defaults to "Y".
#' @param qlearners_list A list of learner configurations created by add_qlearner.
#' @param g_model The pre-specified g-model (propensity score model) for polle.
#'
#' @return A list containing:
#' \itemize{
#'   \item `learner_names`: Names of the processed learners.
#'   \item `expert_policies`: Matrix of predictions for the `test` data.
#'   \item `expert_policies_new`: Matrix of predictions for `new` data (if provided).
#' }
#' @export
expert_fit_predict <- function(train, test, new = NULL,
                               covariates = c("x1", "x2"),
                               treatment_name = "A", outcome_name = "Y",
                               qlearners_list, g_model) {
  
  train <- data.table::as.data.table(train)
  test  <- data.table::as.data.table(test)
  if (!is.null(new)) new <- data.table::as.data.table(new)
  
  train[[treatment_name]] <- as.factor(train[[treatment_name]])
  test[[treatment_name]]  <- as.factor(test[[treatment_name]])
  if (!is.null(new)) new[[treatment_name]] <- as.factor(new[[treatment_name]])
  
  pd_train <- polle::policy_data(train, action = treatment_name,
                                 covariates = covariates, utility = outcome_name)
  pd_test  <- polle::policy_data(test,  action = treatment_name, 
                                 covariates = covariates, utility = outcome_name)
  pd_new   <- if (!is.null(new))
    polle::policy_data(new, action = treatment_name, 
                       covariates = covariates, utility = outcome_name)
  else NULL
  
  expert_policies     <- list()
  expert_policies_new <- list()
  
  for (learner in qlearners_list) {
    name <- learner$name
    cat("Processing:", name, "\n")
    
    if (learner$type == "policytree") {
      # Train policytree model
      forest <- grf::multi_arm_causal_forest(
        X = as.matrix(train[, ..covariates]),
        Y = train[[outcome_name]],
        W = train[[treatment_name]]
      )
      dr_scores <- policytree::double_robust_scores(forest)
      
      tree <- if (isTRUE(learner$hybrid)) {
        # hybrid tree
        policytree::hybrid_policy_tree(as.matrix(train[, ..covariates]),
                                       Gamma=dr_scores, depth = learner$depth)
      } else {
        # regular tree
        policytree::policy_tree(as.matrix(train[, ..covariates]), 
                                Gamma=dr_scores, depth = learner$depth)
      }
      
      # Make policytree prediction 
      expert_policies[[name]] <- as.numeric( # test 
        predict(tree, as.matrix(test[, ..covariates]))
        )
      if (!is.null(new))  # new data  
        expert_policies_new[[name]] <- as.numeric(
          predict(tree, as.matrix(new[, ..covariates]))
          )
    
    # polle learners 
    } else if (learner$type == "ql") {
      # policy learner (classic Q-learning)
      l_obj <- polle::policy_learn(type = "ql", cross_fit_g_models = FALSE)
      # policy optimization
      po <- l_obj(policy_data = pd_train, 
                  q_models = learner$qv_model, g_models = g_model)
      
    } else if (learner$type == "drql") {
      # policy learner (DR Q-learning)
      l_obj <- polle::policy_learn(
        type    = "drql",
        control = polle::control_drql(qv_models = learner$qv_model), 
        cross_fit_g_models = FALSE
      )
      # policy optimization
      po <- l_obj(
        policy_data = pd_train, 
        q_models = learner$qv_model,
        g_models = g_model)
    } else if (learner$type == "ptl") {
      # policy learner (tree-based)
      l_obj <- polle::policy_learn(
        type    = "ptl",
        control = polle::control_ptl(depth = learner$depth, 
                                     hybrid = learner$hybrid),
        cross_fit_g_models = FALSE)
      # policy optimization
      po <- l_obj(policy_data = pd_train, 
                  q_models = list(learner$qv_model), g_models = g_model)
    }
    
    # Make polle predictions 
    if (learner$type != "policytree") {
      # test
      expert_policies[[name]] <- as.numeric(polle::get_policy(po)(pd_test)$d)
      if (!is.null(new)) #new 
        expert_policies_new[[name]] <- as.numeric(polle::get_policy(po)(pd_new)$d)
    }
  }
  
  # Output list
  list(
    learner_names = names(expert_policies),
    expert_policies = do.call(cbind, expert_policies),
    expert_policies_new = if (length(expert_policies_new) > 0)
      do.call(cbind, expert_policies_new) else NULL
  )
}




#' Weighted Probability Distribution from Experts
#'
#' Computes the probability scores associated with covariates across a set of 
#' experts given a vector of expert weights.
#'
#'
#' @param df_pred Data frame containing the covariates for prediction.
#' @param fitted_experts Matrix where rows are observations and columns are 
#'   individual expert policy predictions.
#' @param weights Numeric vector of weights assigned to each expert. 
#'   Length must match the number of columns in `fitted_experts`.
#' @param levels Vector of possible treatment/action levels. Defaults to `1:5`.
#'
#' @return A numeric `matrix` of dimensions `nrow(df_pred)` by `length(levels)` 
#'   where each cell `(i, a)` represents the weighted probability of choosing 
#'   action `a` for observation `i`.
#' @export
weighted_probs_experts <- function(df_pred, fitted_experts, weights, levels=1:5){
  l <- length(levels)
  n <- nrow(df_pred)
  one_hot_experts <- outer(fitted_experts, 1:l, "==") + 0
  
  scores <- matrix(rep(0,n*l),nrow=n,ncol=l)
  for(a in levels) {
      pi_ax <- one_hot_experts[,,a]
      scores[,a] <- as.vector(weights%*%t(pi_ax))
    }
  return(scores)
}
