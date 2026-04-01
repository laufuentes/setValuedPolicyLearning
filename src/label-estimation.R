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
    q_model  <- make_q_model(q_func,
                             covariates   = covariates,
                             sl_library   = sl_library,
                             action_name  = action_name)
    qv_model <- make_q_model(q_func,
                             covariates   = covariates,
                             sl_library   = sl_library,
                             action_name  = action_name)
  } else {
    q_model  <- NULL
    qv_model <- NULL
  }
  
  learners[[name]] <- list(name     = name,
                           type     = type,
                           q_func   = q_func,
                           q_model  = q_model,
                           qv_model = qv_model, 
                           depth    = depth,
                           hybrid   = hybrid)
  message("+ '", name, "' added.")
  learners
}

#' Builder of G-model for Policy Learning
#'
#' A wrapper function to build G-models required for training policies using the 
#' polle package.
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
                  q_models = list(learner$q_model), g_models = g_model)
      
    } else if (learner$type == "drql") {
      # policy learner (DR Q-learning)
      l_obj <- polle::policy_learn(
        type    = "drql",
        control = polle::control_drql(qv_models = list(learner$qv_model)), 
        cross_fit_g_models = FALSE)
      # policy optimization
      po <- l_obj(policy_data = pd_train, 
                  q_models = list(learner$qv_model), g_models = g_model)
      
    } else if (learner$type == "ptl") {
      # policy learner (tree-based)
      l_obj <- polle::policy_learn(
        type    = "ptl",
        control = polle::control_ptl(depth = learner$depth, 
                                     hybrid = learner$hybrid),
        cross_fit_g_models = FALSE)
      # policy optimization
      po <- l_obj(policy_data = pd_train, 
                  q_models = list(learner$q_model), g_models = g_model)
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

#' Cross-Validation Risk Function (Internal)
#'
#' Evaluates the risk of policy combinations using a specific fold of cross-validation.
#' **Note:** This function is adapted from package `SL.ODTR` from 
#' Montoya, L. M., van der Laan, M. J., Luedtke, A. R., Skeem, J. L., Coyle, J. R., 
#' & Petersen, M. L. (2023). The optimal dynamic treatment rule superlearner: 
#' considerations, performance, and application to criminal justice interventions.
#' 
#' Relies on global variables `X`, `A`, `Y`, `folds`, `simplex.grid`, and `m` being defined.
#'
#' @param i Integer indicating the index of the fold to use from `folds`.
#'
#' @return A list containing `risk.combos.test` (negative TMLE psi) and 
#'   `risk.var.combos.test` (variance of the influence curve).
#' @export
CV.risk_fun = function(i){
  test_ind = folds[[i]]
  train_ind <- setdiff(1:n, test_ind)
  g.reg.train <- randomForest::randomForest(x = X[train_ind,, drop=F], y = A[train_ind])
  
  gAW.pred = predict(g.reg.train,newdata = data.frame(X[test_ind,]), type = "prob")
  gAW_bounded <- pmax(gAW.pred, 0.01)
  
  Xtr <- data.frame(A, X)[train_ind,]
  cols <- colnames(Xtr)
  QAW.reg.train = SuperLearner(Y = Y[train_ind], 
                               X = Xtr, 
                               SL.library=c("SL.mean","SL.glm", "SL.xgboost",
                                            "SL.ranger", "SL.randomForest", "SL.ksvm"), 
                               family = family)
  experts <- expert_fit_predict(data.frame(X,A,Y)[train_ind,],data.frame(X,A,Y)[test_ind,])
  experts_name <- experts[[1]]
  
  expert_policies <- experts[[2]]
  experts_one_hot <- (outer(expert_policies, 1:m, "==") + 0) 
  results_df <- array(0, dim = c(length(test_ind), grid.size, m))
  n_combos <- nrow(simplex.grid)
  for(l in 1:m) {
    results_df[,,l] <- experts_one_hot[, ,l]%*%t(simplex.grid)
  }
  dopt.combos.test <- apply(results_df, c(1, 2), which.max)
  
  Qdopt.combos.test = sapply(1:nrow(simplex.grid), function(x){
    newdata_temp <- data.frame(X[test_ind,,drop=F], A = factor(dopt.combos.test[,x], levels=1:m))
    newdata_temp <- newdata_temp[, cols]
    predict(QAW.reg.train, newdata = newdata_temp, type = "response")$pred})
  tmle.obj.test = lapply(1:nrow(simplex.grid), function(x) tmle.d.fun(A = A[test_ind], Y = matrix(Y[test_ind]), d = dopt.combos.test[,x], Qd = Qdopt.combos.test[,x], gAW = gAW_bounded[cbind(1:length(test_ind),dopt.combos.test[,x])], ab = ab))
  risk.combos.test = -unlist(lapply(tmle.obj.test, function(x) x$psi))
  risk.var.combos.test = unlist(lapply(tmle.obj.test, function(x) var(x$IC)))
  toreturn = list(risk.combos.test = risk.combos.test, risk.var.combos.test = risk.var.combos.test)
  return(toreturn)
}

#' Cross-Validation Risk Function on one fold (Internal)
#'
#' Evaluates the risk of policy combinations using a specific fold of cross-validation.
#' **Note:** This function is adapted from package `SL.ODTR` and relies on 
#' global variables `X`, `A`, `Y`, `folds`, `simplex.grid`, and `m` being defined.
#'
#' @param i Integer indicating the index of the fold to use from `folds`.
#'
#' @return A list containing `risk.combos.test` (negative TMLE psi) and 
#'   `risk.var.combos.test` (variance of the influence curve).
#' @export
CV.risk_fun_one_fold = function(train_ind, test_ind){
  g.reg.train <- randomForest::randomForest(x = X[train_ind,, drop=F], y = A[train_ind])
  
  gAW.pred = predict(g.reg.train,newdata = data.frame(X[test_ind,]), type = "prob")
  gAW_bounded <- pmax(gAW.pred, 0.01)
  
  Xtr <- data.frame(A, X)[train_ind,]
  cols <- colnames(Xtr)
  QAW.reg.train = SuperLearner(Y = Y[train_ind], 
                               X = Xtr, 
                               SL.library=SL.library.nuisance, 
                               family = family)
  
  train_df <- data.frame(X,A,Y)[train_ind,] 
  colnames(train_df) <- c(colnames(X), treatment_name, outcome_name)
  
  test_df <- data.frame(X,A,Y)[test_ind,]
  colnames(test_df) <- c(colnames(X), treatment_name, outcome_name)
  
  experts <- expert_fit_predict(train_df,
                                test_df, 
                                covariates = covariates_name, 
                                treatment_name = treatment_name,
                                outcome_name = outcome_name, 
                                qlearners_list = q_learners, 
                                g_model= g_learner)
  experts_name <- experts[[1]]
  
  expert_policies <- experts[[2]]
  experts_one_hot <- (outer(expert_policies, 1:m, "==") + 0) 
  results_df <- array(0, dim = c(length(test_ind), grid.size, m))
  n_combos <- nrow(simplex.grid)
  for(l in 1:m) {
    results_df[,,l] <- experts_one_hot[, ,l]%*%t(simplex.grid)
  }
  dopt.combos.test <- apply(results_df, c(1, 2), which.max)
  
  Qdopt.combos.test = sapply(1:nrow(simplex.grid), function(x){
    newdata_temp <- data.frame(X[test_ind,,drop=F], A = factor(dopt.combos.test[,x], levels=1:m))
    newdata_temp <- newdata_temp[, cols]
    predict(QAW.reg.train, newdata = newdata_temp, type = "response")$pred})
  tmle.obj.test = lapply(1:nrow(simplex.grid), function(x) tmle.d.fun(A = A[test_ind], Y = matrix(Y[test_ind]), d = dopt.combos.test[,x], Qd = Qdopt.combos.test[,x], gAW = gAW_bounded[cbind(1:length(test_ind),dopt.combos.test[,x])], ab = ab))
  risk.combos.test = -unlist(lapply(tmle.obj.test, function(x) x$psi))
  risk.var.combos.test = unlist(lapply(tmle.obj.test, function(x) var(x$IC)))
  toreturn = list(risk.combos.test = risk.combos.test, risk.var.combos.test = risk.var.combos.test)
  return(toreturn)
}

#' Exponential Weights for Expert Aggregation
#'
#' Computes the updated probability distribution (q2) over a set of experts 
#' using an exponential weighting scheme based on estimated losses.
#'
#' @param librarydoptFactorPredict Matrix of expert policy predictions.
#' @param train Data frame used to train the loss estimator.
#' @param test Data frame used to evaluate experts and calculate potential outcomes.
#' @param experts_library List of experts being evaluated.
#' @param m Numeric value for the number of treatment levels. 
#' @param covariates_name Character vector of covariate names. Defaults to `c("x1", "x2")`.
#' @param treatment_name String indicating the treatment variable. Defaults to "A".
#' @param outcome_name String indicating the outcome variable. Defaults to "Y".
#' @param eta Numeric learning rate (step size) for the exponential update. Defaults to 5.
#' @param levels_A Vector of possible treatment levels. Defaults to `1:m`.
#' @param SL.library String vector of libraries for estimating potential outcomes 
#'   using SuperLearner.
#'
#' @return A numeric vector representing the updated probability distribution (`q2`) 
#'   over the experts.
#' @export
exponential_weights <- function(librarydoptFactorPredict, train, test, 
                                experts_library, m=5, 
                                covariates_name=c("x1","x2"), treatment_name="A", 
                                outcome_name="Y", eta=5, levels_A = 1:m,
                                SL.library){
  # q1 is the uniform distribution over {1,...,J} experts
  # Step 1: Get advice from experts (trained on fold 1)
  experts_test_01 <- outer(librarydoptFactorPredict, 1:m, "==") + 0
  
  # Step 2: Draw an arm from the probability distribution 
  # p_{i}= E_{j\sim q_{1}}[\xi^{j}_{i}] with i \in \{1,...,m\}
  # 2.1 - Define the probability distribution 
  p <- apply(experts_test_01, c(1, 3), function(x) sum(x * rep(1/length(experts_library))))
  p1<- pmax(p, 1e-5)
  # 2.2 - Draw arm 
  # Step 3: Compute estimated loss (trained on fold 2)
  ## The chosen estimated loss is l(x,a) = -Q(a,x)
  mod_sl <- SuperLearner(Y=train[,outcome_name], X = train[,c(covariates_name,treatment_name)], 
                         family = gaussian(),
                         SL.library = SL.library)
  potential_outcomes <- lapply(levels_A, function(val) {
    new_data <- test[, c(covariates_name, treatment_name)]
    new_data[, treatment_name] <- factor(val, levels=levels_A)
    predict.SuperLearner(mod_sl, newdata = new_data, onlySL = TRUE)$pred})
  
  loss_potential_outcomes <- -do.call(cbind, potential_outcomes)
  
  # Step 4: Compute loss for each expert 
  weighted_loss_mat <- loss_potential_outcomes/p1
  y_experts <- sapply(1:dim(experts_test_01)[2], function(j) {
    # Element-wise multiply by weighted losses and average over everything
    mean(rowSums(experts_test_01[, j, ] * weighted_loss_mat))
  })
  
  # Step 5: Update the estimated cumulative loss for each expert (expert itself since t=1)
  # Step 6: Compute the new experts distribution 
  ## q2 = (q_{2,1}, ..., q_{2,J}) with q_{2,j} = exp(-y[,j])/(\sum_{i=1}^{J}exp(y[,i]))
  q2 <- exp(-eta*y_experts)/sum(exp(-eta*y_experts))
  return(q2)
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
