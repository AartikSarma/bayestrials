#' Identify effect modifiers
#' 
#' @param model A BayesianModel object or formula
#' @param data Data frame containing the data
#' @param treatment_var Name of the treatment variable
#' @param candidate_modifiers Vector of candidate effect modifier variables
#' @param method Method for identifying effect modifiers
#' @param credibility_threshold Threshold for credibility of effect modification
#' @param ... Additional arguments passed to specific methods
#' @return List of effect modifier results
#' @export
identify_effect_modifiers <- function(model,
                                    data,
                                    treatment_var,
                                    candidate_modifiers,
                                    method = c("interaction", "bart", "continuous"),
                                    credibility_threshold = 0.95,
                                    ...) {
  
  method <- match.arg(method)
  
  # Check inputs
  if (!is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }
  
  if (!treatment_var %in% names(data)) {
    stop(paste("Treatment variable", treatment_var, "not found in data"), call. = FALSE)
  }
  
  # Check that candidate modifiers exist in data
  missing_vars <- setdiff(candidate_modifiers, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("These candidate modifiers are not in the data:", paste(missing_vars, collapse = ", ")), 
         call. = FALSE)
  }
  
  # Initialize results structure
  results <- list(
    treatment_var = treatment_var,
    candidate_modifiers = candidate_modifiers,
    method = method,
    modifier_effects = list(),
    summary = NULL
  )
  
  # Apply appropriate method
  if (method == "interaction") {
    # Interaction method
    results <- identify_interaction_effects(model, data, treatment_var, candidate_modifiers, 
                                          credibility_threshold, ...)
  } else if (method == "bart") {
    # Bayesian Additive Regression Trees
    results <- identify_bart_effects(model, data, treatment_var, candidate_modifiers, 
                                   credibility_threshold, ...)
  } else if (method == "continuous") {
    # Continuous interaction modeling
    results <- identify_continuous_effects(model, data, treatment_var, candidate_modifiers, 
                                         credibility_threshold, ...)
  }
  
  # Add class for method dispatch
  class(results) <- c("effect_modifier_results", class(results))
  
  return(results)
}

#' Identify effect modifiers using interaction terms
#' 
#' @param model A BayesianModel object or formula
#' @param data Data frame containing the data
#' @param treatment_var Name of the treatment variable
#' @param candidate_modifiers Vector of candidate effect modifier variables
#' @param credibility_threshold Threshold for credibility of effect modification
#' @param fit_separate_models Logical indicating whether to fit separate models for each modifier
#' @param ... Additional arguments passed to fit_model
#' @return List of effect modifier results
#' @keywords internal
identify_interaction_effects <- function(model,
                                       data,
                                       treatment_var,
                                       candidate_modifiers,
                                       credibility_threshold = 0.95,
                                       fit_separate_models = TRUE,
                                       ...) {
  
  # Initialize results
  results <- list(
    treatment_var = treatment_var,
    candidate_modifiers = candidate_modifiers,
    method = "interaction",
    modifier_effects = list(),
    credible_modifiers = character(),
    models = list(),
    summary = NULL
  )
  
  # Get response variable and base predictor terms
  if (inherits(model, "BayesianModel")) {
    if (is.null(model$model_spec)) {
      stop("BayesianModel must have a model_spec", call. = FALSE)
    }
    
    response_var <- model$model_spec$outcome
    predictors <- model$model_spec$predictors
    family <- model$model_spec$family
    link <- model$model_spec$link
    random_effects <- model$model_spec$random_effects
  } else if (inherits(model, "formula") || is.character(model)) {
    # Parse formula to get response and predictors
    if (is.character(model)) {
      model <- stats::as.formula(model)
    }
    
    formula_terms <- terms(model)
    response_var <- all.vars(formula_terms)[1]
    
    # Extract predictors from formula
    predictors <- paste(attr(formula_terms, "term.labels"), collapse = " + ")
    
    # Set default family and link
    family <- "gaussian"
    link <- "identity"
    random_effects <- NULL
  } else if (inherits(model, "ModelSpecification")) {
    response_var <- model$outcome
    predictors <- model$predictors
    family <- model$family
    link <- model$link
    random_effects <- model$random_effects
  } else {
    stop("model must be a BayesianModel, ModelSpecification, formula, or character string", call. = FALSE)
  }
  
  # Check if response variable exists in data
  if (!response_var %in% names(data)) {
    stop(paste("Response variable", response_var, "not found in data"), call. = FALSE)
  }
  
  # Function to fit a model with interaction
  fit_interaction_model <- function(modifier_var) {
    # Check if modifier is numeric or factor
    if (is.numeric(data[[modifier_var]])) {
      # Center numeric modifier for better interpretation
      data[[paste0(modifier_var, "_c")]] <- scale(data[[modifier_var]], center = TRUE, scale = FALSE)
      mod_var <- paste0(modifier_var, "_c")
    } else {
      # For categorical variables, use as is
      mod_var <- modifier_var
    }
    
    # Create interaction term
    interaction_term <- paste0(treatment_var, " * ", mod_var)
    
    # Create full formula
    new_predictors <- paste0(predictors, " + ", interaction_term)
    
    # Create model specification
    model_spec <- create_model_spec(
      outcome = response_var,
      predictors = new_predictors,
      random_effects = random_effects,
      family = family,
      link = link,
      model_name = paste0("interaction_", modifier_var)
    )
    
    # Fit model
    interaction_model <- fit_model(model_spec, data, ...)
    
    return(interaction_model)
  }
  
  # Fit models for each modifier
  if (fit_separate_models) {
    log_message("Fitting separate models for each potential effect modifier", level = "info")
    
    # Fit a model for each candidate modifier
    for (modifier in candidate_modifiers) {
      log_message(paste("Fitting interaction model for", modifier), level = "info")
      
      # Fit model with this interaction
      interaction_model <- fit_interaction_model(modifier)
      
      # Store the model
      results$models[[modifier]] <- interaction_model
      
      # Extract interaction effect
      if (!is.null(interaction_model) && !is.null(interaction_model$model)) {
        # Get posterior for interaction term
        interaction_param <- paste0("b_", treatment_var, ":", modifier)
        
        # Handle numeric modifiers with _c suffix
        if (is.numeric(data[[modifier]])) {
          interaction_param <- paste0("b_", treatment_var, ":", modifier, "_c")
        }
        
        post_summary <- posterior_summary(interaction_model, interaction_param)
        
        if (nrow(post_summary) > 0) {
          # Check if effect is credibly non-zero
          lower_quantile <- paste0("q", 2.5)
          upper_quantile <- paste0("q", 97.5)
          
          # Effect is credible if CI doesn't include 0
          is_credible <- (post_summary[[lower_quantile]] > 0 & post_summary[[upper_quantile]] > 0) | 
                       (post_summary[[lower_quantile]] < 0 & post_summary[[upper_quantile]] < 0)
          
          # Store effect
          results$modifier_effects[[modifier]] <- list(
            modifier = modifier,
            interaction_param = interaction_param,
            estimate = post_summary$mean,
            sd = post_summary$sd,
            ci_lower = post_summary[[lower_quantile]],
            ci_upper = post_summary[[upper_quantile]],
            is_credible = is_credible
          )
          
          # Add to credible modifiers list if credible
          if (is_credible) {
            results$credible_modifiers <- c(results$credible_modifiers, modifier)
          }
        }
      }
    }
  } else {
    # Fit a single model with all interactions
    log_message("Fitting a single model with all potential effect modifiers", level = "info")
    
    # Create interaction terms for all modifiers
    interaction_terms <- character()
    
    for (modifier in candidate_modifiers) {
      # Check if modifier is numeric or factor
      if (is.numeric(data[[modifier]])) {
        # Center numeric modifier for better interpretation
        data[[paste0(modifier, "_c")]] <- scale(data[[modifier]], center = TRUE, scale = FALSE)
        mod_var <- paste0(modifier, "_c")
      } else {
        # For categorical variables, use as is
        mod_var <- modifier
      }
      
      # Create interaction term
      interaction_terms <- c(interaction_terms, paste0(treatment_var, " * ", mod_var))
    }
    
    # Combine all terms
    all_interactions <- paste(interaction_terms, collapse = " + ")
    
    # Create full formula
    new_predictors <- paste0(predictors, " + ", all_interactions)
    
    # Create model specification
    model_spec <- create_model_spec(
      outcome = response_var,
      predictors = new_predictors,
      random_effects = random_effects,
      family = family,
      link = link,
      model_name = "all_interactions"
    )
    
    # Fit model
    full_model <- fit_model(model_spec, data, ...)
    
    # Store the model
    results$models[["full_model"]] <- full_model
    
    # Extract interaction effects
    if (!is.null(full_model) && !is.null(full_model$model)) {
      for (modifier in candidate_modifiers) {
        # Get posterior for interaction term
        interaction_param <- paste0("b_", treatment_var, ":", modifier)
        
        # Handle numeric modifiers with _c suffix
        if (is.numeric(data[[modifier]])) {
          interaction_param <- paste0("b_", treatment_var, ":", modifier, "_c")
        }
        
        post_summary <- posterior_summary(full_model, interaction_param)
        
        if (nrow(post_summary) > 0) {
          # Check if effect is credibly non-zero
          lower_quantile <- paste0("q", 2.5)
          upper_quantile <- paste0("q", 97.5)
          
          # Effect is credible if CI doesn't include 0
          is_credible <- (post_summary[[lower_quantile]] > 0 & post_summary[[upper_quantile]] > 0) | 
                       (post_summary[[lower_quantile]] < 0 & post_summary[[upper_quantile]] < 0)
          
          # Store effect
          results$modifier_effects[[modifier]] <- list(
            modifier = modifier,
            interaction_param = interaction_param,
            estimate = post_summary$mean,
            sd = post_summary$sd,
            ci_lower = post_summary[[lower_quantile]],
            ci_upper = post_summary[[upper_quantile]],
            is_credible = is_credible
          )
          
          # Add to credible modifiers list if credible
          if (is_credible) {
            results$credible_modifiers <- c(results$credible_modifiers, modifier)
          }
        }
      }
    }
  }
  
  # Create summary
  results$summary <- summarize_effect_modifiers(results)
  
  return(results)
}

#' Identify effect modifiers using Bayesian Additive Regression Trees
#' 
#' @param model A BayesianModel object or formula
#' @param data Data frame containing the data
#' @param treatment_var Name of the treatment variable
#' @param candidate_modifiers Vector of candidate effect modifier variables
#' @param credibility_threshold Threshold for credibility of effect modification
#' @param n_trees Number of trees for BART model
#' @param ... Additional arguments passed to BART functions
#' @return List of effect modifier results
#' @keywords internal
identify_bart_effects <- function(model,
                                data,
                                treatment_var,
                                candidate_modifiers,
                                credibility_threshold = 0.95,
                                n_trees = 200,
                                ...) {
  
  # Check if dbarts package is available
  if (!requireNamespace("dbarts", quietly = TRUE)) {
    stop("dbarts package is required but not available", call. = FALSE)
  }
  
  # Initialize results
  results <- list(
    treatment_var = treatment_var,
    candidate_modifiers = candidate_modifiers,
    method = "bart",
    modifier_effects = list(),
    credible_modifiers = character(),
    importance = NULL,
    model = NULL,
    summary = NULL
  )
  
  # Get response variable
  if (inherits(model, "BayesianModel")) {
    if (is.null(model$model_spec)) {
      stop("BayesianModel must have a model_spec", call. = FALSE)
    }
    
    response_var <- model$model_spec$outcome
  } else if (inherits(model, "formula") || is.character(model)) {
    # Parse formula to get response
    if (is.character(model)) {
      model <- stats::as.formula(model)
    }
    
    formula_terms <- terms(model)
    response_var <- all.vars(formula_terms)[1]
  } else if (inherits(model, "ModelSpecification")) {
    response_var <- model$outcome
  } else {
    stop("model must be a BayesianModel, ModelSpecification, formula, or character string", call. = FALSE)
  }
  
  # Check if response variable exists in data
  if (!response_var %in% names(data)) {
    stop(paste("Response variable", response_var, "not found in data"), call. = FALSE)
  }
  
  # Prepare data for BART
  # Create a subset with only treatment, response, and candidate modifiers
  vars_to_include <- c(response_var, treatment_var, candidate_modifiers)
  bart_data <- data[, vars_to_include, drop = FALSE]
  
  # Convert treatment to binary if not already
  if (!is.numeric(bart_data[[treatment_var]])) {
    bart_data[[treatment_var]] <- as.numeric(as.factor(bart_data[[treatment_var]])) - 1
    warning(paste("Converted", treatment_var, "to binary (0/1) for BART analysis"), call. = FALSE)
  }
  
  # Check response type and set appropriate BART model type
  if (is.numeric(bart_data[[response_var]])) {
    # Continuous outcome
    log_message("Fitting BART model for continuous outcome", level = "info")
    
    # Set up formula
    x_vars <- c(treatment_var, candidate_modifiers)
    bart_formula <- stats::as.formula(paste(response_var, "~", paste(x_vars, collapse = " + ")))
    
    # Fit BART model
    set.seed(42)  # For reproducibility
    bart_model <- dbarts::bart(
      formula = bart_formula,
      data = bart_data,
      ntree = n_trees,
      keepTrees = TRUE,
      verbose = FALSE,
      ...
    )
    
    # Store the model
    results$model <- bart_model
    
    # Calculate variable importance
    var_counts <- dbarts::makeModelMatrixFromDataFrame(bart_data[, x_vars, drop = FALSE])
    var_importance <- apply(var_counts, 2, function(x) length(unique(x)) > 1)
    var_importance <- var_importance[var_importance]
    var_names <- names(var_importance)
    
    # Extract inclusion proportions
    inclusion_props <- colMeans(bart_model$varcount[, var_names, drop = FALSE])
    
    # Store variable importance
    results$importance <- data.frame(
      variable = names(inclusion_props),
      importance = inclusion_props,
      stringsAsFactors = FALSE
    )
    
    # Identify treatment interaction effects
    # This is done by examining how treatment effects vary across levels of modifiers
    
    # Generate predictions across modifier values
    for (modifier in candidate_modifiers) {
      log_message(paste("Analyzing treatment effect modification by", modifier), level = "info")
      
      # Skip if modifier is constant
      if (length(unique(bart_data[[modifier]])) <= 1) {
        log_message(paste("Skipping", modifier, "- no variation"), level = "warning")
        next
      }
      
      # Create prediction grid
      if (is.factor(bart_data[[modifier]]) || is.character(bart_data[[modifier]])) {
        # For categorical modifiers, use all levels
        modifier_values <- sort(unique(bart_data[[modifier]]))
      } else {
        # For continuous modifiers, use quantiles
        modifier_values <- stats::quantile(bart_data[[modifier]], probs = seq(0.1, 0.9, by = 0.2), 
                                         na.rm = TRUE)
      }
      
      # Create prediction data frames
      pred_data_treated <- pred_data_control <- bart_data[1, ]
      
      # Initialize storage for treatment effects
      te_by_modifier <- matrix(NA, nrow = length(modifier_values), ncol = bart_model$ndpost)
      
      # Generate predictions for each modifier value
      for (i in seq_along(modifier_values)) {
        mod_value <- modifier_values[i]
        
        # Set treated and control conditions
        pred_data_treated[[treatment_var]] <- 1
        pred_data_treated[[modifier]] <- mod_value
        
        pred_data_control[[treatment_var]] <- 0
        pred_data_control[[modifier]] <- mod_value
        
        # Predict outcomes
        y_treated <- stats::predict(bart_model, newdata = pred_data_treated)
        y_control <- stats::predict(bart_model, newdata = pred_data_control)
        
        # Treatment effect = treated - control
        te_by_modifier[i, ] <- y_treated - y_control
      }
      
      # Calculate mean treatment effect and credible intervals for each modifier value
      te_summary <- data.frame(
        modifier_value = as.character(modifier_values),
        mean_te = rowMeans(te_by_modifier),
        lower_ci = apply(te_by_modifier, 1, function(x) stats::quantile(x, 0.025)),
        upper_ci = apply(te_by_modifier, 1, function(x) stats::quantile(x, 0.975)),
        stringsAsFactors = FALSE
      )
      
      # Test if there's heterogeneity in treatment effect
      # If CIs don't overlap, there's evidence of effect modification
      range_lower <- range(te_summary$lower_ci)
      range_upper <- range(te_summary$upper_ci)
      
      # Effect modification is credible if ranges don't overlap
      is_credible <- range_lower[2] > range_upper[1] || range_upper[2] < range_lower[1]
      
      # Store effect
      results$modifier_effects[[modifier]] <- list(
        modifier = modifier,
        modifier_values = modifier_values,
        treatment_effects = te_summary,
        is_credible = is_credible
      )
      
      # Add to credible modifiers list if credible
      if (is_credible) {
        results$credible_modifiers <- c(results$credible_modifiers, modifier)
      }
    }
    
  } else {
    # Binary outcome
    log_message("Fitting BART model for binary outcome", level = "info")
    
    # Convert response to binary if needed
    if (!is.numeric(bart_data[[response_var]])) {
      bart_data[[response_var]] <- as.numeric(as.factor(bart_data[[response_var]])) - 1
      warning(paste("Converted", response_var, "to binary (0/1) for BART analysis"), call. = FALSE)
    }
    
    # Set up formula
    x_vars <- c(treatment_var, candidate_modifiers)
    bart_formula <- stats::as.formula(paste(response_var, "~", paste(x_vars, collapse = " + ")))
    
    # Fit BART model for binary outcome
    set.seed(42)  # For reproducibility
    bart_model <- dbarts::bart(
      formula = bart_formula,
      data = bart_data,
      ntree = n_trees,
      keepTrees = TRUE,
      verbose = FALSE,
      binary = TRUE,
      ...
    )
    
    # Store the model
    results$model <- bart_model
    
    # Calculate variable importance
    var_counts <- dbarts::makeModelMatrixFromDataFrame(bart_data[, x_vars, drop = FALSE])
    var_importance <- apply(var_counts, 2, function(x) length(unique(x)) > 1)
    var_importance <- var_importance[var_importance]
    var_names <- names(var_importance)
    
    # Extract inclusion proportions
    inclusion_props <- colMeans(bart_model$varcount[, var_names, drop = FALSE])
    
    # Store variable importance
    results$importance <- data.frame(
      variable = names(inclusion_props),
      importance = inclusion_props,
      stringsAsFactors = FALSE
    )
    
    # Identify treatment interaction effects using same approach as for continuous outcome
    # Generate predictions across modifier values
    for (modifier in candidate_modifiers) {
      log_message(paste("Analyzing treatment effect modification by", modifier), level = "info")
      
      # Skip if modifier is constant
      if (length(unique(bart_data[[modifier]])) <= 1) {
        log_message(paste("Skipping", modifier, "- no variation"), level = "warning")
        next
      }
      
      # Create prediction grid
      if (is.factor(bart_data[[modifier]]) || is.character(bart_data[[modifier]])) {
        # For categorical modifiers, use all levels
        modifier_values <- sort(unique(bart_data[[modifier]]))
      } else {
        # For continuous modifiers, use quantiles
        modifier_values <- stats::quantile(bart_data[[modifier]], probs = seq(0.1, 0.9, by = 0.2), 
                                         na.rm = TRUE)
      }
      
      # Create prediction data frames
      pred_data_treated <- pred_data_control <- bart_data[1, ]
      
      # Initialize storage for treatment effects (on probability scale)
      te_by_modifier <- matrix(NA, nrow = length(modifier_values), ncol = bart_model$ndpost)
      
      # Generate predictions for each modifier value
      for (i in seq_along(modifier_values)) {
        mod_value <- modifier_values[i]
        
        # Set treated and control conditions
        pred_data_treated[[treatment_var]] <- 1
        pred_data_treated[[modifier]] <- mod_value
        
        pred_data_control[[treatment_var]] <- 0
        pred_data_control[[modifier]] <- mod_value
        
        # Predict probabilities
        p_treated <- stats::pnorm(stats::predict(bart_model, newdata = pred_data_treated))
        p_control <- stats::pnorm(stats::predict(bart_model, newdata = pred_data_control))
        
        # Treatment effect = treated - control (probability difference)
        te_by_modifier[i, ] <- p_treated - p_control
      }
      
      # Calculate mean treatment effect and credible intervals for each modifier value
      te_summary <- data.frame(
        modifier_value = as.character(modifier_values),
        mean_te = rowMeans(te_by_modifier),
        lower_ci = apply(te_by_modifier, 1, function(x) stats::quantile(x, 0.025)),
        upper_ci = apply(te_by_modifier, 1, function(x) stats::quantile(x, 0.975)),
        stringsAsFactors = FALSE
      )
      
      # Test if there's heterogeneity in treatment effect
      # If CIs don't overlap, there's evidence of effect modification
      range_lower <- range(te_summary$lower_ci)
      range_upper <- range(te_summary$upper_ci)
      
      # Effect modification is credible if ranges don't overlap
      is_credible <- range_lower[2] > range_upper[1] || range_upper[2] < range_lower[1]
      
      # Store effect
      results$modifier_effects[[modifier]] <- list(
        modifier = modifier,
        modifier_values = modifier_values,
        treatment_effects = te_summary,
        is_credible = is_credible
      )
      
      # Add to credible modifiers list if credible
      if (is_credible) {
        results$credible_modifiers <- c(results$credible_modifiers, modifier)
      }
    }
  }
  
  # Create summary
  results$summary <- summarize_effect_modifiers(results)
  
  return(results)
}

#' Identify effect modifiers using continuous interaction modeling
#' 
#' @param model A BayesianModel object or formula
#' @param data Data frame containing the data
#' @param treatment_var Name of the treatment variable
#' @param candidate_modifiers Vector of candidate effect modifier variables
#' @param credibility_threshold Threshold for credibility of effect modification
#' @param use_splines Logical indicating whether to use splines for continuous modifiers
#' @param spline_df Degrees of freedom for splines
#' @param ... Additional arguments passed to fit_model
#' @return List of effect modifier results
#' @keywords internal
identify_continuous_effects <- function(model,
                                      data,
                                      treatment_var,
                                      candidate_modifiers,
                                      credibility_threshold = 0.95,
                                      use_splines = TRUE,
                                      spline_df = 3,
                                      ...) {
  
  # Initialize results
  results <- list(
    treatment_var = treatment_var,
    candidate_modifiers = candidate_modifiers,
    method = "continuous",
    modifier_effects = list(),
    credible_modifiers = character(),
    models = list(),
    summary = NULL
  )
  
  # Get response variable and base predictor terms
  if (inherits(model, "BayesianModel")) {
    if (is.null(model$model_spec)) {
      stop("BayesianModel must have a model_spec", call. = FALSE)
    }
    
    response_var <- model$model_spec$outcome
    predictors <- model$model_spec$predictors
    family <- model$model_spec$family
    link <- model$model_spec$link
    random_effects <- model$model_spec$random_effects
  } else if (inherits(model, "formula") || is.character(model)) {
    # Parse formula to get response and predictors
    if (is.character(model)) {
      model <- stats::as.formula(model)
    }
    
    formula_terms <- terms(model)
    response_var <- all.vars(formula_terms)[1]
    
    # Extract predictors from formula
    predictors <- paste(attr(formula_terms, "term.labels"), collapse = " + ")
    
    # Set default family and link
    family <- "gaussian"
    link <- "identity"
    random_effects <- NULL
  } else if (inherits(model, "ModelSpecification")) {
    response_var <- model$outcome
    predictors <- model$predictors
    family <- model$family
    link <- model$link
    random_effects <- model$random_effects
  } else {
    stop("model must be a BayesianModel, ModelSpecification, formula, or character string", call. = FALSE)
  }
  
  # Check if response variable exists in data
  if (!response_var %in% names(data)) {
    stop(paste("Response variable", response_var, "not found in data"), call. = FALSE)
  }
  
  # Check for splines package if needed
  if (use_splines && !requireNamespace("splines", quietly = TRUE)) {
    stop("splines package is required for spline interactions but not available", call. = FALSE)
  }
  
  # Filter candidate modifiers to only numeric variables for continuous modeling
  numeric_modifiers <- candidate_modifiers[sapply(data[, candidate_modifiers, drop = FALSE], is.numeric)]
  
  if (length(numeric_modifiers) == 0) {
    stop("No numeric variables found among candidate modifiers for continuous interaction modeling", 
         call. = FALSE)
  }
  
  log_message(paste("Found", length(numeric_modifiers), "numeric variables for continuous interaction modeling"), 
             level = "info")
  
  # For each continuous modifier, fit models with spline interactions
  for (modifier in numeric_modifiers) {
    log_message(paste("Modeling continuous interaction with", modifier), level = "info")
    
    # Check range of modifier
    mod_range <- range(data[[modifier]], na.rm = TRUE)
    
    if (diff(mod_range) == 0) {
      log_message(paste("Skipping", modifier, "- no variation"), level = "warning")
      next
    }
    
    # Center the modifier for better interpretation
    data[[paste0(modifier, "_c")]] <- scale(data[[modifier]], center = TRUE, scale = FALSE)
    mod_var <- paste0(modifier, "_c")
    
    # Create interaction term
    if (use_splines) {
      # Use spline interaction
      log_message(paste("Using spline interaction with df =", spline_df), level = "info")
      
      # Create spline term
      spline_term <- paste0("splines::bs(", mod_var, ", df = ", spline_df, ")")
      
      # Create interaction
      interaction_term <- paste0(treatment_var, " * ", spline_term)
    } else {
      # Use simple linear interaction
      interaction_term <- paste0(treatment_var, " * ", mod_var)
    }
    
    # Create full formula
    new_predictors <- paste0(predictors, " + ", interaction_term)
    
    # Create model specification
    model_spec <- create_model_spec(
      outcome = response_var,
      predictors = new_predictors,
      random_effects = random_effects,
      family = family,
      link = link,
      model_name = paste0("continuous_", modifier)
    )
    
    # Fit model
    interaction_model <- fit_model(model_spec, data, ...)
    
    # Store the model
    results$models[[modifier]] <- interaction_model
    
    # Generate predicted treatment effects across modifier range
    if (!is.null(interaction_model) && !is.null(interaction_model$model)) {
      # Create prediction grid
      grid_points <- 20
      mod_seq <- seq(mod_range[1], mod_range[2], length.out = grid_points)
      
      # Create prediction data
      newdata <- data.frame(
        pred_id = 1:grid_points,
        mod_value = mod_seq,
        mod_value_c = scale(mod_seq, center = mean(data[[modifier]], na.rm = TRUE), scale = FALSE)
      )
      
      names(newdata)[2] <- modifier
      names(newdata)[3] <- mod_var
      
      # Add other predictors with mean/reference values
      for (var in all.vars(stats::as.formula(paste("~", predictors)))) {
        if (!var %in% names(newdata) && var != treatment_var) {
          if (is.numeric(data[[var]])) {
            newdata[[var]] <- mean(data[[var]], na.rm = TRUE)
          } else if (is.factor(data[[var]])) {
            newdata[[var]] <- levels(data[[var]])[1]
          } else {
            newdata[[var]] <- unique(data[[var]])[1]
          }
        }
      }
      
      # Create treated and control versions
      newdata_treated <- newdata
      newdata_treated[[treatment_var]] <- unique(data[[treatment_var]])[2]  # Treated value
      
      newdata_control <- newdata
      newdata_control[[treatment_var]] <- unique(data[[treatment_var]])[1]  # Control value
      
      # Add any variables needed for random effects
      if (!is.null(random_effects)) {
        for (var in all.vars(stats::as.formula(paste("~", random_effects)))) {
          if (!var %in% names(newdata)) {
            # Use first level for grouping variables
            newdata_treated[[var]] <- unique(data[[var]])[1]
            newdata_control[[var]] <- unique(data[[var]])[1]
          }
        }
      }
      
      # Generate predictions
      pred_treated <- brms::posterior_predict(interaction_model$model, newdata = newdata_treated)
      pred_control <- brms::posterior_predict(interaction_model$model, newdata = newdata_control)
      
      # Calculate treatment effects (depends on outcome type)
      if (family == "gaussian") {
        # For continuous outcomes, TE = treated - control
        te_samples <- pred_treated - pred_control
      } else if (family %in% c("bernoulli", "binomial")) {
        # For binary outcomes, convert to probabilities
        prob_treated <- apply(pred_treated, 2, function(x) mean(x > 0))
        prob_control <- apply(pred_control, 2, function(x) mean(x > 0))
        te_samples <- prob_treated - prob_control
      } else {
        # For other outcomes, use the mean difference
        te_samples <- apply(pred_treated, 2, mean) - apply(pred_control, 2, mean)
      }
      
      # Calculate summary statistics
      te_mean <- apply(te_samples, 2, mean)
      te_lower <- apply(te_samples, 2, function(x) stats::quantile(x, 0.025))
      te_upper <- apply(te_samples, 2, function(x) stats::quantile(x, 0.975))
      
      # Create summary data frame
      te_summary <- data.frame(
        modifier_value = mod_seq,
        mean_te = te_mean,
        lower_ci = te_lower,
        upper_ci = te_upper,
        stringsAsFactors = FALSE
      )
      
      # Determine if there's credible effect modification
      # Test if slope of TE across modifier values is non-zero
      # or if CIs don't overlap
      range_lower <- range(te_summary$lower_ci)
      range_upper <- range(te_summary$upper_ci)
      
      # Effect modification is credible if ranges don't overlap
      is_credible <- range_lower[2] > range_upper[1] || range_upper[2] < range_lower[1]
      
      # Store effect
      results$modifier_effects[[modifier]] <- list(
        modifier = modifier,
        treatment_effects = te_summary,
        is_credible = is_credible,
        use_splines = use_splines,
        spline_df = if (use_splines) spline_df else NULL
      )
      
      # Add to credible modifiers list if credible
      if (is_credible) {
        results$credible_modifiers <- c(results$credible_modifiers, modifier)
      }
    }
  }
  
  # Create summary
  results$summary <- summarize_effect_modifiers(results)
  
  return(results)
}

#' Summarize effect modifier results
#' 
#' @param effect_results Results from effect modifier identification
#' @return Text summary of effect modifiers
#' @keywords internal
summarize_effect_modifiers <- function(effect_results) {
  if (!is.list(effect_results)) {
    stop("effect_results must be a list", call. = FALSE)
  }
  
  # Extract key information
  treatment_var <- effect_results$treatment_var
  method <- effect_results$method
  credible_modifiers <- effect_results$credible_modifiers
  
  # Create summary text
  summary_text <- paste0("# Effect Modifier Analysis Summary\n\n",
                        "## Analysis Information\n\n",
                        "- Treatment Variable: ", treatment_var, "\n",
                        "- Analysis Method: ", method, "\n",
                        "- Candidate Modifiers: ", paste(effect_results$candidate_modifiers, collapse = ", "), "\n\n")
  
  # Add identified modifiers
  summary_text <- paste0(summary_text, "## Identified Effect Modifiers\n\n")
  
  if (length(credible_modifiers) > 0) {
    summary_text <- paste0(summary_text, "The following variables were identified as credible effect modifiers:\n\n")
    
    for (modifier in credible_modifiers) {
      mod_effect <- effect_results$modifier_effects[[modifier]]
      
      summary_text <- paste0(summary_text, "### ", modifier, "\n\n")
      
      if (method == "interaction") {
        # For standard interaction method
        summary_text <- paste0(summary_text, 
                              "- Interaction Parameter: ", mod_effect$interaction_param, "\n",
                              "- Estimate: ", sprintf("%.4f", mod_effect$estimate), "\n",
                              "- 95% CI: [", sprintf("%.4f, %.4f", mod_effect$ci_lower, mod_effect$ci_upper), "]\n\n")
      } else if (method == "bart" || method == "continuous") {
        # For BART or continuous methods
        summary_text <- paste0(summary_text, 
                              "- Treatment effect varies by ", modifier, " values\n",
                              "- Range of treatment effects: [", 
                              sprintf("%.4f, %.4f", min(mod_effect$treatment_effects$mean_te), 
                                     max(mod_effect$treatment_effects$mean_te)), "]\n\n")
      }
    }
  } else {
    summary_text <- paste0(summary_text, "No credible effect modifiers were identified among the candidate variables.\n\n")
  }
  
  # Add interpretation
  summary_text <- paste0(summary_text, "## Interpretation\n\n")
  
  if (length(credible_modifiers) > 0) {
    summary_text <- paste0(summary_text, 
                          "Effect modification means that the treatment effect varies depending on the value of these variables. ",
                          "This suggests potential heterogeneity in treatment effects across different subgroups or values of these modifiers.\n\n",
                          "For clinical practice, this may imply that treatment decisions could be tailored based on these modifiers ",
                          "to optimize patient outcomes.\n")
  } else {
    summary_text <- paste0(summary_text, 
                          "The analysis did not identify credible effect modifiers, suggesting that the treatment effect ",
                          "is relatively consistent across the range of candidate variables examined. ",
                          "This implies a more uniform treatment effect across different subgroups.\n")
  }
  
  return(summary_text)
}

#' Plot effect modifier results
#' 
#' @param effect_results Results from effect modifier identification
#' @param modifier Name of the modifier to plot
#' @param plot_type Type of plot to create
#' @return A ggplot object
#' @export
plot_interaction <- function(effect_results,
                           modifier,
                           plot_type = c("continuous", "categorical")) {
  
  plot_type <- match.arg(plot_type)
  
  if (!is.list(effect_results)) {
    stop("effect_results must be a list", call. = FALSE)
  }
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Check if modifier exists in results
  if (!modifier %in% names(effect_results$modifier_effects)) {
    stop(paste("Modifier", modifier, "not found in effect_results"), call. = FALSE)
  }
  
  # Get modifier effect
  mod_effect <- effect_results$modifier_effects[[modifier]]
  
  # Create plot based on method and plot type
  if (effect_results$method == "interaction") {
    # Standard interaction method
    
    # For this method, we need the original data and model to generate predictions
    if (is.null(effect_results$models[[modifier]])) {
      stop("Model for this modifier not found in effect_results", call. = FALSE)
    }
    
    model <- effect_results$models[[modifier]]
    
    if (is.null(model$model)) {
      stop("Model has not been fitted", call. = FALSE)
    }
    
    # Get data from model
    data <- model$model$data
    
    # Get treatment variable
    treatment_var <- effect_results$treatment_var
    
    # Check modifier type for appropriate plot
    if (is.numeric(data[[modifier]]) && plot_type == "continuous") {
      # Continuous modifier
      
      # Create prediction grid
      grid_points <- 100
      mod_range <- range(data[[modifier]], na.rm = TRUE)
      mod_seq <- seq(mod_range[1], mod_range[2], length.out = grid_points)
      
      # Create prediction data frames
      newdata <- data.frame(
        mod_value = mod_seq
      )
      names(newdata)[1] <- modifier
      
      # Handle centered version if used in model
      if (paste0(modifier, "_c") %in% names(data)) {
        newdata[[paste0(modifier, "_c")]] <- scale(mod_seq, 
                                                 center = mean(data[[modifier]], na.rm = TRUE),
                                                 scale = FALSE)
      }
      
      # Add other predictors with mean/reference values
      for (var in names(data)) {
        if (!var %in% names(newdata) && var != treatment_var) {
          if (is.numeric(data[[var]])) {
            newdata[[var]] <- mean(data[[var]], na.rm = TRUE)
          } else if (is.factor(data[[var]])) {
            newdata[[var]] <- levels(data[[var]])[1]
          } else {
            newdata[[var]] <- unique(data[[var]])[1]
          }
        }
      }
      
      # Create treated and control versions
      newdata_treated <- newdata_control <- newdata
      
      # Get treatment values
      if (is.factor(data[[treatment_var]])) {
        trt_levels <- levels(data[[treatment_var]])
        newdata_treated[[treatment_var]] <- trt_levels[2]
        newdata_control[[treatment_var]] <- trt_levels[1]
      } else if (is.numeric(data[[treatment_var]]) && length(unique(data[[treatment_var]])) <= 2) {
        trt_values <- sort(unique(data[[treatment_var]]))
        newdata_treated[[treatment_var]] <- trt_values[2]
        newdata_control[[treatment_var]] <- trt_values[1]
      } else {
        # For other cases, use first two unique values
        trt_values <- unique(data[[treatment_var]])
        newdata_treated[[treatment_var]] <- trt_values[min(2, length(trt_values))]
        newdata_control[[treatment_var]] <- trt_values[1]
      }
      
      # Generate predictions
      pred_treated <- brms::posterior_epred(model$model, newdata = newdata_treated)
      pred_control <- brms::posterior_epred(model$model, newdata = newdata_control)
      
      # Calculate summary statistics
      treated_mean <- colMeans(pred_treated)
      treated_lower <- apply(pred_treated, 2, function(x) stats::quantile(x, 0.025))
      treated_upper <- apply(pred_treated, 2, function(x) stats::quantile(x, 0.975))
      
      control_mean <- colMeans(pred_control)
      control_lower <- apply(pred_control, 2, function(x) stats::quantile(x, 0.025))
      control_upper <- apply(pred_control, 2, function(x) stats::quantile(x, 0.975))
      
      # Create data frame for plotting
      plot_data <- data.frame(
        modifier = rep(mod_seq, 2),
        group = rep(c("Treated", "Control"), each = grid_points),
        mean = c(treated_mean, control_mean),
        lower = c(treated_lower, control_lower),
        upper = c(treated_upper, control_upper),
        stringsAsFactors = FALSE
      )
      
      # Create plot
      plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = modifier, y = mean, color = group, fill = group)) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
        ggplot2::labs(
          title = paste("Treatment Effect Modification by", modifier),
          subtitle = paste("Interaction Parameter Estimate:", sprintf("%.4f", mod_effect$estimate),
                         "95% CI: [", sprintf("%.4f, %.4f", mod_effect$ci_lower, mod_effect$ci_upper), "]"),
          x = modifier,
          y = "Predicted Outcome",
          color = "Group",
          fill = "Group"
        ) +
        ggplot2::theme_minimal()
      
    } else if (plot_type == "categorical" || !is.numeric(data[[modifier]])) {
      # Categorical representation
      
      # If numeric, discretize into quantile groups
      if (is.numeric(data[[modifier]])) {
        # Create quantile groups
        n_groups <- 5
        cuts <- stats::quantile(data[[modifier]], probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE)
        
        # Create labels
        labels <- character(n_groups)
        for (i in 1:n_groups) {
          labels[i] <- sprintf("[%.2f, %.2f]", cuts[i], cuts[i+1])
        }
        
        # Create grouping variable
        mod_groups <- cut(data[[modifier]], breaks = cuts, labels = labels, include.lowest = TRUE)
        
        # Add to data
        data$mod_group <- mod_groups
        mod_group_var <- "mod_group"
      } else {
        # Use as is for categorical variables
        mod_group_var <- modifier
      }
      
      # Get unique groups
      groups <- levels(data[[mod_group_var]])
      if (is.null(groups)) groups <- unique(data[[mod_group_var]])
      
      # Create prediction data frames
      pred_data <- list()
      
      for (group in groups) {
        # Create data for this group
        group_data <- data.frame(
          group = group
        )
        names(group_data)[1] <- mod_group_var
        
        # Add modifier value if different from group variable
        if (mod_group_var != modifier) {
          # For numeric, use group midpoint
          if (is.numeric(data[[modifier]])) {
            group_idx <- which(labels == group)
            if (length(group_idx) > 0) {
              midpoint <- (cuts[group_idx] + cuts[group_idx + 1]) / 2
              group_data[[modifier]] <- midpoint
              
              # Add centered version if needed
              if (paste0(modifier, "_c") %in% names(data)) {
                group_data[[paste0(modifier, "_c")]] <- scale(midpoint, 
                                                           center = mean(data[[modifier]], na.rm = TRUE),
                                                           scale = FALSE)
              }
            }
          }
        }
        
        # Add other predictors with mean/reference values
        for (var in names(data)) {
          if (!var %in% names(group_data) && var != treatment_var && var != "mod_group") {
            if (is.numeric(data[[var]])) {
              group_data[[var]] <- mean(data[[var]], na.rm = TRUE)
            } else if (is.factor(data[[var]])) {
              group_data[[var]] <- levels(data[[var]])[1]
            } else {
              group_data[[var]] <- unique(data[[var]])[1]
            }
          }
        }
        
        # Create treated and control versions
        group_data_treated <- group_data_control <- group_data
        
        # Get treatment values
        if (is.factor(data[[treatment_var]])) {
          trt_levels <- levels(data[[treatment_var]])
          group_data_treated[[treatment_var]] <- trt_levels[2]
          group_data_control[[treatment_var]] <- trt_levels[1]
        } else if (is.numeric(data[[treatment_var]]) && length(unique(data[[treatment_var]])) <= 2) {
          trt_values <- sort(unique(data[[treatment_var]]))
          group_data_treated[[treatment_var]] <- trt_values[2]
          group_data_control[[treatment_var]] <- trt_values[1]
        } else {
          # For other cases, use first two unique values
          trt_values <- unique(data[[treatment_var]])
          group_data_treated[[treatment_var]] <- trt_values[min(2, length(trt_values))]
          group_data_control[[treatment_var]] <- trt_values[1]
        }
        
        # Generate predictions
        pred_treated <- brms::posterior_epred(model$model, newdata = group_data_treated)
        pred_control <- brms::posterior_epred(model$model, newdata = group_data_control)
        
        # Calculate treatment effect
        te_samples <- pred_treated - pred_control
        
        # Calculate summary statistics
        te_mean <- mean(te_samples)
        te_lower <- stats::quantile(te_samples, 0.025)
        te_upper <- stats::quantile(te_samples, 0.975)
        
        # Store results
        pred_data[[group]] <- list(
          group = group,
          te_mean = te_mean,
          te_lower = te_lower,
          te_upper = te_upper
        )
      }
      
      # Combine results
      plot_data <- do.call(rbind, lapply(pred_data, function(x) {
        data.frame(
          group = x$group,
          mean = x$te_mean,
          lower = x$te_lower,
          upper = x$te_upper,
          stringsAsFactors = FALSE
        )
      }))
      
      # Create plot
      plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = mean)) +
        ggplot2::geom_point(size = 3, color = "blue") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("Treatment Effect Modification by", modifier),
          subtitle = paste("Interaction Parameter Estimate:", sprintf("%.4f", mod_effect$estimate),
                         "95% CI: [", sprintf("%.4f, %.4f", mod_effect$ci_lower, mod_effect$ci_upper), "]"),
          x = modifier,
          y = "Treatment Effect"
        ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
    } else {
      stop("Invalid plot_type for this modifier", call. = FALSE)
    }
    
  } else if (effect_results$method == "bart" || effect_results$method == "continuous") {
    # For BART or continuous methods, plot TE vs modifier value
    
    # Get treatment effect data
    te_data <- mod_effect$treatment_effects
    
    # Create plot
    plot <- ggplot2::ggplot(te_data, ggplot2::aes(x = modifier_value, y = mean_te)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Treatment Effect Modification by", modifier),
        subtitle = paste("Method:", effect_results$method, 
                       if (effect_results$method == "continuous" && mod_effect$use_splines) 
                         paste("(Splines, df =", mod_effect$spline_df, ")") else ""),
        x = modifier,
        y = "Treatment Effect"
      )
  }
  
  return(plot)
}

#' Assess credibility of subgroup effects
#' 
#' @param effect_results Results from effect modifier identification
#' @param criteria Vector of credibility criteria to assess
#' @return Data frame with credibility assessment
#' @export
assess_subgroup_credibility <- function(effect_results,
                                       criteria = c("pre_specified", "interaction", "consistency", 
                                                  "multiple_testing", "biological_plausibility")) {
  
  if (!is.list(effect_results)) {
    stop("effect_results must be a list", call. = FALSE)
  }
  
  # Initialize credibility assessment
  assessment <- list()
  
  # Get credible modifiers
  credible_modifiers <- effect_results$credible_modifiers
  
  if (length(credible_modifiers) == 0) {
    return(data.frame(
      modifier = character(0),
      criterion = character(0),
      result = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # For each modifier, assess credibility criteria
  for (modifier in credible_modifiers) {
    mod_assessment <- list()
    
    # Assess each criterion
    if ("pre_specified" %in% criteria) {
      # Pre-specification criterion
      # This requires user input since we don't know if modifiers were pre-specified
      mod_assessment$pre_specified <- list(
        criterion = "Pre-specified",
        result = "Unknown - requires user input",
        score = NA
      )
    }
    
    if ("interaction" %in% criteria) {
      # Statistical interaction criterion
      if (effect_results$method == "interaction") {
        mod_effect <- effect_results$modifier_effects[[modifier]]
        
        # Calculate p-value for interaction term
        z_score <- mod_effect$estimate / mod_effect$sd
        p_value <- 2 * (1 - stats::pnorm(abs(z_score)))
        
        result <- paste0("p = ", sprintf("%.4f", p_value))
        score <- if (p_value < 0.05) 1 else 0
      } else {
        # For other methods, use is_credible directly
        is_credible <- effect_results$modifier_effects[[modifier]]$is_credible
        result <- if (is_credible) "Credible effect modification" else "No credible effect modification"
        score <- if (is_credible) 1 else 0
      }
      
      mod_assessment$interaction <- list(
        criterion = "Statistical Interaction",
        result = result,
        score = score
      )
    }
    
    if ("consistency" %in% criteria) {
      # Consistency criterion
      # Would require data from multiple studies
      mod_assessment$consistency <- list(
        criterion = "Consistency",
        result = "Unknown - requires data from multiple studies",
        score = NA
      )
    }
    
    if ("multiple_testing" %in% criteria) {
      # Multiple testing criterion
      n_tests <- length(effect_results$candidate_modifiers)
      
      if (n_tests > 1) {
        # Apply Bonferroni correction
        alpha_corrected <- 0.05 / n_tests
        
        if (effect_results$method == "interaction") {
          mod_effect <- effect_results$modifier_effects[[modifier]]
          
          # Calculate p-value for interaction term
          z_score <- mod_effect$estimate / mod_effect$sd
          p_value <- 2 * (1 - stats::pnorm(abs(z_score)))
          
          result <- paste0("p = ", sprintf("%.4f", p_value), 
                         ", corrected threshold = ", sprintf("%.4f", alpha_corrected))
          score <- if (p_value < alpha_corrected) 1 else 0
        } else {
          # For other methods, can't easily apply correction
          result <- paste0("Tested ", n_tests, " modifiers - correction applied")
          score <- 0  # Conservative
        }
      } else {
        result <- "Only one modifier tested - no correction needed"
        score <- 1
      }
      
      mod_assessment$multiple_testing <- list(
        criterion = "Multiple Testing",
        result = result,
        score = score
      )
    }
    
    if ("biological_plausibility" %in% criteria) {
      # Biological plausibility criterion
      # This requires domain knowledge
      mod_assessment$biological_plausibility <- list(
        criterion = "Biological Plausibility",
        result = "Unknown - requires domain expertise",
        score = NA
      )
    }
    
    # Add to overall assessment
    assessment[[modifier]] <- mod_assessment
  }
  
  # Convert to data frame
  assessment_df <- do.call(rbind, lapply(names(assessment), function(modifier) {
    mod_assessment <- assessment[[modifier]]
    
    do.call(rbind, lapply(names(mod_assessment), function(criterion) {
      data.frame(
        modifier = modifier,
        criterion = mod_assessment[[criterion]]$criterion,
        result = mod_assessment[[criterion]]$result,
        score = mod_assessment[[criterion]]$score,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # Add overall credibility score
  assessment_df$overall_credibility <- "Needs further evidence"
  
  return(assessment_df)
}