#' BayesianModel class
#' 
#' @description S3 class for storing fitted Bayesian models and metadata
#' 
#' @export
BayesianModel <- function(brms_model, model_spec, prior_spec = NULL, data_info = NULL) {
  structure(
    list(
      model = brms_model,
      model_spec = model_spec,
      prior_spec = prior_spec,
      data_info = data_info,
      diagnostics = NULL,
      sensitivity = NULL
    ),
    class = "BayesianModel"
  )
}

#' Print method for BayesianModel
#' 
#' @param x A BayesianModel object
#' @param ... Additional arguments passed to print
#' @return x invisibly
#' @export
print.BayesianModel <- function(x, ...) {
  cat("Bayesian Clinical Trial Model\n")
  cat("----------------------------\n")
  
  # Model specification
  if (!is.null(x$model_spec)) {
    cat("Model Specification:", ifelse(is.null(x$model_spec$model_name), 
                                      "[Unnamed]", x$model_spec$model_name), "\n")
    cat("  Formula: ", build_formula_string(x$model_spec), "\n", sep = "")
    cat("  Family: ", x$model_spec$family, "\n", sep = "")
    cat("  Link: ", x$model_spec$link, "\n", sep = "")
  } else {
    cat("Model Specification: [Not Available]\n")
  }
  
  # Model summary
  if (!is.null(x$model)) {
    if (inherits(x$model, "brmsfit")) {
      cat("\nModel Summary:\n")
      cat("  Algorithm: ", x$model$algorithm, "\n", sep = "")
      cat("  Samples: ", x$model$fit@sim$iter - x$model$fit@sim$warmup, "\n", sep = "")
      cat("  Chains: ", x$model$fit@sim$chains, "\n", sep = "")
      
      # Print treatment effect if available
      treatment_effect <- extract_treatment_effect(x)
      if (!is.null(treatment_effect)) {
        cat("\nTreatment Effect:\n")
        cat("  Estimate: ", sprintf("%.4f", treatment_effect$estimate), "\n", sep = "")
        cat("  95% CI: [", sprintf("%.4f, %.4f", treatment_effect$ci_lower, treatment_effect$ci_upper), "]\n", sep = "")
      }
    } else {
      cat("\nModel Object: [Non-brms Model]\n")
    }
  } else {
    cat("\nModel: [Not Fitted]\n")
  }
  
  # Diagnostics
  if (!is.null(x$diagnostics)) {
    cat("\nDiagnostics:\n")
    if (!is.null(x$diagnostics$warnings) && length(x$diagnostics$warnings) > 0) {
      cat("  Warnings: ", length(x$diagnostics$warnings), "\n", sep = "")
    }
    if (!is.null(x$diagnostics$rhat) && any(x$diagnostics$rhat > 1.01)) {
      cat("  Rhat issues: ", sum(x$diagnostics$rhat > 1.01), " parameters\n", sep = "")
    }
    if (!is.null(x$diagnostics$n_eff) && any(x$diagnostics$n_eff < 1000)) {
      cat("  ESS issues: ", sum(x$diagnostics$n_eff < 1000), " parameters\n", sep = "")
    }
    if (!is.null(x$diagnostics$divergences) && x$diagnostics$divergences > 0) {
      cat("  Divergences: ", x$diagnostics$divergences, "\n", sep = "")
    }
  }
  
  invisible(x)
}

#' Extract a summary of the treatment effect
#' 
#' @param model A BayesianModel object
#' @param parameter Name of the treatment parameter (default: "treatment")
#' @return Data frame with treatment effect estimate and credible interval
#' @keywords internal
extract_treatment_effect <- function(model, parameter = "treatment") {
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    return(NULL)
  }
  
  # Try to extract posterior summary
  if (requireNamespace("brms", quietly = TRUE) && inherits(model$model, "brmsfit")) {
    # Look for treatment parameter
    parameter_pattern <- paste0("^b_", parameter, "$")
    posterior_summary <- brms::posterior_summary(model$model)
    
    # Find parameter in posterior summary
    param_row <- which(grepl(parameter_pattern, rownames(posterior_summary)))
    
    if (length(param_row) > 0) {
      estimate <- posterior_summary[param_row, "Estimate"]
      ci_lower <- posterior_summary[param_row, "Q2.5"]
      ci_upper <- posterior_summary[param_row, "Q97.5"]
      
      return(data.frame(
        parameter = parameter,
        estimate = estimate,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(NULL)
}

#' Fit a Bayesian model using brms
#' 
#' @param model_spec A ModelSpecification object
#' @param data Data frame containing the data
#' @param chains Number of MCMC chains
#' @param cores Number of cores to use
#' @param iter Number of iterations per chain
#' @param seed Random seed for reproducibility
#' @param save_model Logical indicating whether to save model object
#' @param control List of control parameters for Stan
#' @return A BayesianModel object
#' @export
fit_model <- function(model_spec, 
                     data,
                     chains = 4,
                     cores = parallel::detectCores(),
                     iter = 2000,
                     seed = NULL,
                     save_model = TRUE,
                     control = list(adapt_delta = 0.9)) {
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Check model spec
  if (!inherits(model_spec, "ModelSpecification")) {
    stop("model_spec must be a ModelSpecification object", call. = FALSE)
  }
  
  # Validate model spec
  validate_model_spec(model_spec)
  
  # Set seed for reproducibility and ensure it's a single integer
  if (!is.null(seed)) {
    # Try to coerce to a single integer
    tryCatch({
      seed <- as.integer(seed)[1]
      set_seed(seed)
    }, error = function(e) {
      stop("seed must be coercible to a single integer value", call. = FALSE)
    }, warning = function(w) {
      # If there's a warning, coerce anyway but log it
      seed <- as.integer(seed)[1]
      log_message(paste("Warning when coercing seed:", w$message), level = "warning")
      set_seed(seed)
    })
  }
  
  # Get priors if available
  prior_spec <- attr(model_spec, "priors")
  brms_priors <- NULL
  
  if (!is.null(prior_spec)) {
    if (inherits(prior_spec, "PriorSpecification")) {
      brms_priors <- to_brms_prior(prior_spec, model_spec)
    } else {
      warning("Invalid prior_spec attribute, ignoring", call. = FALSE)
    }
  }
  
  # Build formula
  formula <- build_formula(model_spec)
  
  # Build family
  # Use the brms family constructor functions rather than building a string
  if (model_spec$family == "bernoulli") {
    family <- brms::bernoulli(link = model_spec$link)
  } else if (model_spec$family == "binomial") {
    family <- brms::binomial(link = model_spec$link)
  } else if (model_spec$family == "gaussian") {
    family <- brms::gaussian(link = model_spec$link)
  } else if (model_spec$family == "poisson") {
    family <- brms::poisson(link = model_spec$link)
  } else if (model_spec$family == "negbinomial") {
    family <- brms::negbinomial(link = model_spec$link)
  } else if (model_spec$family == "gamma") {
    family <- brms::Gamma(link = model_spec$link)
  } else if (model_spec$family == "weibull") {
    family <- brms::weibull(link = model_spec$link)
  } else {
    # For other families, try the string approach as fallback
    family <- paste0(model_spec$family, "(link = '", model_spec$link, "')")
  }
  
  # Store data info
  data_info <- list(
    n_obs = nrow(data),
    variables = names(data),
    timestamp = Sys.time()
  )
  
  # Define default control parameters if not provided
  if (is.null(control)) {
    # Use higher adapt_delta to help prevent divergent transitions
    control <- list(adapt_delta = 0.95, max_treedepth = 12)
  }
  
  # Fit model
  log_message("Fitting model...", level = "info")
  
  # Use a different approach to capture warnings without failing
  # Initialize warning messages container
  warning_messages <- character(0)
  
  # Check for build tools - try a basic Stan/brms operation that will fail if tools aren't installed
  build_tools_check <- try({
    # This is a minimally small Stan model that should compile quickly if tools are installed
    test_model <- "
    data {
      int<lower=0> N;
    }
    parameters {
      real mu;
    }
    model {
      mu ~ normal(0, 1);
    }
    "
    rstan::stan_model(model_code = test_model, model_name = "test_model", verbose = FALSE)
  }, silent = TRUE)
  
  # Check if we had a build tools error
  if (inherits(build_tools_check, "try-error")) {
    # General pattern matching for missing build tools
    has_build_error <- grepl("tools for compilation", build_tools_check[1]) || 
                       grepl("C++ compiler", build_tools_check[1]) ||
                       grepl("cannot compile", build_tools_check[1]) ||
                       grepl("'cmath' file not found", build_tools_check[1]) ||
                       grepl("make: \\*\\*\\* .*Error", build_tools_check[1])
    
    if (has_build_error) {
      # Detect OS for more specific instructions
      is_mac <- Sys.info()["sysname"] == "Darwin"
      is_windows <- .Platform$OS.type == "windows"
      
      if (is_mac) {
        error_msg <- paste(
          "macOS: Missing C++ compiler tools needed by Stan/brms.",
          "To install required tools, open Terminal and run: xcode-select --install",
          "Then follow the prompts to install the Xcode Command Line Tools.",
          "After installation completes, restart R/RStudio and try again.",
          "For detailed help, see: https://mac.r-project.org/tools/"
        )
      } else if (is_windows) {
        r_version <- paste0(R.Version()$major, ".", substr(R.Version()$minor, 1, 1))
        error_msg <- paste(
          "Windows: Missing Rtools needed by Stan/brms.",
          paste0("For R ", r_version, ", download and install Rtools from:"),
          "https://cran.r-project.org/bin/windows/Rtools/",
          "Ensure you check the 'Add Rtools to system PATH' option during installation.",
          "After installation completes, restart R/RStudio and try again."
        )
      } else {
        error_msg <- paste(
          "Linux: Missing C++ compiler tools needed by Stan/brms.",
          "For Ubuntu/Debian, run: sudo apt-get install build-essential",
          "For Fedora/RHEL, run: sudo dnf install gcc-c++ make",
          "After installation completes, restart R/RStudio and try again.",
          "See https://mc-stan.org/docs/stan-users-guide/prereqs.html for more details."
        )
      }
      
      log_message(error_msg, level = "error")
      return(NULL)
    } else {
      # Other non-build related error, just log it
      log_message(paste("Error in test compilation:", build_tools_check[1]), level = "error")
      # But don't exit, maybe brm will still work
    }
  }
  
  # If we passed the build tools check, continue with model fitting
  # Set up warning handler
  withCallingHandlers({
    model <- brms::brm(
      formula = formula,
      data = data,
      family = family,
      prior = brms_priors,
      chains = chains,
      cores = cores,
      iter = iter,
      # seed parameter removed
      control = control,
      silent = 2
    )
  }, warning = function(w) {
    # Capture the warning
    warning_messages <<- c(warning_messages, w$message)
    log_message(paste("Warning in model fitting:", w$message), level = "warning")
    # Keep going, don't stop on warnings
    invokeRestart("muffleWarning")
  }, error = function(e) {
    # Check specifically for cmath errors that might appear during the actual model fitting
    if (grepl("'cmath' file not found", e$message) || 
        grepl("make: \\*\\*\\* .*Error", e$message)) {
      log_message(paste(
        "macOS compiler error: Missing Xcode Command Line Tools.",
        "Open Terminal and run: xcode-select --install",
        "Then restart R/RStudio and try again."
      ), level = "error")
    } else {
      log_message(paste("Error fitting model:", e$message), level = "error")
    }
    return(NULL)
  })
  
  # If no model was created but we have warnings, set model to NULL with message
  if (!exists("model") && length(warning_messages) > 0) {
    log_message("Model fitting failed due to warnings", level = "error")
    model <- NULL
  }
  
  # Check if model fitting was successful
  if (is.null(model)) {
    log_message("Model fitting failed, returning NULL", level = "error")
    return(NULL)
  }
  
  # Create BayesianModel object
  bayes_model <- BayesianModel(
    brms_model = model,
    model_spec = model_spec,
    prior_spec = prior_spec,
    data_info = data_info
  )
  
  # Check diagnostics
  bayes_model$diagnostics <- check_diagnostics(bayes_model, quiet = TRUE)
  
  return(bayes_model)
}

#' Update an existing model
#' 
#' @param model A BayesianModel object
#' @param formula New formula (optional)
#' @param priors New prior specifications (optional)
#' @param refit Logical indicating whether to refit the model
#' @param ... Additional arguments passed to brms::update
#' @return Updated BayesianModel object
#' @export
update_model <- function(model, formula = NULL, priors = NULL, refit = TRUE, ...) {
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot update a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Update model specification if formula provided
  if (!is.null(formula)) {
    if (inherits(formula, "ModelSpecification")) {
      model$model_spec <- formula
    } else {
      # Try to update ModelSpecification
      if (!is.null(model$model_spec)) {
        model$model_spec$predictors <- formula
      } else {
        warning("Cannot update model specification, no ModelSpecification object found", call. = FALSE)
      }
    }
  }
  
  # Update priors if provided
  if (!is.null(priors)) {
    if (inherits(priors, "PriorSpecification")) {
      model$prior_spec <- priors
    } else {
      warning("priors must be a PriorSpecification object", call. = FALSE)
    }
  }
  
  # If refit requested, update the model
  if (refit) {
    # Prepare formula
    update_formula <- NULL
    if (!is.null(formula)) {
      if (inherits(formula, "ModelSpecification")) {
        update_formula <- build_formula(formula)
      } else {
        update_formula <- formula
      }
    }
    
    # Prepare priors
    brms_priors <- NULL
    if (!is.null(priors) && inherits(priors, "PriorSpecification")) {
      brms_priors <- to_brms_prior(priors, update_formula)
    }
    
    # Update the model
    updated_brms_model <- brms::update(
      model$model, 
      formula = update_formula,
      prior = brms_priors,
      ...
    )
    
    # Update the BayesianModel object
    model$model <- updated_brms_model
    model$diagnostics <- check_diagnostics(model, quiet = TRUE)
  }
  
  return(model)
}

#' Fit multiple models
#' 
#' @param model_specs List of ModelSpecification objects
#' @param data Data frame containing the data
#' @param parallel Logical indicating whether to use parallel processing
#' @param ... Additional arguments passed to fit_model
#' @return List of BayesianModel objects
#' @export
fit_models <- function(model_specs, data, parallel = TRUE, ...) {
  if (!is.list(model_specs)) {
    stop("model_specs must be a list of ModelSpecification objects", call. = FALSE)
  }
  
  # Check that all elements are ModelSpecification objects
  if (!all(sapply(model_specs, inherits, "ModelSpecification"))) {
    stop("All elements of model_specs must be ModelSpecification objects", call. = FALSE)
  }
  
  # Use names from the list or create them
  if (is.null(names(model_specs))) {
    names(model_specs) <- paste0("model", seq_along(model_specs))
  }
  
  # Define function to fit a single model
  fit_one_model <- function(spec, data, ...) {
    model_name <- spec$model_name
    if (is.null(model_name) || model_name == "") {
      model_name <- names(model_specs)[match(spec, model_specs)]
    }
    
    log_message(paste("Fitting model:", model_name), level = "info")
    
    # Fit the model
    model <- fit_model(spec, data, ...)
    
    return(model)
  }
  
  # Fit models (in parallel if requested)
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    log_message("Using parallel processing to fit models", level = "info")
    
    # Determine number of cores to use
    cores <- min(length(model_specs), parallel::detectCores() - 1)
    cores <- max(cores, 1)  # Ensure at least 1 core
    
    # Create cluster
    cl <- parallel::makeCluster(cores)
    
    # Export needed functions and data
    parallel::clusterExport(cl, varlist = c("fit_model", "validate_model_spec", 
                                          "build_formula", "build_formula_string", 
                                          "log_message", "set_seed", 
                                          "check_diagnostics", "BayesianModel"), 
                           envir = environment())
    
    # Ensure required packages are loaded on each worker
    parallel::clusterEvalQ(cl, {
      library(brms)
      library(methods)
    })
    
    # Fit models in parallel
    models <- parallel::parLapply(cl, model_specs, fit_one_model, data = data, ...)
    
    # Stop cluster
    parallel::stopCluster(cl)
  } else {
    # Fit models sequentially
    models <- lapply(model_specs, fit_one_model, data = data, ...)
  }
  
  # Name the models
  names(models) <- names(model_specs)
  
  return(models)
}

#' Extract posterior samples from a model
#' 
#' @param model A BayesianModel object
#' @param variable Name of the variable to extract
#' @param transform Function to transform samples (optional)
#' @return Data frame of posterior samples
#' @export
extract_posterior <- function(model, variable, transform = NULL) {
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot extract posterior from a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Extract posterior draws
  if (inherits(model$model, "brmsfit")) {
    posterior_draws <- brms::posterior_samples(model$model)
    
    # Select variables matching pattern
    if (is.null(variable)) {
      # Return all variables
      selected_draws <- posterior_draws
    } else {
      # Get columns matching variable pattern
      var_pattern <- paste0("^", variable)
      var_cols <- grep(var_pattern, colnames(posterior_draws))
      
      if (length(var_cols) == 0) {
        stop(paste("No variables matching pattern:", variable), call. = FALSE)
      }
      
      selected_draws <- posterior_draws[, var_cols, drop = FALSE]
    }
    
    # Apply transformation if provided
    if (!is.null(transform) && is.function(transform)) {
      selected_draws <- as.data.frame(lapply(selected_draws, transform))
    }
    
    return(selected_draws)
  } else {
    stop("Unsupported model type for posterior extraction", call. = FALSE)
  }
}

#' Summarize posterior samples
#' 
#' @param model A BayesianModel object
#' @param variables Vector of variable names to summarize (NULL for all)
#' @param probs Numeric vector of quantiles to compute
#' @param transform Function to transform samples (optional)
#' @return Data frame of posterior summaries
#' @export
posterior_summary <- function(model, variables = NULL, probs = c(0.025, 0.5, 0.975), transform = NULL) {
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  # Extract posterior draws for requested variables
  posterior_draws <- extract_posterior(model, variables, transform)
  
  # Initialize results list
  summary_list <- list()
  
  # Compute summary statistics for each parameter
  for (param in colnames(posterior_draws)) {
    samples <- posterior_draws[[param]]
    
    # Compute statistics
    mean_val <- mean(samples)
    sd_val <- sd(samples)
    quantiles <- quantile(samples, probs = probs)
    
    # Create summary row
    summary_row <- c(
      parameter = param,
      mean = mean_val,
      sd = sd_val,
      setNames(as.numeric(quantiles), paste0("q", probs * 100))
    )
    
    summary_list[[param]] <- summary_row
  }
  
  # Convert to data frame
  summary_df <- as.data.frame(do.call(rbind, summary_list), stringsAsFactors = FALSE)
  
  # Convert character columns to numeric where appropriate
  numeric_cols <- c("mean", "sd", paste0("q", probs * 100))
  summary_df[, numeric_cols] <- lapply(summary_df[, numeric_cols], as.numeric)
  
  return(summary_df)
}

#' Compare models using information criteria
#' 
#' @param ... BayesianModel objects to compare
#' @param criterion Criterion for comparison (loo, waic, kfold)
#' @param seed Random seed for cross-validation
#' @return Comparison results from brms
#' @export
compare_models <- function(..., criterion = c("loo", "waic", "kfold"), seed = NULL) {
  criterion <- match.arg(criterion)
  
  # Extract models
  models <- list(...)
  
  # Check if models are BayesianModel objects and extract brms models
  brms_models <- list()
  model_names <- character(length(models))
  
  for (i in seq_along(models)) {
    model <- models[[i]]
    
    if (inherits(model, "BayesianModel")) {
      if (is.null(model$model)) {
        stop(paste("Model", i, "has not been fitted"), call. = FALSE)
      }
      
      brms_models[[i]] <- model$model
      
      # Get model name
      if (!is.null(model$model_spec) && !is.null(model$model_spec$model_name)) {
        model_names[i] <- model$model_spec$model_name
      } else {
        model_names[i] <- paste0("model", i)
      }
    } else if (inherits(model, "brmsfit")) {
      brms_models[[i]] <- model
      model_names[i] <- paste0("model", i)
    } else {
      stop(paste("Model", i, "must be a BayesianModel or brmsfit object"), call. = FALSE)
    }
  }
  
  # Set names
  names(brms_models) <- model_names
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set_seed(seed)
  }
  
  # Compare models using the specified criterion
  if (criterion == "loo") {
    # Check if brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("brms package is required but not available", call. = FALSE)
    }
    
    # Calculate LOO for each model
    loo_list <- lapply(brms_models, brms::loo)
    
    # Compare LOO values
    comparison <- brms::loo_compare(loo_list)
    
    # Add model names
    rownames(comparison) <- model_names
    
    return(comparison)
  } else if (criterion == "waic") {
    # Check if brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("brms package is required but not available", call. = FALSE)
    }
    
    # Calculate WAIC for each model
    waic_list <- lapply(brms_models, brms::waic)
    
    # Compare WAIC values
    comparison <- brms::loo_compare(waic_list)
    
    # Add model names
    rownames(comparison) <- model_names
    
    return(comparison)
  } else if (criterion == "kfold") {
    # Check if brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("brms package is required but not available", call. = FALSE)
    }
    
    # Calculate K-fold CV for each model
    kfold_list <- lapply(brms_models, brms::kfold_cv, folds = 10)
    
    # Compare K-fold values
    comparison <- brms::loo_compare(kfold_list)
    
    # Add model names
    rownames(comparison) <- model_names
    
    return(comparison)
  }
}

#' Extract trial-specific effects from a hierarchical model
#' 
#' @param model A BayesianModel object with hierarchical structure
#' @param parameter Name of the parameter to extract
#' @param transform Function to transform samples (optional)
#' @return Data frame of trial-specific effects
#' @export
extract_trial_effects <- function(model, parameter, transform = NULL) {
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot extract effects from a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Extract posterior draws
  if (inherits(model$model, "brmsfit")) {
    # Get overall effect
    overall_pattern <- paste0("^b_", parameter, "$")
    overall_cols <- grep(overall_pattern, colnames(brms::posterior_samples(model$model)))
    
    if (length(overall_cols) == 0) {
      warning("No overall effect found for parameter:", parameter, call. = FALSE)
    }
    
    # Get trial-specific effects (assuming r_trial[*,parameter] naming convention)
    trial_pattern <- paste0("^r_trial\\[.+,", parameter, "\\]$")
    trial_cols <- grep(trial_pattern, colnames(brms::posterior_samples(model$model)))
    
    if (length(trial_cols) == 0) {
      stop(paste("No trial-specific effects found for parameter:", parameter), call. = FALSE)
    }
    
    # Extract posterior samples
    posterior_draws <- brms::posterior_samples(model$model)
    
    # Extract trial names from column names
    trial_names <- gsub(paste0("^r_trial\\[(.+),", parameter, "\\]$"), "\\1", 
                       colnames(posterior_draws)[trial_cols])
    
    # Combine overall effect with trial-specific deviations
    trial_effects <- data.frame(posterior_draws[, trial_cols, drop = FALSE])
    
    # Rename columns to trial names
    colnames(trial_effects) <- trial_names
    
    # Add overall effect if available
    if (length(overall_cols) > 0) {
      overall_effect <- posterior_draws[, overall_cols[1], drop = FALSE]
      colnames(overall_effect) <- "overall"
      
      # For each trial, add overall effect to the trial-specific deviation
      for (trial in trial_names) {
        trial_effects[[trial]] <- trial_effects[[trial]] + overall_effect[[1]]
      }
      
      # Include overall effect in the results
      trial_effects <- cbind(overall_effect, trial_effects)
    }
    
    # Apply transformation if provided
    if (!is.null(transform) && is.function(transform)) {
      trial_effects <- as.data.frame(lapply(trial_effects, transform))
    }
    
    # Calculate summary statistics
    summary_list <- list()
    
    for (trial in colnames(trial_effects)) {
      samples <- trial_effects[[trial]]
      
      # Compute statistics
      mean_val <- mean(samples)
      sd_val <- sd(samples)
      ci_lower <- quantile(samples, 0.025)
      ci_upper <- quantile(samples, 0.975)
      
      summary_list[[trial]] <- data.frame(
        trial = trial,
        parameter = parameter,
        estimate = mean_val,
        sd = sd_val,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        stringsAsFactors = FALSE
      )
    }
    
    # Combine summaries
    summary_df <- do.call(rbind, summary_list)
    rownames(summary_df) <- NULL
    
    # Add posterior samples as attribute
    attr(summary_df, "posterior_samples") <- trial_effects
    
    return(summary_df)
  } else {
    stop("Unsupported model type for effect extraction", call. = FALSE)
  }
}