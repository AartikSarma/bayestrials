#' Check model diagnostics
#' 
#' @param model A BayesianModel object
#' @param diagnostics Vector of diagnostics to check
#' @param plot Logical indicating whether to generate diagnostic plots
#' @param quiet Logical indicating whether to suppress messages
#' @return List of diagnostic results
#' @export
check_diagnostics <- function(model, 
                             diagnostics = c("rhat", "ess", "divergences", "ppc"),
                             plot = TRUE,
                             quiet = FALSE) {
  
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot check diagnostics for a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check if brms and posterior are available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("posterior package is required but not available", call. = FALSE)
  }
  
  # Initialize diagnostics list
  diag_results <- list(
    warnings = character(),
    rhat = NULL,
    ess = NULL,
    divergences = NULL,
    ppc = NULL,
    plots = list()
  )
  
  # Extract model name
  model_name <- if (!is.null(model$model_spec) && !is.null(model$model_spec$model_name)) {
    model$model_spec$model_name
  } else {
    "unnamed"
  }
  
  # Check Rhat values
  if ("rhat" %in% diagnostics) {
    # Get rhat values
    rhat_values <- brms::rhat(model$model)
    diag_results$rhat <- rhat_values
    
    # Flag problematic parameters
    if (any(rhat_values > 1.01, na.rm = TRUE)) {
      problematic <- names(rhat_values)[rhat_values > 1.01]
      warning_msg <- paste0(
        "High Rhat detected for ", length(problematic), " parameters in model '", 
        model_name, "'. Max Rhat = ", round(max(rhat_values, na.rm = TRUE), 3)
      )
      
      if (!quiet) log_message(warning_msg, level = "warning")
      diag_results$warnings <- c(diag_results$warnings, warning_msg)
      
      # Add plot if requested
      if (plot && requireNamespace("bayesplot", quietly = TRUE)) {
        diag_results$plots$rhat <- bayesplot::mcmc_rhat(rhat_values) + 
          ggplot2::labs(title = "Rhat Values", subtitle = model_name)
      }
    } else {
      if (!quiet) log_message(paste0("All Rhat values look good for model '", model_name, "'"), level = "info")
    }
  }
  
  # Check effective sample size
  if ("ess" %in% diagnostics) {
    # Get ESS values
    ess_bulk <- brms::neff_ratio(model$model, robust = TRUE)
    diag_results$ess <- ess_bulk
    
    # Flag problematic parameters
    if (any(ess_bulk < 0.1, na.rm = TRUE)) {
      problematic <- names(ess_bulk)[ess_bulk < 0.1]
      warning_msg <- paste0(
        "Low ESS ratio (<10%) detected for ", length(problematic), " parameters in model '", 
        model_name, "'. Min ESS ratio = ", round(min(ess_bulk, na.rm = TRUE), 3)
      )
      
      if (!quiet) log_message(warning_msg, level = "warning")
      diag_results$warnings <- c(diag_results$warnings, warning_msg)
      
      # Add plot if requested
      if (plot && requireNamespace("bayesplot", quietly = TRUE)) {
        diag_results$plots$ess <- bayesplot::mcmc_neff(ess_bulk) + 
          ggplot2::labs(title = "Effective Sample Size Ratios", subtitle = model_name)
      }
    } else {
      if (!quiet) log_message(paste0("All ESS values look good for model '", model_name, "'"), level = "info")
    }
  }
  
  # Check for divergent transitions
  if ("divergences" %in% diagnostics) {
    # Extract divergence information
    if (inherits(model$model, "brmsfit")) {
      sampler_params <- brms::nuts_params(model$model)
      divergent <- which(sampler_params$Parameter == "divergent__")
      divergences <- sum(sampler_params$Value[divergent])
      diag_results$divergences <- divergences
      
      if (divergences > 0) {
        warning_msg <- paste0(
          divergences, " divergent transitions detected in model '", model_name, 
          "'. Consider increasing adapt_delta."
        )
        
        if (!quiet) log_message(warning_msg, level = "warning")
        diag_results$warnings <- c(diag_results$warnings, warning_msg)
        
        # Add plot if requested
        if (plot && requireNamespace("bayesplot", quietly = TRUE)) {
          # Extract all nuts parameters
          np <- brms::nuts_params(model$model)
          diag_results$plots$pairs <- bayesplot::mcmc_pairs(
            brms::posterior_samples(model$model)[, 1:min(4, ncol(brms::posterior_samples(model$model)))],
            np = np, pars = NULL,
            off_diag_args = list(size = 0.7)
          )
        }
      } else {
        if (!quiet) log_message(paste0("No divergent transitions detected in model '", model_name, "'"), level = "info")
      }
    }
  }
  
  # Posterior predictive checks
  if ("ppc" %in% diagnostics && plot) {
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      tryCatch({
        # Generate posterior predictive samples
        ppfit <- brms::pp_check(model$model, type = "dens_overlay", nsamples = 50)
        diag_results$plots$ppc_dens <- ppfit
        
        ppfit_scatter <- brms::pp_check(model$model, type = "scatter_avg")
        diag_results$plots$ppc_scatter <- ppfit_scatter
      }, error = function(e) {
        warning_msg <- paste0("Could not generate posterior predictive checks: ", e$message)
        if (!quiet) log_message(warning_msg, level = "warning")
        diag_results$warnings <- c(diag_results$warnings, warning_msg)
      })
    }
  }
  
  # Set overall diagnostics status
  if (length(diag_results$warnings) > 0) {
    diag_results$status <- "issues"
  } else {
    diag_results$status <- "good"
  }
  
  return(diag_results)
}

#' Plot model diagnostics
#' 
#' @param model A BayesianModel object
#' @param type Type of diagnostic plot to generate
#' @param parameters Vector of parameters to include in plots (NULL for all)
#' @return A ggplot object
#' @export
plot_diagnostics <- function(model, 
                            type = c("trace", "rhat", "ess", "pairs", "ppc", "intervals"),
                            parameters = NULL) {
  
  type <- match.arg(type)
  
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot plot diagnostics for a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check for required packages
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("bayesplot package is required but not available", call. = FALSE)
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Extract model name
  model_name <- if (!is.null(model$model_spec) && !is.null(model$model_spec$model_name)) {
    model$model_spec$model_name
  } else {
    "unnamed"
  }
  
  # Extract posterior draws
  posterior <- brms::posterior_samples(model$model)
  
  # Filter parameters if provided
  if (!is.null(parameters)) {
    # Convert to regex pattern for partial matching
    pattern <- paste0("^", parameters, collapse = "|")
    selected_params <- grep(pattern, colnames(posterior), value = TRUE)
    
    if (length(selected_params) == 0) {
      stop("No parameters match the provided patterns", call. = FALSE)
    }
    
    posterior <- posterior[, selected_params, drop = FALSE]
  }
  
  # Generate requested plot
  if (type == "trace") {
    # Trace plots
    if (ncol(posterior) > 6) {
      warning("Too many parameters for trace plot, showing first 6", call. = FALSE)
      posterior <- posterior[, 1:6, drop = FALSE]
    }
    
    plot <- bayesplot::mcmc_trace(
      as.matrix(posterior),
      facet_args = list(ncol = 2)
    ) + 
      ggplot2::labs(title = paste("Trace Plots -", model_name))
    
  } else if (type == "rhat") {
    # Rhat plots
    rhat_values <- brms::rhat(model$model)
    
    if (!is.null(parameters)) {
      pattern <- paste0("^", parameters, collapse = "|")
      selected_rhat <- rhat_values[grep(pattern, names(rhat_values))]
      
      if (length(selected_rhat) == 0) {
        stop("No parameters match the provided patterns for Rhat", call. = FALSE)
      }
      
      rhat_values <- selected_rhat
    }
    
    plot <- bayesplot::mcmc_rhat(rhat_values) + 
      ggplot2::labs(title = paste("Rhat Values -", model_name))
    
  } else if (type == "ess") {
    # ESS plots
    ess_values <- brms::neff_ratio(model$model)
    
    if (!is.null(parameters)) {
      pattern <- paste0("^", parameters, collapse = "|")
      selected_ess <- ess_values[grep(pattern, names(ess_values))]
      
      if (length(selected_ess) == 0) {
        stop("No parameters match the provided patterns for ESS", call. = FALSE)
      }
      
      ess_values <- selected_ess
    }
    
    plot <- bayesplot::mcmc_neff(ess_values) + 
      ggplot2::labs(title = paste("Effective Sample Size Ratios -", model_name))
    
  } else if (type == "pairs") {
    # Pairs plots
    if (ncol(posterior) > 6) {
      warning("Too many parameters for pairs plot, showing first 4", call. = FALSE)
      posterior <- posterior[, 1:4, drop = FALSE]
    }
    
    # Extract nuts parameters
    np <- brms::nuts_params(model$model)
    
    plot <- bayesplot::mcmc_pairs(
      as.matrix(posterior),
      np = np, pars = NULL,
      off_diag_args = list(size = 0.7)
    )
    
  } else if (type == "ppc") {
    # Posterior predictive checks
    plot <- brms::pp_check(model$model, type = "dens_overlay", nsamples = 50) +
      ggplot2::labs(title = paste("Posterior Predictive Check -", model_name))
    
  } else if (type == "intervals") {
    # Credible intervals
    if (ncol(posterior) > 20) {
      warning("Too many parameters for interval plot, showing first 20", call. = FALSE)
      posterior <- posterior[, 1:20, drop = FALSE]
    }
    
    plot <- bayesplot::mcmc_intervals(as.matrix(posterior)) +
      ggplot2::labs(title = paste("Posterior Intervals -", model_name))
  }
  
  return(plot)
}

#' Extract model fit statistics
#' 
#' @param model A BayesianModel object
#' @param type Type of fit statistic to compute (loo, waic)
#' @return Object containing model fit statistics
#' @export
model_fit_stats <- function(model, type = c("loo", "waic")) {
  type <- match.arg(type)
  
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot compute fit statistics for a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Compute requested statistic
  if (type == "loo") {
    loo_result <- tryCatch({
      brms::loo(model$model)
    }, error = function(e) {
      warning(paste("Error computing LOO:", e$message), call. = FALSE)
      return(NULL)
    })
    
    return(loo_result)
  } else if (type == "waic") {
    waic_result <- tryCatch({
      brms::waic(model$model)
    }, error = function(e) {
      warning(paste("Error computing WAIC:", e$message), call. = FALSE)
      return(NULL)
    })
    
    return(waic_result)
  }
}

#' Compute posterior predictive p-value
#' 
#' @param model A BayesianModel object
#' @param test_statistic Function computing test statistic from observed and replicated data
#' @param nsamples Number of posterior samples to use
#' @return Posterior predictive p-value
#' @export
posterior_predictive_pvalue <- function(model, test_statistic = NULL, nsamples = 1000) {
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot compute p-value for a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Define default test statistic if not provided
  if (is.null(test_statistic)) {
    # Default: use mean as test statistic
    test_statistic <- function(y_obs, y_rep) {
      mean(y_obs, na.rm = TRUE)
    }
  }
  
  # Get posterior predictive samples
  y_rep <- brms::posterior_predict(model$model, nsamples = nsamples)
  
  # Get observed data
  y_obs <- model$model$data[[model$model$formula$resp]]
  
  # Compute test statistic for observed data
  t_obs <- test_statistic(y_obs, NULL)
  
  # Compute test statistic for each replicated dataset
  t_rep <- apply(y_rep, 1, function(y) {
    test_statistic(NULL, y)
  })
  
  # Compute p-value as proportion of replications with test statistic >= observed
  p_value <- mean(t_rep >= t_obs, na.rm = TRUE)
  
  return(list(
    p_value = p_value,
    t_obs = t_obs,
    t_rep = t_rep
  ))
}

#' Parameter recovery test
#' 
#' @param model A BayesianModel object
#' @param true_values Named vector of true parameter values
#' @param parameters Vector of parameter names to assess (NULL for all)
#' @return Data frame with parameter recovery statistics
#' @export
parameter_recovery <- function(model, true_values, parameters = NULL) {
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot assess parameter recovery for a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check true_values
  if (!is.vector(true_values) || is.null(names(true_values))) {
    stop("true_values must be a named vector", call. = FALSE)
  }
  
  # Get posterior summary
  post_summary <- posterior_summary(model, parameters = parameters)
  
  # Find parameters with true values
  common_params <- intersect(post_summary$parameter, names(true_values))
  
  if (length(common_params) == 0) {
    stop("No parameters in common between model and true_values", call. = FALSE)
  }
  
  # Restrict to common parameters
  post_summary <- post_summary[post_summary$parameter %in% common_params, ]
  
  # Add true values
  post_summary$true_value <- true_values[post_summary$parameter]
  
  # Compute bias and relative error
  post_summary$bias <- post_summary$mean - post_summary$true_value
  post_summary$relative_error <- abs(post_summary$bias) / abs(post_summary$true_value)
  
  # Compute coverage
  post_summary$covered <- post_summary$true_value >= post_summary$q2.5 & 
                         post_summary$true_value <= post_summary$q97.5
  
  # Compute z-score
  post_summary$z_score <- (post_summary$true_value - post_summary$mean) / post_summary$sd
  
  return(post_summary)
}

#' Check residuals for a Bayesian model
#' 
#' @param model A BayesianModel object
#' @param type Type of residuals (default, pearson, response)
#' @param ndraws Number of posterior draws to use
#' @param plot Logical indicating whether to create a plot
#' @return Data frame of residuals or a plot if requested
#' @export
check_residuals <- function(model, type = c("default", "pearson", "response"), ndraws = 100, plot = TRUE) {
  type <- match.arg(type)
  
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Cannot check residuals for a model that hasn't been fitted", call. = FALSE)
  }
  
  # Check if brms is available
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not available", call. = FALSE)
  }
  
  # Extract response variable from model
  resp <- model$model$formula$resp
  
  # Extract observed response
  y_obs <- model$model$data[[resp]]
  
  # Compute fitted values
  fitted_values <- brms::fitted(model$model, scale = "response", summary = FALSE)
  
  # Select a subset of draws
  if (nrow(fitted_values) > ndraws) {
    set.seed(42)  # For reproducibility
    draw_indices <- sample(nrow(fitted_values), ndraws)
    fitted_subset <- fitted_values[draw_indices, ]
  } else {
    fitted_subset <- fitted_values
  }
  
  # Compute mean fitted values
  fitted_mean <- colMeans(fitted_subset)
  
  # Compute residuals
  if (type == "default") {
    # Default: observed - fitted
    residuals <- y_obs - fitted_mean
  } else if (type == "pearson") {
    # Pearson: (observed - fitted) / sqrt(var(fitted))
    fitted_var <- apply(fitted_subset, 2, var)
    residuals <- (y_obs - fitted_mean) / sqrt(fitted_var)
  } else if (type == "response") {
    # Response: observed
    residuals <- y_obs
  }
  
  # Create data frame
  residual_df <- data.frame(
    fitted = fitted_mean,
    residuals = residuals,
    observation = seq_along(y_obs)
  )
  
  # Create plot if requested
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    # Determine model name
    model_name <- if (!is.null(model$model_spec) && !is.null(model$model_spec$model_name)) {
      model$model_spec$model_name
    } else {
      "unnamed"
    }
    
    # Residuals vs fitted plot
    p1 <- ggplot2::ggplot(residual_df, ggplot2::aes(x = fitted, y = residuals)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "blue", se = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Residuals vs Fitted -", model_name),
        x = "Fitted values",
        y = "Residuals"
      )
    
    # Histogram of residuals
    p2 <- ggplot2::ggplot(residual_df, ggplot2::aes(x = residuals)) +
      ggplot2::geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Distribution of Residuals -", model_name),
        x = "Residuals",
        y = "Count"
      )
    
    # Combine plots if patchwork is available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      plot <- p1 + p2 + patchwork::plot_layout(ncol = 1)
    } else {
      # Return just the first plot
      plot <- p1
    }
    
    return(plot)
  } else {
    return(residual_df)
  }
}