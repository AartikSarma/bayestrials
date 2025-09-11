#' Extract Model Estimates with Credible Intervals
#'
#' Extracts parameter estimates from fitted Bayesian models with credible intervals.
#'
#' @param results BayesianModel object or list of results
#' @param parameters Vector of parameter names to extract (NULL for all)
#' @param credible_interval Width of credible interval (default: 0.95)
#' @param include_diagnostics Whether to include convergence diagnostics
#'
#' @return A tibble with parameter estimates and credible intervals
#' @export
extract_estimates <- function(results, 
                             parameters = NULL,
                             credible_interval = 0.95,
                             include_diagnostics = FALSE) {
  
  # Handle different input types
  if (inherits(results, "BayesianModel")) {
    results <- list(model1 = results)
  } else if (inherits(results, "list") && "results" %in% names(results)) {
    # Handle results from fit_bayesian_models
    results <- results$results
  }
  
  if (!is.list(results)) {
    cli::cli_abort("results must be a BayesianModel object or list of models")
  }
  
  # Initialize output
  all_estimates <- tibble::tibble()
  
  # Process each model
  for (model_name in names(results)) {
    model <- results[[model_name]]
    
    if (!inherits(model, "BayesianModel") || is.null(model$model)) {
      cli::cli_warn("Skipping {model_name}: not a fitted BayesianModel")
      next
    }
    
    # Get posterior summary
    post_summary <- brms::posterior_summary(model$model)
    
    # Filter parameters if specified
    if (!is.null(parameters)) {
      param_rows <- sapply(parameters, function(p) {
        which(grepl(p, rownames(post_summary)))
      })
      param_rows <- unique(unlist(param_rows))
      
      if (length(param_rows) > 0) {
        post_summary <- post_summary[param_rows, , drop = FALSE]
      } else {
        cli::cli_warn("No matching parameters found in {model_name}")
        next
      }
    }
    
    # Calculate credible interval bounds
    alpha <- 1 - credible_interval
    lower_col <- paste0("Q", round(alpha/2 * 100, 1))
    upper_col <- paste0("Q", round((1 - alpha/2) * 100, 1))
    
    # Create estimates tibble
    model_estimates <- tibble::tibble(
      model_name = model_name,
      parameter = rownames(post_summary),
      mean = post_summary[, "Estimate"],
      sd = post_summary[, "Est.Error"],
      !!lower_col := post_summary[, "Q2.5"],
      !!upper_col := post_summary[, "Q97.5"]
    )
    
    # Add model-specific information
    if (!is.null(model$model_spec)) {
      model_estimates$family <- model$model_spec$family
      model_estimates$prior_name <- "default"  # Could be enhanced
    }
    
    # Add diagnostics if requested
    if (include_diagnostics) {
      # Get Rhat values
      rhat_values <- brms::rhat(model$model)
      if (!is.null(rhat_values)) {
        rhat_df <- tibble::tibble(
          parameter = names(rhat_values),
          rhat = as.numeric(rhat_values)
        )
        model_estimates <- dplyr::left_join(model_estimates, rhat_df, by = "parameter")
      }
      
      # Get effective sample sizes
      neff_values <- brms::neff_ratio(model$model)
      if (!is.null(neff_values)) {
        neff_df <- tibble::tibble(
          parameter = names(neff_values),
          neff_ratio = as.numeric(neff_values)
        )
        model_estimates <- dplyr::left_join(model_estimates, neff_df, by = "parameter")
      }
    }
    
    all_estimates <- dplyr::bind_rows(all_estimates, model_estimates)
  }
  
  return(all_estimates)
}

#' Extract Convergence Diagnostics
#'
#' Extracts convergence diagnostics from fitted Bayesian models.
#'
#' @param results BayesianModel object or list of results
#' @param parameters Vector of parameter names to check (NULL for all)
#'
#' @return A tibble with convergence diagnostics
#' @export
extract_diagnostics <- function(results, parameters = NULL) {
  
  # Handle different input types
  if (inherits(results, "BayesianModel")) {
    results <- list(model1 = results)
  } else if (inherits(results, "list") && "results" %in% names(results)) {
    results <- results$results
  }
  
  if (!is.list(results)) {
    cli::cli_abort("results must be a BayesianModel object or list of models")
  }
  
  # Initialize output
  all_diagnostics <- tibble::tibble()
  
  # Process each model
  for (model_name in names(results)) {
    model <- results[[model_name]]
    
    # Initialize model diagnostics
    model_diag <- tibble::tibble(
      model_name = model_name,
      converged = FALSE,
      has_error = FALSE,
      error_message = NA_character_,
      rhat_max = NA_real_,
      rhat_mean = NA_real_,
      all_rhat_ok = FALSE,
      n_eff_min = NA_real_,
      n_eff_mean = NA_real_,
      all_neff_ok = FALSE,
      n_divergent = NA_integer_,
      no_divergences = FALSE,
      max_treedepth_hit = NA_integer_,
      bfmi_low = NA_integer_
    )
    
    # Check if model exists and is fitted
    if (!inherits(model, "BayesianModel") || is.null(model$model)) {
      model_diag$has_error <- TRUE
      model_diag$error_message <- "Model not fitted or invalid"
      all_diagnostics <- dplyr::bind_rows(all_diagnostics, model_diag)
      next
    }
    
    # Extract diagnostics
    tryCatch({
      # Rhat diagnostics
      rhat_values <- brms::rhat(model$model)
      if (!is.null(rhat_values)) {
        # Filter parameters if specified
        if (!is.null(parameters)) {
          param_indices <- sapply(parameters, function(p) {
            which(grepl(p, names(rhat_values)))
          })
          param_indices <- unique(unlist(param_indices))
          
          if (length(param_indices) > 0) {
            rhat_values <- rhat_values[param_indices]
          }
        }
        
        model_diag$rhat_max <- max(rhat_values, na.rm = TRUE)
        model_diag$rhat_mean <- mean(rhat_values, na.rm = TRUE)
        model_diag$all_rhat_ok <- all(rhat_values < 1.1, na.rm = TRUE)
      }
      
      # Effective sample size diagnostics
      neff_values <- brms::neff_ratio(model$model)
      if (!is.null(neff_values)) {
        # Filter parameters if specified
        if (!is.null(parameters)) {
          param_indices <- sapply(parameters, function(p) {
            which(grepl(p, names(neff_values)))
          })
          param_indices <- unique(unlist(param_indices))
          
          if (length(param_indices) > 0) {
            neff_values <- neff_values[param_indices]
          }
        }
        
        model_diag$n_eff_min <- min(neff_values, na.rm = TRUE)
        model_diag$n_eff_mean <- mean(neff_values, na.rm = TRUE)
        model_diag$all_neff_ok <- all(neff_values > 0.1, na.rm = TRUE)  # 10% of total samples
      }
      
      # HMC-specific diagnostics
      if (inherits(model$model$fit, "stanfit")) {
        sampler_params <- rstan::get_sampler_params(model$model$fit, inc_warmup = FALSE)
        
        # Count divergent transitions
        divergent <- do.call(rbind, sampler_params)[, 'divergent__']
        model_diag$n_divergent <- sum(divergent)
        model_diag$no_divergences <- model_diag$n_divergent == 0
        
        # Check for hitting max treedepth
        treedepth <- do.call(rbind, sampler_params)[, 'treedepth__']
        max_td <- model$model$fit@stan_args[[1]]$control$max_treedepth %||% 10
        model_diag$max_treedepth_hit <- sum(treedepth >= max_td)
        
        # Check for low BFMI (Bayesian Fraction of Missing Information)
        energy <- do.call(rbind, sampler_params)[, 'energy__']
        bfmi <- apply(energy, 2, function(x) {
          var_energy <- stats::var(x)
          mean_energy_var <- mean(diff(x)^2)
          bfmi <- var_energy / (var_energy + mean_energy_var)
          return(bfmi)
        })
        model_diag$bfmi_low <- sum(bfmi < 0.2)
      }
      
      # Overall convergence assessment
      model_diag$converged <- model_diag$all_rhat_ok && 
        model_diag$all_neff_ok && 
        model_diag$no_divergences &&
        (is.na(model_diag$max_treedepth_hit) || model_diag$max_treedepth_hit == 0) &&
        (is.na(model_diag$bfmi_low) || model_diag$bfmi_low == 0)
      
    }, error = function(e) {
      model_diag$has_error <- TRUE
      model_diag$error_message <- as.character(e$message)
    })
    
    all_diagnostics <- dplyr::bind_rows(all_diagnostics, model_diag)
  }
  
  return(all_diagnostics)
}

#' Create Diagnostic Summary
#'
#' Creates a comprehensive diagnostic summary for model results.
#'
#' @param results BayesianModel object or list of results
#' @param include_warnings Whether to include detailed warnings
#'
#' @return A list with diagnostic summaries and recommendations
#' @export
create_diagnostic_summary <- function(results, include_warnings = TRUE) {
  
  # Get diagnostics
  diagnostics <- extract_diagnostics(results)
  
  # Overall summary
  n_models <- nrow(diagnostics)
  n_converged <- sum(diagnostics$converged, na.rm = TRUE)
  n_errors <- sum(diagnostics$has_error, na.rm = TRUE)
  
  # Convergence summary
  convergence_rate <- n_converged / n_models * 100
  
  # Diagnostic issues
  rhat_issues <- sum(diagnostics$rhat_max > 1.1, na.rm = TRUE)
  neff_issues <- sum(diagnostics$n_eff_min < 0.1, na.rm = TRUE)
  divergence_issues <- sum(diagnostics$n_divergent > 0, na.rm = TRUE)
  treedepth_issues <- sum(diagnostics$max_treedepth_hit > 0, na.rm = TRUE)
  bfmi_issues <- sum(diagnostics$bfmi_low > 0, na.rm = TRUE)
  
  # Create summary
  summary <- list(
    overall = list(
      n_models = n_models,
      n_converged = n_converged,
      n_errors = n_errors,
      convergence_rate = convergence_rate
    ),
    issues = list(
      rhat_issues = rhat_issues,
      neff_issues = neff_issues,
      divergence_issues = divergence_issues,
      treedepth_issues = treedepth_issues,
      bfmi_issues = bfmi_issues
    ),
    recommendations = character()
  )
  
  # Generate recommendations
  if (convergence_rate < 80) {
    summary$recommendations <- c(summary$recommendations,
                                "Low convergence rate detected. Consider reviewing model specification.")
  }
  
  if (rhat_issues > 0) {
    summary$recommendations <- c(summary$recommendations,
                                paste("R-hat > 1.1 detected in", rhat_issues, "model(s). Consider increasing iterations."))
  }
  
  if (neff_issues > 0) {
    summary$recommendations <- c(summary$recommendations,
                                paste("Low effective sample size detected in", neff_issues, "model(s). Consider increasing iterations."))
  }
  
  if (divergence_issues > 0) {
    summary$recommendations <- c(summary$recommendations,
                                paste("Divergent transitions detected in", divergence_issues, "model(s). Consider increasing adapt_delta or reparameterizing model."))
  }
  
  if (treedepth_issues > 0) {
    summary$recommendations <- c(summary$recommendations,
                                paste("Maximum treedepth reached in", treedepth_issues, "model(s). Consider increasing max_treedepth."))
  }
  
  if (bfmi_issues > 0) {
    summary$recommendations <- c(summary$recommendations,
                                paste("Low BFMI detected in", bfmi_issues, "model(s). Model may be difficult to sample from."))
  }
  
  if (length(summary$recommendations) == 0) {
    summary$recommendations <- "All models appear to have converged successfully!"
  }
  
  # Add detailed diagnostics
  summary$detailed <- diagnostics
  
  # Add warnings if requested
  if (include_warnings) {
    warning_models <- diagnostics[diagnostics$has_error | !diagnostics$converged, ]
    if (nrow(warning_models) > 0) {
      summary$warnings <- warning_models
    }
  }
  
  class(summary) <- c("bayestrials_diagnostics", "list")
  return(summary)
}

#' Print Diagnostic Summary
#'
#' @param x A bayestrials_diagnostics object
#' @param ... Additional arguments (ignored)
#'
#' @return x invisibly
#' @export
print.bayestrials_diagnostics <- function(x, ...) {
  cli::cli_h1("Bayesian Model Diagnostics Summary")
  
  # Overall results
  cli::cli_h2("Overall Results")
  cli::cli_ul(c(
    "Total models: {x$overall$n_models}",
    "Converged models: {x$overall$n_converged}",
    "Models with errors: {x$overall$n_errors}",
    "Convergence rate: {round(x$overall$convergence_rate, 1)}%"
  ))
  
  # Issues
  if (any(unlist(x$issues) > 0)) {
    cli::cli_h2("Diagnostic Issues")
    if (x$issues$rhat_issues > 0) {
      cli::cli_alert_warning("R-hat issues: {x$issues$rhat_issues} model(s)")
    }
    if (x$issues$neff_issues > 0) {
      cli::cli_alert_warning("Effective sample size issues: {x$issues$neff_issues} model(s)")
    }
    if (x$issues$divergence_issues > 0) {
      cli::cli_alert_warning("Divergent transitions: {x$issues$divergence_issues} model(s)")
    }
    if (x$issues$treedepth_issues > 0) {
      cli::cli_alert_warning("Max treedepth reached: {x$issues$treedepth_issues} model(s)")
    }
    if (x$issues$bfmi_issues > 0) {
      cli::cli_alert_warning("Low BFMI: {x$issues$bfmi_issues} model(s)")
    }
  }
  
  # Recommendations
  cli::cli_h2("Recommendations")
  for (rec in x$recommendations) {
    if (grepl("successfully", rec)) {
      cli::cli_alert_success(rec)
    } else {
      cli::cli_alert_info(rec)
    }
  }
  
  invisible(x)
}