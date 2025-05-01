#' Run sensitivity analyses
#' 
#' @param base_model A BayesianModel object representing the base model
#' @param data Data frame containing the data
#' @param sensitivity_specs List of model specifications for sensitivity analyses
#' @param type Type of sensitivity analysis (prior, model, missing_data)
#' @param parallel Logical indicating whether to use parallel processing
#' @param ... Additional arguments passed to fit_model
#' @return List of sensitivity analysis results
#' @export
run_sensitivity_analyses <- function(base_model,
                                   data,
                                   sensitivity_specs,
                                   type = c("prior", "model", "missing_data"),
                                   parallel = TRUE,
                                   ...) {
  
  type <- match.arg(type)
  
  if (!inherits(base_model, "BayesianModel")) {
    stop("base_model must be a BayesianModel object", call. = FALSE)
  }
  
  if (!is.list(sensitivity_specs)) {
    stop("sensitivity_specs must be a list", call. = FALSE)
  }
  
  # Set up results structure
  results <- list(
    base_model = base_model,
    sensitivity_models = list(),
    comparisons = list(),
    type = type,
    summary = NULL
  )
  
  # Fit sensitivity models
  if (type == "prior") {
    # Prior sensitivity
    log_message("Running prior sensitivity analyses", level = "info")
    
    # Check that sensitivity_specs are all ModelSpecification objects
    if (!all(sapply(sensitivity_specs, inherits, "ModelSpecification"))) {
      stop("For prior sensitivity, all sensitivity_specs must be ModelSpecification objects with different priors", 
           call. = FALSE)
    }
    
    # Fit all models
    sensitivity_models <- fit_models(sensitivity_specs, data, parallel = parallel, ...)
    results$sensitivity_models <- sensitivity_models
    
    # Compare posterior distributions for key parameters
    if (!is.null(base_model$model) && !is.null(sensitivity_models) && length(sensitivity_models) > 0) {
      # Get treatment parameter if available
      treatment_param <- NULL
      if (!is.null(base_model$model_spec)) {
        pred_terms <- strsplit(base_model$model_spec$predictors, "\\+|\\*")[[1]]
        pred_terms <- trimws(pred_terms)
        if (any(grepl("treatment", pred_terms, ignore.case = TRUE))) {
          treatment_param <- "b_treatment"
        } else if (any(grepl("arm", pred_terms, ignore.case = TRUE))) {
          treatment_param <- "b_arm"
        } else if (any(grepl("trt", pred_terms, ignore.case = TRUE))) {
          treatment_param <- "b_trt"
        }
      }
      
      if (is.null(treatment_param)) {
        # Try to find a key parameter
        posterior_names <- colnames(brms::posterior_samples(base_model$model))
        if (any(grepl("^b_", posterior_names))) {
          treatment_param <- posterior_names[grep("^b_", posterior_names)[1]]
        }
      }
      
      if (!is.null(treatment_param)) {
        # Compare posteriors for this parameter
        results$comparisons$treatment <- compare_posteriors(
          c(list(base = base_model), sensitivity_models), 
          parameter = treatment_param
        )
      }
    }
    
    # Create summary
    results$summary <- summarize_sensitivity(results)
    
  } else if (type == "model") {
    # Model structure sensitivity
    log_message("Running model structure sensitivity analyses", level = "info")
    
    # Check that sensitivity_specs are all ModelSpecification objects
    if (!all(sapply(sensitivity_specs, inherits, "ModelSpecification"))) {
      stop("For model sensitivity, all sensitivity_specs must be ModelSpecification objects with different model structures", 
           call. = FALSE)
    }
    
    # Fit all models
    sensitivity_models <- fit_models(sensitivity_specs, data, parallel = parallel, ...)
    results$sensitivity_models <- sensitivity_models
    
    # Compare models using LOO or WAIC
    if (!is.null(base_model$model) && !is.null(sensitivity_models) && length(sensitivity_models) > 0) {
      # Calculate LOO for base model
      base_loo <- brms::loo(base_model$model)
      
      # Calculate LOO for sensitivity models
      sensitivity_loos <- lapply(sensitivity_models, function(model) {
        if (!is.null(model$model)) {
          brms::loo(model$model)
        } else {
          NULL
        }
      })
      
      # Remove NULL results
      sensitivity_loos <- sensitivity_loos[!sapply(sensitivity_loos, is.null)]
      
      if (length(sensitivity_loos) > 0) {
        # Combine all LOOs
        all_loos <- c(list(base = base_loo), sensitivity_loos)
        
        # Compare models
        results$comparisons$loo <- brms::loo_compare(all_loos)
      }
    }
    
    # Create summary
    results$summary <- summarize_sensitivity(results)
    
  } else if (type == "missing_data") {
    # Missing data sensitivity
    log_message("Running missing data sensitivity analyses", level = "info")
    
    # For missing data sensitivity, sensitivity_specs should be a list of modified datasets
    # or a list of functions that modify the data
    if (all(sapply(sensitivity_specs, is.data.frame))) {
      # List of modified datasets
      sensitivity_data <- sensitivity_specs
    } else if (all(sapply(sensitivity_specs, is.function))) {
      # List of functions to modify data
      sensitivity_data <- lapply(sensitivity_specs, function(func) func(data))
    } else {
      stop("For missing data sensitivity, sensitivity_specs must be either data frames or functions that return data frames", 
           call. = FALSE)
    }
    
    # Get base model specification
    base_spec <- base_model$model_spec
    if (is.null(base_spec)) {
      stop("Base model must have a model_spec to run missing data sensitivity", call. = FALSE)
    }
    
    # Fit models with different data
    sensitivity_models <- list()
    
    for (i in seq_along(sensitivity_data)) {
      name <- names(sensitivity_data)[i]
      if (is.null(name) || name == "") name <- paste0("missing_data_", i)
      
      # Fit model with this dataset
      log_message(paste("Fitting model for", name), level = "info")
      
      sensitivity_models[[name]] <- fit_model(base_spec, sensitivity_data[[i]], ...)
    }
    
    results$sensitivity_models <- sensitivity_models
    
    # Compare treatment effects
    if (!is.null(base_model$model) && !is.null(sensitivity_models) && length(sensitivity_models) > 0) {
      # Get treatment parameter if available
      treatment_param <- NULL
      if (!is.null(base_model$model_spec)) {
        pred_terms <- strsplit(base_model$model_spec$predictors, "\\+|\\*")[[1]]
        pred_terms <- trimws(pred_terms)
        if (any(grepl("treatment", pred_terms, ignore.case = TRUE))) {
          treatment_param <- "b_treatment"
        } else if (any(grepl("arm", pred_terms, ignore.case = TRUE))) {
          treatment_param <- "b_arm"
        } else if (any(grepl("trt", pred_terms, ignore.case = TRUE))) {
          treatment_param <- "b_trt"
        }
      }
      
      if (is.null(treatment_param)) {
        # Try to find a key parameter
        posterior_names <- colnames(brms::posterior_samples(base_model$model))
        if (any(grepl("^b_", posterior_names))) {
          treatment_param <- posterior_names[grep("^b_", posterior_names)[1]]
        }
      }
      
      if (!is.null(treatment_param)) {
        # Compare posteriors for this parameter
        results$comparisons$treatment <- compare_posteriors(
          c(list(base = base_model), sensitivity_models), 
          parameter = treatment_param
        )
      }
    }
    
    # Create summary
    results$summary <- summarize_sensitivity(results)
  }
  
  # Add class for method dispatch
  class(results) <- c("sensitivity_results", class(results))
  
  return(results)
}

#' Compare posterior distributions across models
#' 
#' @param models List of BayesianModel objects
#' @param parameter Name of the parameter to compare
#' @param transform Function to transform samples (optional)
#' @return Data frame with posterior comparison
#' @keywords internal
compare_posteriors <- function(models, parameter, transform = NULL) {
  if (!is.list(models) || length(models) < 2) {
    stop("models must be a list with at least two models", call. = FALSE)
  }
  
  # Check that all models are fitted BayesianModel objects
  valid_models <- sapply(models, function(m) {
    inherits(m, "BayesianModel") && !is.null(m$model)
  })
  
  if (!all(valid_models)) {
    warning("Some models are not valid BayesianModel objects or haven't been fitted",
           call. = FALSE)
    models <- models[valid_models]
  }
  
  if (length(models) < 2) {
    stop("Need at least two valid models for comparison", call. = FALSE)
  }
  
  # Extract posteriors for parameter
  posteriors <- list()
  
  for (name in names(models)) {
    model <- models[[name]]
    
    # Extract posterior samples
    post_samples <- extract_posterior(model, parameter, transform)
    
    # Check if parameter exists
    if (ncol(post_samples) > 0) {
      # Get first matching column
      param_col <- grep(parameter, colnames(post_samples), value = TRUE)[1]
      posteriors[[name]] <- post_samples[[param_col]]
    } else {
      warning(paste("Parameter", parameter, "not found in model", name), call. = FALSE)
    }
  }
  
  # Calculate posterior summaries
  summary_list <- list()
  
  for (name in names(posteriors)) {
    samples <- posteriors[[name]]
    
    # Calculate summary statistics
    summary_list[[name]] <- data.frame(
      model = name,
      parameter = parameter,
      mean = mean(samples),
      sd = sd(samples),
      q2.5 = quantile(samples, 0.025),
      q25 = quantile(samples, 0.25),
      q50 = quantile(samples, 0.5),
      q75 = quantile(samples, 0.75),
      q97.5 = quantile(samples, 0.975),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all summaries
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  # Calculate differences from base model
  base_mean <- summary_df$mean[summary_df$model == "base"]
  
  if (length(base_mean) > 0) {
    summary_df$diff_from_base <- summary_df$mean - base_mean
    summary_df$pct_diff_from_base <- (summary_df$mean - base_mean) / abs(base_mean) * 100
  }
  
  # Store raw posteriors as attribute
  attr(summary_df, "posteriors") <- posteriors
  
  return(summary_df)
}

#' Update model with different missing data mechanism
#' 
#' @param model_spec A ModelSpecification object
#' @param data Data frame containing the data
#' @param missing_mechanism Type of missing data mechanism (mar, mnar)
#' @param missing_vars Variables with missing values to handle
#' @param method Method for handling missing data (imputation, ...)
#' @return Updated data frame or model specification
#' @export
update_model_missing <- function(model_spec, 
                                data, 
                                missing_mechanism = c("mar", "mnar"),
                                missing_vars = NULL,
                                method = c("imputation", "model_based", "complete_case")) {
  
  missing_mechanism <- match.arg(missing_mechanism)
  method <- match.arg(method)
  
  if (!inherits(model_spec, "ModelSpecification")) {
    stop("model_spec must be a ModelSpecification object", call. = FALSE)
  }
  
  # Identify variables with missing values if not specified
  if (is.null(missing_vars)) {
    missing_counts <- sapply(data, function(x) sum(is.na(x)))
    missing_vars <- names(missing_counts)[missing_counts > 0]
  }
  
  if (length(missing_vars) == 0) {
    warning("No missing values found in data", call. = FALSE)
    return(list(model_spec = model_spec, data = data))
  }
  
  # Handle missing data based on method
  if (method == "imputation") {
    # Simple imputation
    imputed_data <- data
    
    for (var in missing_vars) {
      if (var %in% names(data)) {
        if (is.numeric(data[[var]])) {
          # Impute with mean for numeric variables
          imputed_data[[var]][is.na(data[[var]])] <- mean(data[[var]], na.rm = TRUE)
        } else if (is.factor(data[[var]]) || is.character(data[[var]])) {
          # Impute with mode for categorical variables
          mode_val <- names(sort(table(data[[var]]), decreasing = TRUE)[1])
          imputed_data[[var]][is.na(data[[var]])] <- mode_val
        }
      }
    }
    
    return(list(model_spec = model_spec, data = imputed_data))
    
  } else if (method == "model_based") {
    # Model-based imputation using brms
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("brms package is required but not available", call. = FALSE)
    }
    
    # For MAR, we can use brms's built-in missing data handling
    if (missing_mechanism == "mar") {
      # Just return the original specs, as brms will handle MAR appropriately
      return(list(model_spec = model_spec, data = data))
    } else if (missing_mechanism == "mnar") {
      # For MNAR, we need to model the missing data mechanism explicitly
      # This is a complex topic and here we just implement a simple version
      # where we add indicators for missingness to the model
      
      # Create missingness indicators
      for (var in missing_vars) {
        if (var %in% names(data)) {
          indicator_name <- paste0("missing_", var)
          data[[indicator_name]] <- as.numeric(is.na(data[[var]]))
        }
      }
      
      # Update model to include missingness indicators
      new_predictors <- model_spec$predictors
      
      for (var in missing_vars) {
        indicator_name <- paste0("missing_", var)
        if (indicator_name %in% names(data)) {
          new_predictors <- paste0(new_predictors, " + ", indicator_name)
        }
      }
      
      new_model_spec <- model_spec
      new_model_spec$predictors <- new_predictors
      
      # Set model name to indicate MNAR sensitivity
      if (!is.null(model_spec$model_name) && model_spec$model_name != "") {
        new_model_spec$model_name <- paste0(model_spec$model_name, "_mnar")
      } else {
        new_model_spec$model_name <- "mnar_sensitivity"
      }
      
      return(list(model_spec = new_model_spec, data = data))
    }
  } else if (method == "complete_case") {
    # Complete case analysis
    complete_data <- stats::na.omit(data)
    
    # Warn if too many cases are dropped
    dropped_pct <- (nrow(data) - nrow(complete_data)) / nrow(data) * 100
    if (dropped_pct > 10) {
      warning(paste0("Complete case analysis drops ", round(dropped_pct, 1), 
                    "% of observations"), call. = FALSE)
    }
    
    return(list(model_spec = model_spec, data = complete_data))
  }
  
  # Default return
  warning("Invalid method specified, returning original data and model", call. = FALSE)
  return(list(model_spec = model_spec, data = data))
}

#' Summarize sensitivity analyses
#' 
#' @param sensitivity_results Results from run_sensitivity_analyses
#' @param parameters Vector of parameter names to summarize (NULL for default)
#' @param format Output format (markdown, html, latex)
#' @return Summary of sensitivity analyses
#' @export
summarize_sensitivity <- function(sensitivity_results,
                                parameters = NULL,
                                format = c("markdown", "html", "latex")) {
  
  format <- match.arg(format)
  
  if (!inherits(sensitivity_results, "sensitivity_results")) {
    stop("sensitivity_results must be an object returned by run_sensitivity_analyses", call. = FALSE)
  }
  
  # Extract base model and sensitivity models
  base_model <- sensitivity_results$base_model
  sensitivity_models <- sensitivity_results$sensitivity_models
  analysis_type <- sensitivity_results$type
  
  # Default parameter is treatment effect if not specified
  if (is.null(parameters)) {
    # Try to find treatment parameter
    if (!is.null(base_model$model_spec)) {
      pred_terms <- strsplit(base_model$model_spec$predictors, "\\+|\\*")[[1]]
      pred_terms <- trimws(pred_terms)
      if (any(grepl("treatment", pred_terms, ignore.case = TRUE))) {
        parameters <- "b_treatment"
      } else if (any(grepl("arm", pred_terms, ignore.case = TRUE))) {
        parameters <- "b_arm"
      } else if (any(grepl("trt", pred_terms, ignore.case = TRUE))) {
        parameters <- "b_trt"
      }
    }
    
    if (is.null(parameters) && !is.null(base_model$model)) {
      # Try to find a key parameter
      posterior_names <- colnames(brms::posterior_samples(base_model$model))
      if (any(grepl("^b_", posterior_names))) {
        parameters <- posterior_names[grep("^b_", posterior_names)[1]]
      }
    }
    
    if (is.null(parameters)) {
      # Default to first parameter
      if (!is.null(base_model$model)) {
        parameters <- colnames(brms::posterior_samples(base_model$model))[1]
      } else {
        stop("Cannot determine default parameter for sensitivity summary", call. = FALSE)
      }
    }
  }
  
  # Create summary table for each parameter
  summary_tables <- list()
  
  for (param in parameters) {
    # Compare posteriors for this parameter
    if (!is.null(base_model$model)) {
      posterior_comparison <- compare_posteriors(
        c(list(base = base_model), sensitivity_models),
        parameter = param
      )
      
      # Create summary table
      summary_table <- data.frame(
        Model = posterior_comparison$model,
        Parameter = posterior_comparison$parameter,
        Mean = posterior_comparison$mean,
        SD = posterior_comparison$sd,
        `2.5%` = posterior_comparison$q2.5,
        `50%` = posterior_comparison$q50,
        `97.5%` = posterior_comparison$q97.5,
        `Diff from Base` = posterior_comparison$diff_from_base,
        `% Diff` = posterior_comparison$pct_diff_from_base,
        stringsAsFactors = FALSE
      )
      
      # Format numeric columns
      summary_table$Mean <- sprintf("%.4f", summary_table$Mean)
      summary_table$SD <- sprintf("%.4f", summary_table$SD)
      summary_table$`2.5%` <- sprintf("%.4f", summary_table$`2.5%`)
      summary_table$`50%` <- sprintf("%.4f", summary_table$`50%`)
      summary_table$`97.5%` <- sprintf("%.4f", summary_table$`97.5%`)
      summary_table$`Diff from Base` <- sprintf("%.4f", summary_table$`Diff from Base`)
      summary_table$`% Diff` <- sprintf("%.2f%%", summary_table$`% Diff`)
      
      # Replace NA with empty string for base model diff
      summary_table$`Diff from Base`[summary_table$Model == "base"] <- ""
      summary_table$`% Diff`[summary_table$Model == "base"] <- ""
      
      summary_tables[[param]] <- summary_table
    }
  }
  
  # Create comparison tables if available
  comparison_tables <- list()
  
  if (!is.null(sensitivity_results$comparisons$loo)) {
    loo_compare <- sensitivity_results$comparisons$loo
    
    # Create LOO comparison table
    loo_table <- as.data.frame(loo_compare)
    
    # Add model names
    loo_table$Model <- rownames(loo_compare)
    
    # Reorder columns
    loo_cols <- c("Model", colnames(loo_compare))
    loo_table <- loo_table[, loo_cols]
    
    # Format numeric columns
    for (col in setdiff(loo_cols, "Model")) {
      loo_table[[col]] <- sprintf("%.2f", loo_table[[col]])
    }
    
    comparison_tables$loo <- loo_table
  }
  
  # Create summary text
  summary_text <- paste0("# Sensitivity Analysis Summary\n\n",
                        "## Analysis Type: ", analysis_type, "\n\n")
  
  # Add model descriptions
  summary_text <- paste0(summary_text, "## Models\n\n")
  
  # Base model
  if (!is.null(base_model$model_spec)) {
    base_name <- if (!is.null(base_model$model_spec$model_name)) {
      base_model$model_spec$model_name
    } else {
      "Base Model"
    }
    
    summary_text <- paste0(summary_text, "### ", base_name, "\n\n",
                          "- Formula: ", build_formula_string(base_model$model_spec), "\n",
                          "- Family: ", base_model$model_spec$family, "(", base_model$model_spec$link, ")\n\n")
  }
  
  # Sensitivity models
  for (name in names(sensitivity_models)) {
    model <- sensitivity_models[[name]]
    
    if (!is.null(model$model_spec)) {
      summary_text <- paste0(summary_text, "### ", name, "\n\n",
                            "- Formula: ", build_formula_string(model$model_spec), "\n",
                            "- Family: ", model$model_spec$family, "(", model$model_spec$link, ")\n\n")
    }
  }
  
  # Add parameter summaries
  summary_text <- paste0(summary_text, "## Parameter Estimates\n\n")
  
  for (param in names(summary_tables)) {
    summary_text <- paste0(summary_text, "### Parameter: ", param, "\n\n")
    
    # Format table based on output format
    if (format == "markdown") {
      # Create markdown table
      table_text <- knitr::kable(summary_tables[[param]], format = "markdown")
      summary_text <- paste0(summary_text, table_text, "\n\n")
    } else if (format == "html") {
      # Create HTML table
      table_text <- knitr::kable(summary_tables[[param]], format = "html")
      summary_text <- paste0(summary_text, table_text, "\n\n")
    } else if (format == "latex") {
      # Create LaTeX table
      table_text <- knitr::kable(summary_tables[[param]], format = "latex")
      summary_text <- paste0(summary_text, table_text, "\n\n")
    }
  }
  
  # Add model comparisons if available
  if (length(comparison_tables) > 0) {
    summary_text <- paste0(summary_text, "## Model Comparisons\n\n")
    
    if (!is.null(comparison_tables$loo)) {
      summary_text <- paste0(summary_text, "### Leave-One-Out Cross-Validation\n\n",
                            "Lower LOOIC indicates better fit. Negative values in 'elpd_diff' indicate better models ",
                            "compared to the top model.\n\n")
      
      # Format table based on output format
      if (format == "markdown") {
        table_text <- knitr::kable(comparison_tables$loo, format = "markdown")
        summary_text <- paste0(summary_text, table_text, "\n\n")
      } else if (format == "html") {
        table_text <- knitr::kable(comparison_tables$loo, format = "html")
        summary_text <- paste0(summary_text, table_text, "\n\n")
      } else if (format == "latex") {
        table_text <- knitr::kable(comparison_tables$loo, format = "latex")
        summary_text <- paste0(summary_text, table_text, "\n\n")
      }
    }
  }
  
  # Add conclusion
  summary_text <- paste0(summary_text, "## Conclusion\n\n",
                        "This sensitivity analysis examined the robustness of the base model to changes in ",
                        analysis_type, " assumptions.\n\n")
  
  # Find largest effect
  max_diff <- 0
  max_diff_model <- ""
  
  for (param in names(summary_tables)) {
    tbl <- summary_tables[[param]]
    diff_col <- which(colnames(tbl) == "Diff from Base")
    
    # Skip if no diff column
    if (length(diff_col) == 0) next
    
    # Find max absolute difference
    diffs <- as.numeric(gsub("[^0-9\\.-]", "", tbl$`Diff from Base`))
    diffs <- diffs[!is.na(diffs)]
    
    if (length(diffs) > 0) {
      abs_diffs <- abs(diffs)
      max_idx <- which.max(abs_diffs)
      
      if (abs_diffs[max_idx] > max_diff) {
        max_diff <- abs_diffs[max_idx]
        max_diff_model <- tbl$Model[tbl$`Diff from Base` == sprintf("%.4f", diffs[max_idx])]
        max_diff_param <- param
      }
    }
  }
  
  if (max_diff > 0 && max_diff_model != "") {
    summary_text <- paste0(summary_text, "The largest sensitivity was observed in the '", 
                          max_diff_model, "' model for parameter '", max_diff_param, 
                          "', with a difference of ", sprintf("%.4f", max_diff), 
                          " from the base model.\n\n")
  }
  
  # Overall assessment
  if (max_diff <= 0.1) {
    summary_text <- paste0(summary_text, "The results appear to be robust to changes in ",
                          analysis_type, " assumptions, with only minor variations across models.\n")
  } else if (max_diff <= 0.3) {
    summary_text <- paste0(summary_text, "The results show moderate sensitivity to changes in ",
                          analysis_type, " assumptions. Caution should be exercised in interpretation.\n")
  } else {
    summary_text <- paste0(summary_text, "The results show substantial sensitivity to changes in ",
                          analysis_type, " assumptions. Findings should be interpreted with caution.\n")
  }
  
  # Return formatted summary
  return(summary_text)
}

#' Plot sensitivity analysis results
#' 
#' @param sensitivity_results Results from run_sensitivity_analyses
#' @param parameter Name of the parameter to plot
#' @param plot_type Type of plot to create (forest, ridgeline, scatter)
#' @return A ggplot object
#' @export
plot_sensitivity <- function(sensitivity_results,
                           parameter,
                           plot_type = c("forest", "ridgeline", "scatter")) {
  
  plot_type <- match.arg(plot_type)
  
  if (!inherits(sensitivity_results, "sensitivity_results")) {
    stop("sensitivity_results must be an object returned by run_sensitivity_analyses", call. = FALSE)
  }
  
  if (missing(parameter)) {
    stop("parameter must be specified", call. = FALSE)
  }
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Extract base model and sensitivity models
  base_model <- sensitivity_results$base_model
  sensitivity_models <- sensitivity_results$sensitivity_models
  analysis_type <- sensitivity_results$type
  
  # Compare posteriors for this parameter
  posterior_comparison <- compare_posteriors(
    c(list(base = base_model), sensitivity_models),
    parameter = parameter
  )
  
  # Extract raw posteriors
  posteriors <- attr(posterior_comparison, "posteriors")
  
  if (is.null(posteriors) || length(posteriors) == 0) {
    stop(paste("No posterior samples found for parameter:", parameter), call. = FALSE)
  }
  
  # Create plot based on type
  if (plot_type == "forest") {
    # Forest plot of parameter estimates
    
    # Create data frame for plotting
    plot_data <- data.frame(
      model = posterior_comparison$model,
      estimate = as.numeric(posterior_comparison$mean),
      lower = as.numeric(posterior_comparison$q2.5),
      upper = as.numeric(posterior_comparison$q97.5),
      stringsAsFactors = FALSE
    )
    
    # Order by model name
    plot_data$model <- factor(plot_data$model, 
                             levels = c("base", setdiff(plot_data$model, "base")))
    
    # Create forest plot
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = estimate, y = model)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = 0.2) +
      ggplot2::geom_point(size = 3, color = "blue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Sensitivity Analysis for", parameter),
        subtitle = paste("Analysis type:", analysis_type),
        x = "Estimate",
        y = NULL
      )
    
  } else if (plot_type == "ridgeline") {
    # Ridgeline plot of posterior distributions
    
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("ggridges package is required but not available", call. = FALSE)
    }
    
    # Combine posteriors into a data frame
    combined_data <- data.frame()
    
    for (model_name in names(posteriors)) {
      model_data <- data.frame(
        model = model_name,
        value = posteriors[[model_name]],
        stringsAsFactors = FALSE
      )
      combined_data <- rbind(combined_data, model_data)
    }
    
    # Order models with base first
    combined_data$model <- factor(combined_data$model, 
                                 levels = c("base", setdiff(unique(combined_data$model), "base")))
    
    # Create ridgeline plot
    plot <- ggplot2::ggplot(combined_data, ggplot2::aes(x = value, y = model, fill = model)) +
      ggridges::geom_density_ridges(alpha = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(
        title = paste("Posterior Distributions for", parameter),
        subtitle = paste("Analysis type:", analysis_type),
        x = "Parameter Value",
        y = NULL
      )
    
  } else if (plot_type == "scatter") {
    # Scatter plot of estimates vs standard deviations
    
    # Create data frame for plotting
    plot_data <- data.frame(
      model = posterior_comparison$model,
      estimate = as.numeric(posterior_comparison$mean),
      sd = as.numeric(posterior_comparison$sd),
      stringsAsFactors = FALSE
    )
    
    # Create scatter plot
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = estimate, y = sd, color = model == "base")) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_text(ggplot2::aes(label = model), vjust = -1) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_manual(values = c("grey40", "red"), guide = "none") +
      ggplot2::labs(
        title = paste("Sensitivity Analysis for", parameter),
        subtitle = paste("Analysis type:", analysis_type),
        x = "Estimate",
        y = "Standard Deviation"
      )
  }
  
  return(plot)
}