#' Generate comprehensive clinical trial analysis reports
#'
#' Creates publication-ready reports from Bayesian clinical trial models with 
#' automated statistical summaries, visualizations, and clinical interpretations.
#' Essential for clinical documentation, investigator communications, and 
#' scientific publications. Supports multiple output formats and customizable
#' report templates.
#'
#' @param model A BayesianModel object from \code{\link{fit_model}} containing
#'   fitted results. Model must have completed MCMC sampling with convergence.
#'   The report will automatically extract:
#'   \itemize{
#'     \item Model specification and formula
#'     \item Posterior parameter estimates with credible intervals
#'     \item Treatment effect summaries and probabilities
#'     \item Convergence diagnostics and model fit assessments
#'   }
#'   
#' @param template Character string specifying the report template. Available options:
#'   \itemize{
#'     \item \code{"basic_report"} - Standard analysis report with treatment effects,
#'       posterior summaries, and model diagnostics. Suitable for most clinical trials.
#'     \item \code{"diagnostics_report"} - Detailed convergence and model fit assessment.
#'       Essential for model validation and troubleshooting.
#'     \item \code{"sensitivity_report"} - Template for sensitivity analysis results.
#'       Use with sensitivity analysis outputs.
#'     \item \code{"subgroup_report"} - Template for effect modifier analysis.
#'       Use with subgroup analysis results.
#'   }
#'   Custom templates can be created following R Markdown parameterized report format.
#'   
#' @param output_file Character string specifying output file path (optional).
#'   If not provided, automatically generates timestamped filename:
#'   \itemize{
#'     \item Format: \code{[model_name]_report_[YYYYMMDD_HHMMSS].[extension]}
#'     \item Example: \code{"efficacy_model_report_20241201_143022.html"}
#'   }
#'   Include full path if saving to specific directory.
#'   
#' @param format Character string specifying output format:
#'   \itemize{
#'     \item \code{"html"} - Interactive HTML report with embedded plots (default).
#'       Best for sharing via email or web platforms. Self-contained.
#'     \item \code{"pdf"} - Publication-ready PDF document. Requires LaTeX installation.
#'       Ideal for clinical publications and formal documentation.
#'     \item \code{"word"} - Microsoft Word document for collaborative editing.
#'       Compatible with institutional review processes.
#'   }
#'   
#' @param include_plots Logical indicating whether to include visualizations:
#'   \itemize{
#'     \item \code{TRUE} - Includes posterior distribution plots, model fit diagnostics,
#'       and treatment effect visualizations (recommended)
#'     \item \code{FALSE} - Text-only report with tables and numerical summaries
#'   }
#'   
#' @param include_diagnostics Logical indicating whether to include MCMC diagnostics:
#'   \itemize{
#'     \item \code{TRUE} - Includes R-hat values, effective sample sizes, 
#'       divergent transitions, and convergence assessments (recommended)
#'     \item \code{FALSE} - Omits technical diagnostic information
#'   }
#'
#' @return Character string with path to the generated report file.
#'   File is ready for sharing, submission, or publication.
#'
#' @details
#' \strong{Report Contents:}
#' 
#' \strong{Basic Report Template Includes:}
#' \itemize{
#'   \item Model specification summary (formula, family, link function)
#'   \item Posterior parameter estimates with 95% credible intervals
#'   \item Treatment effect estimates with probability of benefit
#'   \item Model convergence diagnostics (if requested)
#'   \item Posterior distribution visualizations (if requested)
#'   \item Model fit assessment plots
#'   \item Clinical interpretation and conclusions
#' }
#' 
#' \strong{Clinical Communication Features:}
#' \itemize{
#'   \item Automatic calculation of probability of positive treatment effect
#'   \item Clinical significance thresholds and interpretations
#'   \item Credible intervals (not confidence intervals) for Bayesian results
#'   \item Clear differentiation between statistical and clinical significance
#' }
#' 
#' \strong{Quality Assurance:}
#' \itemize{
#'   \item Timestamped generation for audit trails
#'   \item Reproducible analysis documentation
#'   \item Model diagnostic requirements for validation
#'   \item Clear statistical methodology descriptions
#' }
#' 
#' \strong{Technical Requirements:}
#' \itemize{
#'   \item Requires rmarkdown package for report generation
#'   \item PDF output requires LaTeX installation (TinyTeX recommended)
#'   \item HTML output is self-contained with embedded images
#'   \item Word output compatible with Office 2016+
#' }
#'
#' @examples
#' \dontrun{
#' # Fit a clinical trial model
#' data(synthetic_trial_data)
#' priors <- default_priors()
#' model_spec <- create_model_spec(
#'   outcome ~ treatment + age + sex + baseline_score,
#'   family = "gaussian",
#'   model_name = "primary_efficacy"
#' )
#' 
#' model <- fit_model(model_spec, synthetic_trial_data, priors)
#' 
#' # Generate basic HTML report
#' report_path <- generate_report(model)
#' browseURL(report_path)  # Open in browser
#' 
#' # Generate PDF report for formal documentation
#' pdf_report <- generate_report(
#'   model, 
#'   template = "basic_report",
#'   output_file = "formal_efficacy_analysis.pdf",
#'   format = "pdf",
#'   include_plots = TRUE,
#'   include_diagnostics = TRUE
#' )
#' 
#' # Generate diagnostics-focused report
#' diag_report <- generate_report(
#'   model,
#'   template = "diagnostics_report", 
#'   format = "html"
#' )
#' 
#' # Generate Word report for clinical team review
#' clinical_report <- generate_report(
#'   model,
#'   output_file = "clinical_summary.docx",
#'   format = "word",
#'   include_plots = TRUE,
#'   include_diagnostics = FALSE  # Simplified for clinical audience
#' )
#' 
#' # Batch report generation for multiple models
#' models_list <- list(
#'   primary = primary_model,
#'   safety = safety_model,
#'   sensitivity = sensitivity_model
#' )
#' 
#' report_paths <- sapply(names(models_list), function(name) {
#'   generate_report(
#'     models_list[[name]],
#'     output_file = paste0(name, "_analysis.html"),
#'     format = "html"
#'   )
#' })
#' 
#' # Verify all reports generated successfully
#' all(file.exists(report_paths))
#' }
#'
#' @seealso 
#' \code{\link{fit_model}} for fitting Bayesian models,
#' \code{\link{sensitivity_report}} for sensitivity analysis reports,
#' \code{\link{subgroup_report}} for subgroup analysis reports,
#' \code{\link{reproducibility_report}} for workflow documentation
#'
#' @export
generate_report <- function(model,
                           template = "basic_report",
                           output_file = NULL,
                           format = c("html", "pdf", "word"),
                           include_plots = TRUE,
                           include_diagnostics = TRUE) {
  
  format <- match.arg(format)
  
  # Check if rmarkdown is available
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rmarkdown package is required but not available", call. = FALSE)
  }
  
  # Check model
  if (!inherits(model, "BayesianModel")) {
    stop("model must be a BayesianModel object", call. = FALSE)
  }
  
  if (is.null(model$model)) {
    stop("Model has not been fitted", call. = FALSE)
  }
  
  # Set default output file if not provided
  if (is.null(output_file)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    model_name <- if (!is.null(model$model_spec) && !is.null(model$model_spec$model_name)) {
      model$model_spec$model_name
    } else {
      "model"
    }
    
    extension <- switch(format,
                       html = ".html",
                       pdf = ".pdf",
                       word = ".docx")
    
    output_file <- paste0(model_name, "_report_", timestamp, extension)
  }
  
  # Get template file
  template_file <- system.file("templates", paste0(template, ".Rmd"), package = "bayestrials")
  
  if (!file.exists(template_file)) {
    # If not found in package, check current directory
    template_file <- file.path("templates", paste0(template, ".Rmd"))
    
    if (!file.exists(template_file)) {
      # If still not found, check for standard templates
      standard_templates <- c(
        "basic_report",
        "diagnostics_report",
        "sensitivity_report",
        "subgroup_report"
      )
      
      if (template %in% standard_templates) {
        # Use built-in template content
        template_content <- get_template_content(template)
        
        # Create temporary template file
        template_file <- tempfile(pattern = "template_", fileext = ".Rmd")
        writeLines(template_content, template_file)
      } else {
        stop(paste("Template", template, "not found"), call. = FALSE)
      }
    }
  }
  
  # Create parameter list for rendering
  params <- list(
    model = model,
    include_plots = include_plots,
    include_diagnostics = include_diagnostics,
    timestamp = Sys.time()
  )
  
  # Render report
  rmarkdown::render(
    input = template_file,
    output_file = output_file,
    output_format = switch(format,
                          html = rmarkdown::html_document(),
                          pdf = rmarkdown::pdf_document(),
                          word = rmarkdown::word_document()),
    params = params
  )
  
  # Return path to generated report
  return(output_file)
}

#' Get content for standard templates
#' 
#' @param template_name Name of the template
#' @return Template content as character string
#' @keywords internal
get_template_content <- function(template_name) {
  if (template_name == "basic_report") {
    return('---
title: "Bayesian Model Report"
author: "bayestrials"
date: "`r params$timestamp`"
output: html_document
params:
  model: NULL
  include_plots: true
  include_diagnostics: true
  timestamp: !r Sys.time()
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(bayestrials)
library(ggplot2)
library(knitr)

# Ensure model is available
if (is.null(params$model)) {
  stop("No model provided")
}

model <- params$model
```

## Model Overview

```{r model-info}
# Extract model information
model_name <- if (!is.null(model$model_spec) && !is.null(model$model_spec$model_name)) {
  model$model_spec$model_name
} else {
  "Unnamed Model"
}

family <- model$model_spec$family
link <- model$model_spec$link
formula <- build_formula_string(model$model_spec)
```

**Model Name**: `r model_name`

**Family**: `r family` with `r link` link

**Formula**: `r formula`

## Summary Statistics

```{r summary}
# Get posterior summary
post_summary <- posterior_summary(model)

# Display summary table
kable(post_summary, digits = 4)
```

## Treatment Effect

```{r treatment-effect}
# Try to find treatment parameter
treat_param <- NULL
if (!is.null(model$model_spec)) {
  pred_terms <- strsplit(model$model_spec$predictors, "\\+|\\*")[[1]]
  pred_terms <- trimws(pred_terms)
  if (any(grepl("treatment", pred_terms, ignore.case = TRUE))) {
    treat_param <- "b_treatment"
  } else if (any(grepl("arm", pred_terms, ignore.case = TRUE))) {
    treat_param <- "b_arm"
  } else if (any(grepl("trt", pred_terms, ignore.case = TRUE))) {
    treat_param <- "b_trt"
  }
}

if (!is.null(treat_param)) {
  # Filter for treatment parameter
  treat_summary <- post_summary[grep(treat_param, post_summary$parameter), ]
  
  # Display treatment effect table
  kable(treat_summary, digits = 4, caption = "Treatment Effect Estimates")
  
  # Calculate probability of positive effect
  if (nrow(treat_summary) > 0) {
    post_samples <- extract_posterior(model, treat_param)
    prob_positive <- mean(post_samples[[1]] > 0)
    cat("Probability of positive treatment effect:", sprintf("%.3f", prob_positive))
  }
}
```

```{r treatment-plot, eval=params$include_plots && !is.null(treat_param)}
if (params$include_plots && !is.null(treat_param) && nrow(treat_summary) > 0) {
  # Plot posterior distribution
  plot_posterior(post_samples)
}
```

## Model Diagnostics

```{r diagnostics, eval=params$include_diagnostics}
if (params$include_diagnostics) {
  # Check diagnostics
  diag_results <- check_diagnostics(model, quiet = TRUE)
  
  # Display diagnostic summary
  cat("Diagnostic Status:", diag_results$status, "\n\n")
  
  if (length(diag_results$warnings) > 0) {
    cat("Warnings:\n")
    for (warning in diag_results$warnings) {
      cat("- ", warning, "\n")
    }
  } else {
    cat("No warnings detected.\n")
  }
  
  # Show Rhat values
  if (!is.null(diag_results$rhat)) {
    high_rhat <- diag_results$rhat[diag_results$rhat > 1.01]
    if (length(high_rhat) > 0) {
      cat("\nHigh Rhat values (> 1.01):\n")
      print(high_rhat)
    } else {
      cat("\nAll Rhat values look good (≤ 1.01).\n")
    }
  }
  
  # Show ESS values
  if (!is.null(diag_results$ess)) {
    low_ess <- diag_results$ess[diag_results$ess < 0.1]
    if (length(low_ess) > 0) {
      cat("\nLow ESS ratios (< 0.1):\n")
      print(low_ess)
    } else {
      cat("\nAll ESS ratios look good (≥ 0.1).\n")
    }
  }
}
```

```{r diagnostic-plots, eval=params$include_diagnostics && params$include_plots, fig.width=10, fig.height=8}
if (params$include_diagnostics && params$include_plots) {
  # Show diagnostic plots
  if (!is.null(diag_results$plots$rhat)) {
    cat("## Rhat Plot\n\n")
    print(diag_results$plots$rhat)
  }
  
  if (!is.null(diag_results$plots$ess)) {
    cat("## Effective Sample Size Plot\n\n")
    print(diag_results$plots$ess)
  }
  
  if (!is.null(diag_results$plots$ppc_dens)) {
    cat("## Posterior Predictive Check\n\n")
    print(diag_results$plots$ppc_dens)
  }
}
```

## Model Fit

```{r model-fit, eval=params$include_plots, fig.width=10, fig.height=6}
if (params$include_plots) {
  # Generate model fit plot
  fit_plot <- plot_model_fit(model, type = "predicted_vs_observed")
  print(fit_plot)
}
```

## Conclusion

This report presents the results of a Bayesian analysis using the `r family` family with a `r link` link function. The model formula specified was: `r formula`. 

```{r conclusion}
if (!is.null(treat_param) && nrow(treat_summary) > 0) {
  # Generate conclusion about treatment effect
  effect_size <- treat_summary$mean[1]
  ci_lower <- treat_summary$q2.5[1]
  ci_upper <- treat_summary$q97.5[1]
  
  conclusion <- paste0(
    "The estimated treatment effect is ", sprintf("%.3f", effect_size),
    " (95% CI: [", sprintf("%.3f", ci_lower), ", ", sprintf("%.3f", ci_upper), "])."
  )
  
  if (ci_lower > 0) {
    conclusion <- paste0(conclusion, " There is strong evidence that the treatment has a positive effect.")
  } else if (ci_upper < 0) {
    conclusion <- paste0(conclusion, " There is strong evidence that the treatment has a negative effect.")
  } else if (prob_positive > 0.9) {
    conclusion <- paste0(conclusion, " There is suggestive evidence that the treatment has a positive effect (Prob > 0: ", sprintf("%.3f", prob_positive), ").")
  } else if (prob_positive < 0.1) {
    conclusion <- paste0(conclusion, " There is suggestive evidence that the treatment has a negative effect (Prob > 0: ", sprintf("%.3f", prob_positive), ").")
  } else {
    conclusion <- paste0(conclusion, " The evidence does not strongly support either a positive or negative treatment effect.")
  }
  
  cat(conclusion)
}
```

---

Report generated with `bayestrials` on `r format(params$timestamp, "%Y-%m-%d %H:%M:%S")`.
')
  } else if (template_name == "diagnostics_report") {
    return('---
title: "Bayesian Model Diagnostics Report"
author: "bayestrials"
date: "`r params$timestamp`"
output: html_document
params:
  model: NULL
  include_plots: true
  include_diagnostics: true
  timestamp: !r Sys.time()
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(bayestrials)
library(ggplot2)
library(knitr)

# Ensure model is available
if (is.null(params$model)) {
  stop("No model provided")
}

model <- params$model
```

## Model Information

```{r model-info}
# Extract model information
model_name <- if (!is.null(model$model_spec) && !is.null(model$model_spec$model_name)) {
  model$model_spec$model_name
} else {
  "Unnamed Model"
}

family <- model$model_spec$family
link <- model$model_spec$link
formula <- build_formula_string(model$model_spec)
```

**Model Name**: `r model_name`

**Family**: `r family` with `r link` link

**Formula**: `r formula`

**Algorithm**: `r model$model$algorithm`

**Chains**: `r model$model$fit@sim$chains`

**Iterations per chain**: `r model$model$fit@sim$iter`

**Warmup**: `r model$model$fit@sim$warmup`

## MCMC Diagnostics

```{r run-diagnostics}
# Run full diagnostics
diag_results <- check_diagnostics(model, quiet = TRUE)
```

### Convergence Summary

**Diagnostic Status**: `r diag_results$status`

```{r diagnostic-summary}
# Display diagnostic warnings
if (length(diag_results$warnings) > 0) {
  cat("**Warnings**:\n\n")
  for (warning in diag_results$warnings) {
    cat("- ", warning, "\n")
  }
} else {
  cat("No convergence issues detected.\n")
}
```

### Rhat Values

Rhat (or potential scale reduction factor) measures chain convergence. Values close to 1.0 indicate good convergence, while values above 1.01 suggest potential convergence issues.

```{r rhat-summary}
# Show Rhat statistics
if (!is.null(diag_results$rhat)) {
  rhat_stats <- summary(diag_results$rhat)
  kable(t(as.matrix(rhat_stats)), digits = 4, caption = "Rhat Summary Statistics")
  
  high_rhat <- diag_results$rhat[diag_results$rhat > 1.01]
  if (length(high_rhat) > 0) {
    cat("\n**Parameters with high Rhat values (> 1.01)**:\n\n")
    high_rhat_df <- data.frame(
      parameter = names(high_rhat),
      rhat = high_rhat
    )
    kable(high_rhat_df, digits = 4)
  } else {
    cat("\nAll Rhat values look good (≤ 1.01).\n")
  }
}
```

```{r rhat-plot, eval=params$include_plots, fig.width=10, fig.height=6}
if (params$include_plots) {
  plot_diagnostics(model, type = "rhat")
}
```

### Effective Sample Size

Effective Sample Size (ESS) measures the effective number of independent samples, accounting for autocorrelation in the chains. Higher values are better. ESS ratios below 0.1 (10% of total samples) may indicate issues.

```{r ess-summary}
# Show ESS statistics
if (!is.null(diag_results$ess)) {
  ess_stats <- summary(diag_results$ess)
  kable(t(as.matrix(ess_stats)), digits = 4, caption = "ESS Ratio Summary Statistics")
  
  low_ess <- diag_results$ess[diag_results$ess < 0.1]
  if (length(low_ess) > 0) {
    cat("\n**Parameters with low ESS ratios (< 0.1)**:\n\n")
    low_ess_df <- data.frame(
      parameter = names(low_ess),
      ess_ratio = low_ess
    )
    kable(low_ess_df, digits = 4)
  } else {
    cat("\nAll ESS ratios look good (≥ 0.1).\n")
  }
}
```

```{r ess-plot, eval=params$include_plots, fig.width=10, fig.height=6}
if (params$include_plots) {
  plot_diagnostics(model, type = "ess")
}
```

### MCMC Trace Plots

Trace plots show the history of parameter values across iterations. Well-behaved chains should mix well and be stationary.

```{r trace-plots, eval=params$include_plots, fig.width=10, fig.height=8}
if (params$include_plots) {
  # Get a few key parameters for trace plots
  post_samples <- brms::posterior_samples(model$model)
  key_params <- colnames(post_samples)[1:min(6, ncol(post_samples))]
  
  plot_diagnostics(model, type = "trace", parameters = key_params)
}
```

### Posterior Intervals

```{r intervals-plot, eval=params$include_plots, fig.width=10, fig.height=8}
if (params$include_plots) {
  # Show interval plot for key parameters
  plot_diagnostics(model, type = "intervals")
}
```

## Posterior Predictive Checks

Posterior predictive checks assess model fit by comparing the observed data to predictions from the posterior distribution.

```{r ppc-plots, eval=params$include_plots, fig.width=10, fig.height=6}
if (params$include_plots) {
  # Density overlay
  cat("### Density Overlay\n\n")
  print(diag_results$plots$ppc_dens)
  
  # Scatter plot
  cat("\n### Scatter Plot\n\n")
  print(diag_results$plots$ppc_scatter)
}
```

## Model Fit Statistics

```{r fit-stats}
# Calculate LOO
loo_result <- tryCatch({
  model_fit_stats(model, type = "loo")
}, error = function(e) NULL)

# Calculate WAIC
waic_result <- tryCatch({
  model_fit_stats(model, type = "waic")
}, error = function(e) NULL)

# Display results
if (!is.null(loo_result)) {
  cat("### LOO (Leave-One-Out Cross-Validation)\n\n")
  print(loo_result)
}

if (!is.null(waic_result)) {
  cat("\n### WAIC (Widely Applicable Information Criterion)\n\n")
  print(waic_result)
}
```

## Residual Analysis

```{r residual-plots, eval=params$include_plots, fig.width=10, fig.height=8}
if (params$include_plots) {
  # Show residual plots
  cat("### Residual Plots\n\n")
  print(plot_model_fit(model, type = "residuals"))
}
```

## Posterior Predictive P-Value

The posterior predictive p-value assesses model fit by comparing a test statistic calculated on the observed data to the same statistic calculated on data generated from the model.

```{r pp-pvalue}
# Calculate posterior predictive p-value
pp_pvalue <- tryCatch({
  posterior_predictive_pvalue(model)
}, error = function(e) NULL)

if (!is.null(pp_pvalue)) {
  cat("Posterior predictive p-value:", sprintf("%.4f", pp_pvalue$p_value), "\n\n")
  
  if (pp_pvalue$p_value < 0.05 || pp_pvalue$p_value > 0.95) {
    cat("The p-value is extreme (< 0.05 or > 0.95), suggesting potential model misfit.\n")
  } else {
    cat("The p-value is not extreme, suggesting adequate model fit.\n")
  }
}
```

## Divergent Transitions

Divergent transitions occur when the MCMC sampler encounters regions of the posterior distribution that are difficult to explore, often due to complex geometry.

```{r divergences}
# Show divergence information
if (!is.null(diag_results$divergences)) {
  cat("Number of divergent transitions:", diag_results$divergences, "\n\n")
  
  if (diag_results$divergences > 0) {
    cat("The presence of divergent transitions suggests potential issues with the posterior geometry. ",
        "Consider increasing adapt_delta, simplifying the model, or reparameterizing.\n")
  } else {
    cat("No divergent transitions detected, which is good.\n")
  }
}
```

## Conclusion and Recommendations

```{r conclusion}
# Generate recommendations based on diagnostics
if (diag_results$status == "good") {
  cat("The model diagnostics indicate good convergence and no major issues. The model appears to be reliable for inference.\n")
} else {
  cat("The model diagnostics indicate some issues that should be addressed:\n\n")
  
  if (!is.null(diag_results$rhat) && any(diag_results$rhat > 1.01)) {
    cat("- **Convergence Issues**: Some parameters have high Rhat values, suggesting chain convergence problems. Consider running more iterations or adjusting priors.\n")
  }
  
  if (!is.null(diag_results$ess) && any(diag_results$ess < 0.1)) {
    cat("- **Low Effective Sample Size**: Some parameters have low ESS ratios, indicating high autocorrelation. Consider thinning, running more iterations, or reparameterizing the model.\n")
  }
  
  if (!is.null(diag_results$divergences) && diag_results$divergences > 0) {
    cat("- **Divergent Transitions**: The model has divergent transitions, which may lead to biased inference. Try increasing adapt_delta, simplifying the model, or reparameterizing.\n")
  }
  
  if (!is.null(pp_pvalue) && (pp_pvalue$p_value < 0.05 || pp_pvalue$p_value > 0.95)) {
    cat("- **Poor Model Fit**: The posterior predictive p-value suggests potential model misfit. Consider revising the model structure or checking for data issues.\n")
  }
}
```

---

Report generated with `bayestrials` on `r format(params$timestamp, "%Y-%m-%d %H:%M:%S")`.
')
  } else if (template_name == "sensitivity_report") {
    return('---
title: "Sensitivity Analysis Report"
author: "bayestrials"
date: "`r params$timestamp`"
output: html_document
params:
  model: NULL
  include_plots: true
  include_diagnostics: true
  timestamp: !r Sys.time()
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(bayestrials)
library(ggplot2)
library(knitr)

# Ensure model is available
if (is.null(params$model)) {
  stop("No model provided")
}

# Check if model is a sensitivity_results object
if (!inherits(params$model, "sensitivity_results")) {
  stop("Model must be a sensitivity_results object")
}

sens_results <- params$model
```

## Sensitivity Analysis Overview

This report presents the results of a sensitivity analysis to evaluate the robustness of the Bayesian model to variations in assumptions.

```{r analysis-info}
# Extract information
analysis_type <- sens_results$type
base_model <- sens_results$base_model
sensitivity_models <- sens_results$sensitivity_models
```

**Analysis Type**: `r analysis_type`

**Base Model**: `r if (!is.null(base_model$model_spec) && !is.null(base_model$model_spec$model_name)) base_model$model_spec$model_name else "Unnamed Model"`

**Number of Sensitivity Models**: `r length(sensitivity_models)`

## Model Descriptions

### Base Model

```{r base-model-info}
# Extract base model information
base_name <- if (!is.null(base_model$model_spec) && !is.null(base_model$model_spec$model_name)) {
  base_model$model_spec$model_name
} else {
  "Base Model"
}

base_family <- base_model$model_spec$family
base_link <- base_model$model_spec$link
base_formula <- build_formula_string(base_model$model_spec)
```

**Model Name**: `r base_name`

**Family**: `r base_family` with `r base_link` link

**Formula**: `r base_formula`

### Sensitivity Models

```{r sensitivity-model-info}
# Create table of sensitivity models
sens_table <- data.frame(
  Model = character(0),
  Family = character(0),
  Link = character(0),
  Formula = character(0),
  stringsAsFactors = FALSE
)

for (model_name in names(sensitivity_models)) {
  model <- sensitivity_models[[model_name]]
  
  if (!is.null(model) && !is.null(model$model_spec)) {
    sens_table <- rbind(sens_table, data.frame(
      Model = model_name,
      Family = model$model_spec$family,
      Link = model$model_spec$link,
      Formula = build_formula_string(model$model_spec),
      stringsAsFactors = FALSE
    ))
  }
}

# Display table
kable(sens_table)
```

## Parameter Estimates

```{r find-parameters}
# Try to find treatment parameter
treat_param <- NULL
if (!is.null(base_model$model_spec)) {
  pred_terms <- strsplit(base_model$model_spec$predictors, "\\+|\\*")[[1]]
  pred_terms <- trimws(pred_terms)
  if (any(grepl("treatment", pred_terms, ignore.case = TRUE))) {
    treat_param <- "b_treatment"
  } else if (any(grepl("arm", pred_terms, ignore.case = TRUE))) {
    treat_param <- "b_arm"
  } else if (any(grepl("trt", pred_terms, ignore.case = TRUE))) {
    treat_param <- "b_trt"
  }
}

# Use first parameter if treatment parameter not found
if (is.null(treat_param) && !is.null(base_model$model)) {
  post_samples <- brms::posterior_samples(base_model$model)
  if (ncol(post_samples) > 0) {
    treat_param <- colnames(post_samples)[1]
  }
}
```

```{r parameter-comparison, eval=!is.null(treat_param)}
if (!is.null(treat_param)) {
  # Create comparison for this parameter if available in sens_results
  if (is.null(sens_results$comparisons$treatment) && !is.null(base_model$model)) {
    # Create comparison manually
    comparison <- compare_posteriors(
      c(list(base = base_model), sensitivity_models),
      parameter = treat_param
    )
  } else {
    # Use existing comparison
    comparison <- sens_results$comparisons$treatment
  }
  
  if (!is.null(comparison)) {
    # Create formatted table
    comp_table <- data.frame(
      Model = comparison$model,
      Estimate = sprintf("%.4f", comparison$mean),
      SD = sprintf("%.4f", comparison$sd),
      `2.5%` = sprintf("%.4f", comparison$q2.5),
      `97.5%` = sprintf("%.4f", comparison$q97.5),
      `Diff from Base` = ifelse(comparison$model == "base", 
                               "", 
                               sprintf("%.4f", comparison$diff_from_base)),
      `% Diff` = ifelse(comparison$model == "base", 
                       "", 
                       sprintf("%.2f%%", comparison$pct_diff_from_base)),
      stringsAsFactors = FALSE
    )
    
    # Display table
    kable(comp_table, caption = paste("Estimates for", treat_param))
  }
}
```

```{r forest-plot, eval=!is.null(treat_param) && params$include_plots, fig.width=10, fig.height=6}
if (!is.null(treat_param) && params$include_plots && !is.null(comparison)) {
  # Create forest plot
  forest_data <- data.frame(
    trial = comparison$model,
    estimate = comparison$mean,
    ci_lower = comparison$q2.5,
    ci_upper = comparison$q97.5,
    stringsAsFactors = FALSE
  )
  
  plot_forest(forest_data, parameter = treat_param)
}
```

```{r ridgeline-plot, eval=!is.null(treat_param) && params$include_plots, fig.width=10, fig.height=6}
if (!is.null(treat_param) && params$include_plots && !is.null(comparison)) {
  # Create ridgeline plot
  plot_sensitivity(sens_results, treat_param, plot_type = "ridgeline")
}
```

## Model Comparisons

```{r model-comparison}
# Display model comparison metrics if available
if (!is.null(sens_results$comparisons$loo)) {
  cat("### Leave-One-Out Cross-Validation Comparison\n\n")
  kable(sens_results$comparisons$loo, digits = 2)
  cat("\nLower LOOIC values indicate better model fit.\n")
}
```

## Sensitivity Assessment

```{r assessment}
if (!is.null(comparison)) {
  # Calculate max absolute percentage difference
  if ("pct_diff_from_base" %in% names(comparison)) {
    diff_values <- comparison$pct_diff_from_base[!is.na(comparison$pct_diff_from_base)]
    max_diff <- max(abs(diff_values), na.rm = TRUE)
    
    cat("Maximum absolute percentage difference from base model:", sprintf("%.2f%%", max_diff), "\n\n")
    
    # Assess sensitivity
    if (max_diff <= 10) {
      cat("**Assessment**: The results appear to be robust to changes in", analysis_type, 
          "assumptions, with only minor variations across models (< 10% difference).\n")
    } else if (max_diff <= 30) {
      cat("**Assessment**: The results show moderate sensitivity to changes in", analysis_type, 
          "assumptions. Differences of 10-30% suggest caution in interpretation.\n")
    } else {
      cat("**Assessment**: The results show substantial sensitivity to changes in", analysis_type, 
          "assumptions. Differences > 30% indicate strong dependence on", analysis_type, "assumptions.\n")
    }
  }
}
```

## Conclusion

```{r conclusion}
cat("This sensitivity analysis examined the robustness of the base model to changes in", 
    analysis_type, "assumptions. ")

if (!is.null(comparison)) {
  # Generate conclusion based on sensitivity type and results
  if (analysis_type == "prior") {
    if (max_diff <= 10) {
      cat("The results are robust to different prior specifications, suggesting that the data is strongly informative and dominates the prior influence.")
    } else if (max_diff <= 30) {
      cat("The results show moderate sensitivity to prior specifications, suggesting that the prior has some influence on the posterior. Consider using weakly informative priors or collecting more data.")
    } else {
      cat("The results are highly sensitive to prior choices, indicating that the data may not be strongly informative. Interpretation should be done with caution, recognizing the strong influence of prior assumptions.")
    }
  } else if (analysis_type == "model") {
    if (max_diff <= 10) {
      cat("The results are consistent across different model specifications, suggesting that the key findings are robust to modeling decisions.")
    } else if (max_diff <= 30) {
      cat("The results show moderate variation across model specifications, suggesting that some modeling choices influence the conclusions. Consider model averaging or reporting results from multiple models.")
    } else {
      cat("The results vary substantially depending on model specification, indicating that modeling decisions strongly impact conclusions. Results should be interpreted with caution and model uncertainty should be emphasized.")
    }
  } else if (analysis_type == "missing_data") {
    if (max_diff <= 10) {
      cat("The results are robust to different missing data mechanisms, suggesting that missing data handling has minimal impact on conclusions.")
    } else if (max_diff <= 30) {
      cat("The results show moderate sensitivity to missing data assumptions, suggesting that the handling of missing data influences conclusions to some degree. Multiple imputation or modeling the missing data mechanism explicitly may be warranted.")
    } else {
      cat("The results are highly sensitive to missing data assumptions, indicating that conclusions depend strongly on how missing data is handled. This suggests potential bias due to missing data, requiring careful consideration.")
    }
  }
} else {
  cat("Due to limitations in the available results, a detailed assessment could not be provided. Further analysis is recommended.")
}
```

---

Report generated with `bayestrials` on `r format(params$timestamp, "%Y-%m-%d %H:%M:%S")`.
')
  } else if (template_name == "subgroup_report") {
    return('---
title: "Subgroup Analysis Report"
author: "bayestrials"
date: "`r params$timestamp`"
output: html_document
params:
  model: NULL
  include_plots: true
  include_diagnostics: true
  timestamp: !r Sys.time()
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(bayestrials)
library(ggplot2)
library(knitr)

# Ensure model is available
if (is.null(params$model)) {
  stop("No model provided")
}

# Check if model is an effect_modifier_results object
if (!inherits(params$model, "effect_modifier_results")) {
  stop("Model must be an effect_modifier_results object")
}

effect_results <- params$model
```

## Subgroup Analysis Overview

This report presents the results of a subgroup analysis to identify effect modifiers in a Bayesian model.

```{r analysis-info}
# Extract information
method <- effect_results$method
treatment_var <- effect_results$treatment_var
candidate_modifiers <- effect_results$candidate_modifiers
credible_modifiers <- effect_results$credible_modifiers
```

**Analysis Method**: `r method`

**Treatment Variable**: `r treatment_var`

**Candidate Modifiers**: `r paste(candidate_modifiers, collapse = ", ")`

**Number of Candidate Modifiers**: `r length(candidate_modifiers)`

**Number of Credible Modifiers**: `r length(credible_modifiers)`

## Identified Effect Modifiers

```{r identified-modifiers}
if (length(credible_modifiers) > 0) {
  cat("The following variables were identified as credible effect modifiers:\n\n")
  
  for (i in seq_along(credible_modifiers)) {
    modifier <- credible_modifiers[i]
    cat(i, ". **", modifier, "**\n", sep = "")
  }
} else {
  cat("No credible effect modifiers were identified among the candidate variables.\n")
}
```

## Effect Modifier Details

```{r modifier-details}
for (modifier in credible_modifiers) {
  cat("### ", modifier, "\n\n", sep = "")
  
  mod_effect <- effect_results$modifier_effects[[modifier]]
  
  if (method == "interaction") {
    # For standard interaction method
    cat("- **Interaction Parameter**: ", mod_effect$interaction_param, "\n", sep = "")
    cat("- **Estimate**: ", sprintf("%.4f", mod_effect$estimate), "\n", sep = "")
    cat("- **95% CI**: [", sprintf("%.4f, %.4f", mod_effect$ci_lower, mod_effect$ci_upper), "]\n\n", sep = "")
    
    # Calculate effect size interpretation
    if (abs(mod_effect$estimate) < 0.2) {
      effect_size_desc <- "small"
    } else if (abs(mod_effect$estimate) < 0.5) {
      effect_size_desc <- "moderate"
    } else {
      effect_size_desc <- "large"
    }
    
    cat("This represents a ", effect_size_desc, " effect modification. ", sep = "")
    
    if (mod_effect$estimate > 0) {
      cat("The positive interaction indicates that the treatment effect increases as ", modifier, " increases.\n\n", sep = "")
    } else {
      cat("The negative interaction indicates that the treatment effect decreases as ", modifier, " increases.\n\n", sep = "")
    }
    
  } else if (method == "bart" || method == "continuous") {
    # For BART or continuous methods
    te_range <- range(mod_effect$treatment_effects$mean_te)
    
    cat("- **Treatment effect varies by ", modifier, " values**\n", sep = "")
    cat("- **Range of treatment effects**: [", 
        sprintf("%.4f, %.4f", te_range[1], te_range[2]), "]\n\n", sep = "")
    
    # Calculate max difference
    max_diff <- diff(te_range)
    
    # Calculate effect size interpretation
    if (abs(max_diff) < 0.2) {
      effect_size_desc <- "small"
    } else if (abs(max_diff) < 0.5) {
      effect_size_desc <- "moderate"
    } else {
      effect_size_desc <- "large"
    }
    
    cat("This represents a ", effect_size_desc, " variation in treatment effect across levels of ", 
        modifier, ".\n\n", sep = "")
  }
}
```

```{r modifier-plots, eval=params$include_plots && length(credible_modifiers) > 0, fig.width=10, fig.height=6}
if (params$include_plots && length(credible_modifiers) > 0) {
  for (modifier in credible_modifiers) {
    cat("### Plot for ", modifier, "\n\n", sep = "")
    
    # Create plot
    if (method == "interaction") {
      # For standard interaction method, determine plot type
      model <- effect_results$models[[modifier]]
      data <- model$model$data
      
      if (is.numeric(data[[modifier]])) {
        # Continuous plot for numeric modifier
        print(plot_interaction(model, treatment_var, modifier, type = "continuous"))
      } else {
        # Categorical plot for factor modifier
        print(plot_interaction(model, treatment_var, modifier, type = "categorical"))
      }
    } else if (method == "bart" || method == "continuous") {
      # For BART or continuous methods, plot TE vs modifier value
      print(plot_interaction(effect_results, modifier, plot_type = "continuous"))
    }
  }
}
```

## Credibility Assessment

```{r credibility-assessment}
if (length(credible_modifiers) > 0) {
  # Assess credibility of subgroup effects
  credibility <- assess_subgroup_credibility(effect_results)
  
  if (nrow(credibility) > 0) {
    cat("### Credibility Criteria\n\n")
    
    # Display credibility assessment
    kable(credibility[, c("modifier", "criterion", "result")])
    
    cat("\n### Interpretation of Credibility\n\n")
    
    # Generate interpretation based on assessment
    cat("When evaluating subgroup effects, multiple criteria should be considered to assess credibility:\n\n")
    cat("1. **Pre-specification**: Effects that were hypothesized before analysis are more credible\n")
    cat("2. **Statistical interaction**: The strength of statistical evidence for the interaction\n")
    cat("3. **Consistency**: Whether the effect modification is consistent across studies\n")
    cat("4. **Multiple testing**: Adjustment for testing multiple potential modifiers\n")
    cat("5. **Biological plausibility**: Whether the modification has a plausible mechanism\n\n")
    
    cat("The subgroup effects identified in this analysis should be interpreted in light of these criteria. ")
    
    # Get score statistics if available
    if ("score" %in% names(credibility)) {
      scores <- credibility$score[!is.na(credibility$score)]
      if (length(scores) > 0) {
        mean_score <- mean(scores)
        
        if (mean_score > 0.8) {
          cat("Based on available criteria, these subgroup effects appear to have high credibility.")
        } else if (mean_score > 0.5) {
          cat("Based on available criteria, these subgroup effects have moderate credibility and warrant further investigation.")
        } else {
          cat("Based on available criteria, these subgroup effects have limited credibility and should be interpreted cautiously.")
        }
      }
    }
  }
}
```

## Clinical Implications

```{r clinical-implications}
if (length(credible_modifiers) > 0) {
  cat("The identified effect modifiers suggest potential personalization of treatment based on patient characteristics:\n\n")
  
  for (modifier in credible_modifiers) {
    mod_effect <- effect_results$modifier_effects[[modifier]]
    
    cat("- **", modifier, "**: ", sep = "")
    
    if (method == "interaction") {
      # Describe clinical implications for interaction method
      if (mod_effect$estimate > 0) {
        cat("Patients with higher values of ", modifier, " may benefit more from the treatment.\n", sep = "")
      } else {
        cat("Patients with lower values of ", modifier, " may benefit more from the treatment.\n", sep = "")
      }
    } else if (method == "bart" || method == "continuous") {
      # Describe clinical implications for BART or continuous methods
      te_data <- mod_effect$treatment_effects
      
      # Find modifier value with maximum effect
      max_idx <- which.max(te_data$mean_te)
      max_value <- te_data$modifier_value[max_idx]
      
      cat("Treatment effect is maximized at ", modifier, " = ", sprintf("%.2f", max_value), 
          ". Treatment decisions may be optimized by considering this value.\n", sep = "")
    }
  }
} else {
  cat("No credible effect modifiers were identified, suggesting that the treatment effect is relatively consistent across patient subgroups. This implies that treatment decisions may not need to be tailored based on the candidate variables examined.\n")
}
```

## Conclusion

```{r conclusion}
cat("This subgroup analysis aimed to identify patient characteristics that modify the effect of ", 
    treatment_var, " using the ", method, " method.\n\n", sep = "")

if (length(credible_modifiers) > 0) {
  cat("The analysis identified ", length(credible_modifiers), " credible effect modifier(s): ", 
      paste(credible_modifiers, collapse = ", "), ". ", sep = "")
  
  cat("Effect modification means that the treatment effect varies depending on the value of these variables. ",
      "This suggests potential heterogeneity in treatment effects across different subgroups.\n\n", sep = "")
  
  cat("For clinical practice, these findings suggest that treatment decisions could be tailored based on ",
      paste(credible_modifiers, collapse = ", "), " to optimize patient outcomes. ", sep = "")
  
  cat("However, these findings should be validated in further studies before being used to guide clinical decisions.")
} else {
  cat("The analysis did not identify any credible effect modifiers among the ", 
      length(candidate_modifiers), " candidate variables examined. ", sep = "")
  
  cat("This suggests that the treatment effect is relatively consistent across the range of variables examined. ",
      "The absence of effect modification implies a more uniform treatment effect across different patient subgroups.")
}
```

---

Report generated with `bayestrials` on `r format(params$timestamp, "%Y-%m-%d %H:%M:%S")`.
')
  } else {
    stop(paste("Unknown template:", template_name), call. = FALSE)
  }
}

#' Generate comprehensive sensitivity analysis reports
#'
#' Creates detailed reports documenting the robustness of Bayesian clinical trial 
#' results across different modeling assumptions. Essential for clinical documentation
#' to demonstrate that conclusions are not overly dependent on specific prior 
#' specifications, model choices, or analysis decisions.
#'
#' @param sensitivity_results A sensitivity_results object from 
#'   \code{\link{run_sensitivity_analyses}}. Must contain comparison results
#'   across different analysis scenarios:
#'   \itemize{
#'     \item \strong{Prior sensitivity:} Results from different prior specifications
#'     \item \strong{Model sensitivity:} Results from alternative model structures
#'     \item \strong{Missing data sensitivity:} Results from different imputation approaches
#'     \item \strong{Outlier sensitivity:} Results with/without influential observations
#'   }
#'   
#' @param output_file Character string specifying output file path (optional).
#'   If not provided, automatically generates timestamped filename:
#'   \code{"sensitivity_analysis_[YYYYMMDD_HHMMSS].[extension]"}
#'   
#' @param format Character string specifying output format:
#'   \itemize{
#'     \item \code{"html"} - Interactive report with embedded plots (default)
#'     \item \code{"pdf"} - Publication-ready document for clinical publications
#'     \item \code{"word"} - Editable document for collaborative review
#'   }
#'   
#' @param include_plots Logical indicating whether to include visualizations:
#'   \itemize{
#'     \item \code{TRUE} - Forest plots comparing effect estimates across scenarios
#'     \item \code{FALSE} - Text and tables only
#'   }
#'
#' @return Character string with path to generated sensitivity analysis report.
#'   Report includes automated assessment of result robustness based on 
#'   percentage differences between scenarios.
#'
#' @details
#' \strong{Report Content Overview:}
#' 
#' \strong{Sensitivity Assessment Criteria:}
#' \itemize{
#'   \item \strong{Robust (< 10% difference):} Results minimally affected by assumptions
#'   \item \strong{Moderate sensitivity (10-30%):} Some influence of assumptions
#'   \item \strong{High sensitivity (> 30%):} Strong dependence on assumptions
#' }
#' 
#' \strong{Clinical Interpretation Guidelines:}
#' \itemize{
#'   \item Results robust to prior specifications indicate strong data informativeness
#'   \item Model sensitivity highlights importance of structural assumptions
#'   \item Missing data sensitivity reveals potential bias from incomplete observations
#'   \item Consistent directions of effect across scenarios strengthen conclusions
#' }
#' 
#' \strong{Clinical Applications:}
#' \itemize{
#'   \item Sensitivity analyses recommended for key efficacy and safety endpoints
#'   \item Emphasize robustness to prior assumptions in clinical interpretation
#'   \item Document pre-specified sensitivity scenarios in statistical analysis plan
#'   \item Include discussion of clinical meaningfulness of observed variations
#' }
#'
#' @examples
#' \dontrun{
#' # Run comprehensive sensitivity analysis
#' data(synthetic_trial_data)
#' base_model_spec <- create_model_spec(
#'   outcome ~ treatment + age + sex,
#'   family = "gaussian"
#' )
#' 
#' # Define prior sensitivity scenarios
#' prior_sets <- list(
#'   neutral = default_priors("neutral"),
#'   skeptical = default_priors("skeptical"),
#'   optimistic = default_priors("optimistic")
#' )
#' 
#' # Run sensitivity analysis
#' sensitivity_results <- run_sensitivity_analyses(
#'   model_spec = base_model_spec,
#'   data = synthetic_trial_data,
#'   scenarios = list(
#'     prior_sensitivity = prior_sets,
#'     missing_data = c("complete_case", "imputation")
#'   )
#' )
#' 
#' # Generate comprehensive sensitivity report
#' sensitivity_report_path <- sensitivity_report(
#'   sensitivity_results,
#'   output_file = "efficacy_sensitivity_analysis.html",
#'   format = "html",
#'   include_plots = TRUE
#' )
#' 
#' # Open report in browser
#' browseURL(sensitivity_report_path)
#' 
#' # Generate PDF for formal documentation
#' formal_sensitivity <- sensitivity_report(
#'   sensitivity_results,
#'   output_file = "formal_sensitivity_analysis.pdf",
#'   format = "pdf"
#' )
#' 
#' # Generate simplified version for clinical team
#' clinical_sensitivity <- sensitivity_report(
#'   sensitivity_results,
#'   output_file = "clinical_sensitivity_summary.word",
#'   format = "word",
#'   include_plots = FALSE  # Tables only for easier editing
#' )
#' }
#'
#' @seealso 
#' \code{\link{run_sensitivity_analyses}} for conducting sensitivity analyses,
#' \code{\link{generate_report}} for standard analysis reports,
#' \code{\link{compare_models}} for model comparison utilities,
#' \code{\link{plot_sensitivity}} for sensitivity visualization
#'
#' @export
sensitivity_report <- function(sensitivity_results,
                              output_file = NULL,
                              format = c("html", "pdf", "word"),
                              include_plots = TRUE) {
  
  # Check sensitivity_results
  if (!inherits(sensitivity_results, "sensitivity_results")) {
    stop("sensitivity_results must be the result of run_sensitivity_analyses", call. = FALSE)
  }
  
  # Use generate_report with the sensitivity_report template
  report_path <- generate_report(
    model = sensitivity_results,
    template = "sensitivity_report",
    output_file = output_file,
    format = format,
    include_plots = include_plots,
    include_diagnostics = FALSE
  )
  
  return(report_path)
}

#' Generate comprehensive subgroup analysis reports
#'
#' Creates detailed reports for effect modifier analysis and subgroup investigations
#' in clinical trials. Essential for personalized medicine applications and 
#' clinical publications exploring differential treatment effects across patient
#' populations. Includes statistical testing, clinical interpretation, and
#' credibility assessments.
#'
#' @param effect_results An effect_modifier_results object from 
#'   \code{\link{identify_effect_modifiers}} containing:
#'   \itemize{
#'     \item Treatment effect estimates by subgroup
#'     \item Interaction test results and p-values
#'     \item Effect modifier credibility assessments
#'     \item Subgroup sample sizes and clinical characteristics
#'   }
#'   Must include results from systematic effect modifier analysis.
#'   
#' @param output_file Character string specifying output file path (optional).
#'   If not provided, automatically generates timestamped filename:
#'   \code{"subgroup_analysis_[YYYYMMDD_HHMMSS].[extension]"}
#'   
#' @param format Character string specifying output format:
#'   \itemize{
#'     \item \code{"html"} - Interactive report with embedded forest plots (default)
#'     \item \code{"pdf"} - Publication-ready document with statistical tables
#'     \item \code{"word"} - Editable document for collaborative review and clinical comments
#'   }
#'   
#' @param include_plots Logical indicating whether to include visualizations:
#'   \itemize{
#'     \item \code{TRUE} - Forest plots, interaction plots, and effect modifier visualizations
#'     \item \code{FALSE} - Statistical tables and text summaries only
#'   }
#'
#' @return Character string with path to generated subgroup analysis report.
#'   Report includes automated credibility assessments and statistically-compliant
#'   interpretation guidelines.
#'
#' @details
#' \strong{Subgroup Analysis Reporting Standards:}
#' 
#' \strong{Credibility Assessment Framework:}
#' \itemize{
#'   \item \strong{Pre-specification:} Were subgroups defined before data analysis?
#'   \item \strong{Statistical significance:} Are interaction tests statistically significant?
#'   \item \strong{Clinical plausibility:} Is the biological mechanism reasonable?
#'   \item \strong{Consistency:} Are findings consistent across related endpoints?
#'   \item \strong{Effect magnitude:} Is the difference clinically meaningful?
#' }
#' 
#' \strong{Statistical Guidelines:}
#' \itemize{
#'   \item Subgroup analyses should be pre-specified with clinical rationale
#'   \item Multiple testing adjustments required for exploratory analyses
#'   \item Distinguish between confirmatory and exploratory subgroup analyses
#'   \item Report methodology and multiplicity considerations
#' }
#' 
#' \strong{Clinical Interpretation Guidelines:}
#' \itemize{
#'   \item \strong{Qualitative interactions:} Treatment beneficial in some, harmful in others
#'   \item \strong{Quantitative interactions:} Treatment beneficial in all, but magnitude varies
#'   \item \strong{Statistical vs. clinical significance:} Consider both p-values and effect sizes
#'   \item \strong{External validity:} Assess generalizability to broader patient populations
#' }
#' 
#' \strong{Report Content Structure:}
#' \itemize{
#'   \item Executive summary with key findings
#'   \item Statistical methodology and multiple testing approach
#'   \item Subgroup characteristics and baseline comparisons
#'   \item Forest plots with confidence intervals
#'   \item Formal interaction tests and p-values
#'   \item Credibility assessment using established criteria
#'   \item Clinical interpretation and recommendations
#' }
#'
#' @examples
#' \dontrun{
#' # Conduct comprehensive effect modifier analysis
#' data(synthetic_trial_data)
#' model_spec <- create_model_spec(
#'   outcome ~ treatment * age + treatment * sex + baseline_score,
#'   family = "gaussian"
#' )
#' 
#' priors <- default_priors()
#' fitted_model <- fit_model(model_spec, synthetic_trial_data, priors)
#' 
#' # Identify potential effect modifiers
#' effect_modifier_results <- identify_effect_modifiers(
#'   model = fitted_model,
#'   data = synthetic_trial_data,
#'   treatment_var = "treatment",
#'   modifier_vars = c("age", "sex", "baseline_severity")
#' )
#' 
#' # Generate comprehensive subgroup analysis report
#' subgroup_report_path <- subgroup_report(
#'   effect_modifier_results,
#'   output_file = "subgroup_efficacy_analysis.html",
#'   format = "html",
#'   include_plots = TRUE
#' )
#' 
#' # View report
#' browseURL(subgroup_report_path)
#' 
#' # Generate PDF for formal documentation
#' formal_subgroup <- subgroup_report(
#'   effect_modifier_results,
#'   output_file = "formal_subgroup_analysis.pdf",
#'   format = "pdf"
#' )
#' 
#' # Pre-specified vs. exploratory subgroup analyses
#' prespecified_modifiers <- c("age_group", "sex")  # Pre-specified in SAP
#' exploratory_modifiers <- c("biomarker_level", "comorbidity_score")  # Exploratory
#' 
#' # Run separate analyses
#' prespec_results <- identify_effect_modifiers(
#'   fitted_model, data, "treatment", prespecified_modifiers,
#'   analysis_type = "confirmatory"
#' )
#' 
#' exploratory_results <- identify_effect_modifiers(
#'   fitted_model, data, "treatment", exploratory_modifiers,
#'   analysis_type = "exploratory",
#'   adjust_p_values = TRUE  # Multiple testing correction
#' )
#' 
#' # Generate separate reports
#' subgroup_report(prespec_results, "prespecified_subgroups.html")
#' subgroup_report(exploratory_results, "exploratory_subgroups.html")
#' }
#'
#' @seealso 
#' \code{\link{identify_effect_modifiers}} for conducting subgroup analyses,
#' \code{\link{assess_subgroup_credibility}} for credibility assessment,
#' \code{\link{plot_interaction}} for interaction visualizations,
#' \code{\link{plot_forest}} for subgroup forest plots
#'
#' @export
subgroup_report <- function(effect_results,
                           output_file = NULL,
                           format = c("html", "pdf", "word"),
                           include_plots = TRUE) {
  
  # Check effect_results
  if (!inherits(effect_results, "effect_modifier_results")) {
    stop("effect_results must be the result of identify_effect_modifiers", call. = FALSE)
  }
  
  # Use generate_report with the subgroup_report template
  report_path <- generate_report(
    model = effect_results,
    template = "subgroup_report",
    output_file = output_file,
    format = format,
    include_plots = include_plots,
    include_diagnostics = FALSE
  )
  
  return(report_path)
}

#' Generate reproducibility and audit trail reports
#'
#' Creates detailed documentation of analysis workflow, computational environment,
#' and reproducibility information for clinical trial analyses. Essential for 
#' clinical documentation, audit trails, and ensuring analysis transparency.
#' Documents all analysis steps, package versions, and computational parameters.
#'
#' @param workflow AnalysisWorkflow object from \code{\link{capture_workflow}}
#'   or character string with path to saved workflow file. Contains:
#'   \itemize{
#'     \item Complete analysis steps with function calls and parameters
#'     \item Package versions and computational environment details
#'     \item Data checksums and processing steps
#'     \item Model specifications and prior assumptions
#'     \item Results and intermediate outputs
#'   }
#'   
#' @param format Character string specifying output format:
#'   \itemize{
#'     \item \code{"html"} - Self-contained report with embedded metadata (default)
#'     \item \code{"pdf"} - Formal document suitable for clinical archives
#'     \item \code{"word"} - Editable format for collaborative review
#'   }
#'   
#' @param verify_packages Logical indicating whether to verify current package
#'   versions against workflow requirements:
#'   \itemize{
#'     \item \code{TRUE} - Checks version compatibility and flags differences (recommended)
#'     \item \code{FALSE} - Documents workflow without verification
#'   }
#'
#' @return Character string with path to generated reproducibility report.
#'   Report serves as complete audit trail for quality assurance.
#'
#' @details
#' \strong{Reproducibility Documentation Standards:}
#' 
#' \strong{Computational Environment:}
#' \itemize{
#'   \item R version and platform information
#'   \item Complete package versions with checksums
#'   \item Random seed documentation for MCMC sampling
#'   \item Hardware specifications and execution time
#' }
#' 
#' \strong{Analysis Workflow:}
#' \itemize{
#'   \item Step-by-step function calls with all parameters
#'   \item Data preprocessing and transformation steps
#'   \item Model specification and prior assumptions
#'   \item Convergence diagnostics and model validation
#' }
#' 
#' \strong{Quality Assurance:}
#' \itemize{
#'   \item Data integrity checks and validation
#'   \item Model convergence assessments
#'   \item Results consistency verification
#'   \item Error handling and exception documentation
#' }
#' 
#' \strong{Documentation Standards:}
#' \itemize{
#'   \item Electronic records compliance for clinical trials
#'   \item Statistical principles documentation
#'   \item Software validation and assurance guidelines
#'   \item Audit trail requirements for quality assurance
#' }
#'
#' @examples
#' \dontrun{
#' # Initialize reproducibility tracking
#' workflow <- initialize_reproducibility("clinical_trial_analysis")
#' 
#' # Capture analysis steps
#' workflow <- add_step(workflow, "data_prep", "load_clinical_trial", 
#'                     list(file = "trial_data.csv"))
#' workflow <- add_step(workflow, "model_spec", "create_model_spec",
#'                     list(formula = "outcome ~ treatment + age"))
#' workflow <- add_step(workflow, "model_fit", "fit_model",
#'                     list(chains = 4, iter = 2000))
#' 
#' # Generate reproducibility report
#' repro_report <- reproducibility_report(
#'   workflow,
#'   format = "html",
#'   verify_packages = TRUE
#' )
#' 
#' # View reproducibility documentation
#' browseURL(repro_report)
#' 
#' # Generate audit trail for formal documentation
#' audit_trail <- reproducibility_report(
#'   workflow,
#'   format = "pdf",
#'   verify_packages = TRUE
#' )
#' 
#' # Load and document previous analysis
#' saved_workflow <- load_workflow("previous_analysis.json")
#' historical_report <- reproducibility_report(
#'   saved_workflow,
#'   format = "word",
#'   verify_packages = FALSE  # Historical analysis
#' )
#' 
#' # Comprehensive project documentation
#' project_workflow <- capture_workflow(
#'   analysis_steps = list(
#'     data_processing = data_processing_steps,
#'     model_fitting = model_fitting_steps,
#'     sensitivity_analysis = sensitivity_steps
#'   ),
#'   save_path = "complete_analysis_workflow.json"
#' )
#' 
#' final_documentation <- reproducibility_report(
#'   project_workflow,
#'   format = "pdf"
#' )
#' }
#'
#' @seealso 
#' \code{\link{initialize_reproducibility}} for starting workflow tracking,
#' \code{\link{capture_workflow}} for documenting analysis steps,
#' \code{\link{save_workflow}} and \code{\link{load_workflow}} for workflow persistence,
#' \code{\link{restore_analysis}} for reproducing previous analyses
#'
#' @export
reproducibility_report <- function(workflow,
                                 format = c("html", "pdf", "word"),
                                 verify_packages = TRUE) {
  
  format <- match.arg(format)
  
  # Check if rmarkdown is available
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rmarkdown package is required but not available", call. = FALSE)
  }
  
  # Check workflow
  if (is.character(workflow)) {
    # Load workflow from file
    if (!file.exists(workflow)) {
      stop(paste("Workflow file not found:", workflow), call. = FALSE)
    }
    
    workflow_data <- jsonlite::fromJSON(workflow)
  } else if (inherits(workflow, "AnalysisWorkflow")) {
    # Use workflow object directly
    workflow_data <- workflow
  } else {
    stop("workflow must be a file path or AnalysisWorkflow object", call. = FALSE)
  }
  
  # Create content for reproducibility report
  reproducibility_content <- paste0(
    "---\n",
    "title: \"Reproducibility Report\"\n",
    "author: \"bayestrials\"\n",
    "date: \"", format(Sys.time(), "%Y-%m-%d"), "\"\n",
    "output: ", format, "_document\n",
    "---\n\n",
    
    "## Analysis Workflow Overview\n\n",
    
    "This document describes the reproducibility information for an analysis conducted with bayestrials.\n\n",
    
    "## Session Information\n\n",
    
    "R version: ", R.version$version.string, "\n\n",
    
    "Platform: ", R.version$platform, "\n\n",
    
    "bayestrials version: ", get_version(), "\n\n",
    
    "### Package Versions\n\n",
    
    "```{r echo=FALSE}\n",
    "sessionInfo()\n",
    "```\n\n",
    
    "## Analysis Steps\n\n"
  )
  
  # Add workflow steps
  if (!is.null(workflow_data$steps)) {
    for (i in seq_along(workflow_data$steps)) {
      step <- workflow_data$steps[[i]]
      
      reproducibility_content <- paste0(
        reproducibility_content,
        "### Step ", i, ": ", step$name, "\n\n",
        "Function: `", step$function_name, "`\n\n",
        "Arguments:\n\n",
        "```r\n",
        paste(names(step$arguments), "=", sapply(step$arguments, function(x) paste(deparse(x), collapse = "\n")), 
              collapse = ",\n"),
        "\n```\n\n"
      )
    }
  } else {
    reproducibility_content <- paste0(
      reproducibility_content,
      "No analysis steps were recorded in the workflow.\n\n"
    )
  }
  
  # Add verification information
  if (verify_packages && !is.null(workflow_data$packages)) {
    reproducibility_content <- paste0(
      reproducibility_content,
      "## Package Verification\n\n",
      
      "The following table compares the package versions used in the original analysis with those currently installed:\n\n",
      
      "```{r echo=FALSE}\n",
      "original_packages <- list(\n"
    )
    
    for (i in seq_along(workflow_data$packages)) {
      pkg <- workflow_data$packages[[i]]
      reproducibility_content <- paste0(
        reproducibility_content,
        "  ", pkg$name, " = \"", pkg$version, "\"",
        if (i < length(workflow_data$packages)) ",\n" else "\n"
      )
    }
    
    reproducibility_content <- paste0(
      reproducibility_content,
      ")\n\n",
      
      "current_packages <- lapply(names(original_packages), function(pkg) {\n",
      "  if (requireNamespace(pkg, quietly = TRUE)) {\n",
      "    as.character(packageVersion(pkg))\n",
      "  } else {\n",
      "    \"not installed\"\n",
      "  }\n",
      "})\n",
      "names(current_packages) <- names(original_packages)\n\n",
      
      "package_df <- data.frame(\n",
      "  Package = names(original_packages),\n",
      "  Original = unlist(original_packages),\n",
      "  Current = unlist(current_packages),\n",
      "  Match = unlist(original_packages) == unlist(current_packages),\n",
      "  stringsAsFactors = FALSE\n",
      ")\n\n",
      
      "knitr::kable(package_df)\n",
      "```\n\n"
    )
  }
  
  # Add conclusion
  reproducibility_content <- paste0(
    reproducibility_content,
    "## Conclusion\n\n",
    
    "This report provides the information necessary to reproduce the analysis. ",
    "To ensure complete reproducibility, ensure that the same package versions are installed ",
    "and that the analysis steps are followed in the same order.\n\n",
    
    "---\n\n",
    "Generated with `bayestrials` on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ".\n"
  )
  
  # Create temporary Rmd file
  rmd_file <- tempfile(pattern = "reproducibility_", fileext = ".Rmd")
  writeLines(reproducibility_content, rmd_file)
  
  # Set default output file if not provided
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  output_file <- paste0("reproducibility_report_", timestamp, ".", 
                       switch(format, html = "html", pdf = "pdf", word = "docx"))
  
  # Render report
  output_path <- rmarkdown::render(
    input = rmd_file,
    output_file = output_file,
    output_format = switch(format,
                          html = rmarkdown::html_document(),
                          pdf = rmarkdown::pdf_document(),
                          word = rmarkdown::word_document()),
    quiet = TRUE
  )
  
  return(output_path)
}

#' AnalysisWorkflow class
#' 
#' @description S3 class for tracking analysis steps and reproducibility information
#' 
#' @export
AnalysisWorkflow <- function() {
  structure(
    list(
      steps = list(),
      packages = list(),
      timestamp = Sys.time(),
      session_info = get_session_info()
    ),
    class = "AnalysisWorkflow"
  )
}

#' Add a step to an AnalysisWorkflow
#' 
#' @param workflow An AnalysisWorkflow object
#' @param name Name of the step
#' @param function_name Name of the function called
#' @param args Arguments passed to the function
#' @return Updated AnalysisWorkflow object
#' @export
add_step <- function(workflow, name, function_name, args) {
  if (!inherits(workflow, "AnalysisWorkflow")) {
    stop("workflow must be an AnalysisWorkflow object", call. = FALSE)
  }
  
  # Create step
  step <- list(
    name = name,
    function_name = function_name,
    arguments = args,
    timestamp = Sys.time()
  )
  
  # Add to workflow
  workflow$steps[[length(workflow$steps) + 1]] <- step
  
  # Check if we need to update package information
  if (length(workflow$packages) == 0) {
    # Get loaded packages
    loaded_pkgs <- .packages()
    
    # Add package versions
    workflow$packages <- lapply(loaded_pkgs, function(pkg) {
      list(
        name = pkg,
        version = as.character(packageVersion(pkg))
      )
    })
    
    # Set names
    names(workflow$packages) <- loaded_pkgs
  }
  
  return(workflow)
}

#' Save AnalysisWorkflow to a file
#' 
#' @param workflow An AnalysisWorkflow object
#' @param file_path Path to save the workflow
#' @return Invisibly returns TRUE if successful
#' @export
save_workflow <- function(workflow, file_path) {
  if (!inherits(workflow, "AnalysisWorkflow")) {
    stop("workflow must be an AnalysisWorkflow object", call. = FALSE)
  }
  
  # Check if jsonlite is available
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package is required but not available", call. = FALSE)
  }
  
  # Convert to JSON
  json_data <- jsonlite::toJSON(workflow, pretty = TRUE, auto_unbox = TRUE)
  
  # Write to file
  writeLines(json_data, file_path)
  
  invisible(TRUE)
}

#' Load AnalysisWorkflow from a file
#' 
#' @param file_path Path to the workflow file
#' @return An AnalysisWorkflow object
#' @export
load_workflow <- function(file_path) {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path), call. = FALSE)
  }
  
  # Check if jsonlite is available
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package is required but not available", call. = FALSE)
  }
  
  # Read from file
  json_data <- readLines(file_path, warn = FALSE)
  
  # Parse JSON
  workflow <- jsonlite::fromJSON(json_data, simplifyVector = FALSE)
  
  # Set class
  class(workflow) <- "AnalysisWorkflow"
  
  return(workflow)
}

#' Initialize reproducibility settings
#' 
#' @param project_path Path to the project directory
#' @param seed Random seed for reproducibility
#' @param log_level Level of logging detail
#' @return An AnalysisWorkflow object
#' @export
initialize_reproducibility <- function(project_path = ".",
                                    seed = NULL,
                                    log_level = "info") {
  
  # Set seed
  if (!is.null(seed)) {
    set_seed(seed)
  } else {
    seed <- as.integer(Sys.time())
    set_seed(seed)
  }
  
  # Set log level
  set_option("log_level", log_level)
  
  # Create workflow object
  workflow <- AnalysisWorkflow()
  
  # Add initialization step
  workflow <- add_step(
    workflow,
    name = "Initialization",
    function_name = "initialize_reproducibility",
    args = list(
      project_path = project_path,
      seed = seed,
      log_level = log_level
    )
  )
  
  # Try to set up renv if available
  if (requireNamespace("renv", quietly = TRUE)) {
    # Check if renv is already initialized
    if (!file.exists(file.path(project_path, "renv.lock"))) {
      # Initialize renv
      renv::init(project = project_path, restart = FALSE)
    }
    
    # Take snapshot
    renv::snapshot(project = project_path, prompt = FALSE)
    
    # Add renv info to workflow
    workflow$renv <- list(
      status = "initialized",
      lockfile = file.path(project_path, "renv.lock")
    )
  }
  
  # Save workflow
  workflow_path <- file.path(project_path, "bayestrials_workflow.json")
  save_workflow(workflow, workflow_path)
  
  # Log initialization
  log_message(paste("Reproducibility initialized with seed", seed), level = "info")
  log_message(paste("Workflow saved to", workflow_path), level = "info")
  
  return(workflow)
}

#' Capture workflow steps
#' 
#' @param analysis_steps List of analysis step specifications
#' @param save_path Path to save the workflow file
#' @return An AnalysisWorkflow object
#' @export
capture_workflow <- function(analysis_steps, save_path = "workflow.json") {
  # Create workflow object
  workflow <- AnalysisWorkflow()
  
  # Add steps
  for (step in analysis_steps) {
    workflow <- add_step(
      workflow,
      name = step$name,
      function_name = step$function_name,
      args = step$args
    )
  }
  
  # Save workflow
  save_workflow(workflow, save_path)
  
  return(workflow)
}

#' Restore analysis from workflow
#' 
#' @param workflow_path Path to workflow file
#' @param verify_packages Logical indicating whether to verify package versions
#' @return Results of the last step in the workflow
#' @export
restore_analysis <- function(workflow_path, verify_packages = TRUE) {
  # Load workflow
  workflow <- load_workflow(workflow_path)
  
  # Verify packages if requested
  if (verify_packages && !is.null(workflow$packages)) {
    pkg_issues <- character(0)
    
    for (pkg_name in names(workflow$packages)) {
      pkg <- workflow$packages[[pkg_name]]
      
      if (!requireNamespace(pkg$name, quietly = TRUE)) {
        pkg_issues <- c(pkg_issues, paste(pkg$name, "is not installed"))
      } else if (as.character(packageVersion(pkg$name)) != pkg$version) {
        pkg_issues <- c(pkg_issues, 
                       paste(pkg$name, "version mismatch:", 
                            as.character(packageVersion(pkg$name)), "vs", pkg$version))
      }
    }
    
    if (length(pkg_issues) > 0) {
      warning(paste("Package version issues detected:", 
                   paste(pkg_issues, collapse = "; "), 
                   "Results may not be reproducible"), call. = FALSE)
    }
  }
  
  # Restore renv if available
  if (!is.null(workflow$renv) && requireNamespace("renv", quietly = TRUE)) {
    if (file.exists(workflow$renv$lockfile)) {
      renv::restore(lockfile = workflow$renv$lockfile)
    }
  }
  
  # Execute steps
  results <- NULL
  
  for (i in seq_along(workflow$steps)) {
    step <- workflow$steps[[i]]
    
    log_message(paste("Executing step", i, ":", step$name), level = "info")
    
    # Get function
    func <- try(get(step$function_name), silent = TRUE)
    
    if (inherits(func, "try-error")) {
      stop(paste("Function not found:", step$function_name), call. = FALSE)
    }
    
    # Execute function with arguments
    results <- try(do.call(func, step$arguments), silent = TRUE)
    
    if (inherits(results, "try-error")) {
      warning(paste("Error executing step", i, ":", step$name, "-", attr(results, "condition")), 
             call. = FALSE)
      break
    }
  }
  
  return(results)
}