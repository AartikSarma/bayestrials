#' Visualize posterior distributions from Bayesian clinical trial models
#'
#' Creates publication-ready plots of posterior distributions for model parameters.
#' Essential for communicating Bayesian results to clinical and statistical audiences,
#' showing both point estimates and uncertainty quantification. Multiple plot types
#' provide different perspectives on the posterior distributions.
#'
#' @param posterior_data Data frame, matrix, or BayesianModel object containing 
#'   posterior samples. Can be:
#'   \itemize{
#'     \item Data frame in long format with \code{parameter} and \code{value} columns
#'     \item Data frame in wide format with parameters as columns, samples as rows
#'     \item Matrix with parameters as columns, samples as rows
#'     \item BayesianModel object from \code{\link{fit_model}}
#'   }
#'   
#' @param plot_type Character string specifying visualization type:
#'   \itemize{
#'     \item \code{"halfeye"} - Shows full posterior distribution with quantile intervals.
#'       Best for showing complete uncertainty profile. Requires ggdist package.
#'     \item \code{"interval"} - Point estimates with credible intervals (50% and 95%).
#'       Ideal for formal documentation and clinical summaries.
#'     \item \code{"dots"} - Quantile dot plots showing posterior mass distribution.
#'       Good for discrete or multi-modal posteriors. Requires ggdist package.
#'   }
#'   
#' @param facet_by Character string naming a grouping variable for creating 
#'   separate plot panels. Useful for:
#'   \itemize{
#'     \item Different trials in meta-analysis: \code{facet_by = "trial"}
#'     \item Parameter types: \code{facet_by = "parameter_type"}
#'     \item Subgroups: \code{facet_by = "subgroup"}
#'   }
#'   Variable must be present in posterior_data.
#'   
#' @param color_by Character string naming a variable for color coding.
#'   Enables comparison across:
#'   \itemize{
#'     \item Prior specifications: \code{color_by = "prior_type"}
#'     \item Model types: \code{color_by = "model_name"}
#'     \item Sensitivity analyses: \code{color_by = "analysis_type"}
#'   }
#'   Variable must be present in posterior_data.
#'   
#' @param theme Character string specifying plot appearance:
#'   \itemize{
#'     \item \code{"publication"} - Clean theme suitable for journals and reports
#'     \item \code{"minimal"} - Simple theme with minimal visual elements
#'     \item \code{"classic"} - Traditional statistical plot appearance
#'     \item \code{"bw"} - Black and white theme for print publications
#'   }
#'
#' @return A ggplot object that can be further customized or saved. The plot 
#'   includes a reference line at zero for easy interpretation of treatment effects.
#'
#' @details
#' \strong{Clinical Interpretation:}
#' 
#' \strong{For Treatment Effects:}
#' \itemize{
#'   \item Distributions centered around zero suggest no treatment effect
#'   \item Distributions clearly above/below zero suggest beneficial/harmful effects
#'   \item Width of distribution reflects uncertainty in effect size
#'   \item Tail probability beyond clinical significance threshold indicates strength of evidence
#' }
#' 
#' \strong{For Safety Parameters:}
#' \itemize{
#'   \item Narrow distributions suggest predictable safety profile
#'   \item Long tails indicate potential for extreme adverse events
#'   \item Compare across treatment groups for differential safety assessment
#' }
#' 
#' \strong{Best Practices:}
#' \itemize{
#'   \item Include credible intervals (not confidence intervals) for Bayesian results
#'   \item Show posterior probability of clinically meaningful effect
#'   \item Document prior sensitivity through color coding different prior sets
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data and fit model
#' data(synthetic_trial_data)
#' priors <- default_priors()
#' model_spec <- create_model_spec(outcome ~ treatment + age, family = "gaussian")
#' model <- fit_model(model_spec, synthetic_trial_data, priors)
#' 
#' # Basic posterior plot
#' plot_posterior(model, plot_type = "halfeye")
#' 
#' # Interval plot for formal documentation
#' p1 <- plot_posterior(model, plot_type = "interval", theme = "publication")
#' print(p1)
#' 
#' # Compare across different prior specifications
#' # (assuming you fit models with different priors)
#' posteriors_combined <- rbind(
#'   data.frame(extract_posterior(model_neutral), prior_type = "neutral"),
#'   data.frame(extract_posterior(model_skeptical), prior_type = "skeptical")
#' )
#' 
#' # Color by prior type
#' plot_posterior(posteriors_combined, 
#'               plot_type = "halfeye", 
#'               color_by = "prior_type")
#' 
#' # Facet by parameter type for complex models
#' plot_posterior(model, 
#'               plot_type = "interval",
#'               facet_by = "parameter")
#' 
#' # Extract and format posteriors manually
#' post_samples <- extract_posterior(model)
#' post_long <- tidyr::pivot_longer(post_samples, 
#'                                 cols = everything(),
#'                                 names_to = "parameter", 
#'                                 values_to = "value")
#' 
#' # Create publication-ready plot
#' final_plot <- plot_posterior(post_long, 
#'                             plot_type = "interval",
#'                             theme = "publication") +
#'   ggplot2::labs(
#'     title = "Treatment Effects in Clinical Trial",
#'     subtitle = "Posterior distributions with 95% credible intervals",
#'     caption = "Bayesian analysis with neutral priors"
#'   )
#' }
#'
#' @seealso 
#' \code{\link{fit_model}} for fitting Bayesian models,
#' \code{\link{extract_posterior}} for extracting posterior samples,
#' \code{\link{plot_forest}} for comparing effects across trials,
#' \code{\link{plot_model_comparison}} for comparing different models
#'
#' @export
plot_posterior <- function(posterior_data,
                          plot_type = c("halfeye", "interval", "dots"),
                          facet_by = NULL,
                          color_by = NULL,
                          theme = "publication") {
  
  plot_type <- match.arg(plot_type)
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Convert to data frame if matrix
  if (is.matrix(posterior_data)) {
    posterior_data <- as.data.frame(posterior_data)
  }
  
  # Convert to long format if not already
  if (!"parameter" %in% names(posterior_data) || !"value" %in% names(posterior_data)) {
    # Check if it's a BayesianModel object
    if (inherits(posterior_data, "BayesianModel")) {
      if (!is.null(posterior_data$model)) {
        posterior_data <- brms::posterior_samples(posterior_data$model)
      } else {
        stop("BayesianModel has no fitted model", call. = FALSE)
      }
    }
    
    # Convert wide to long format
    posterior_long <- tidyr::pivot_longer(
      posterior_data,
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    )
  } else {
    # Already in long format
    posterior_long <- posterior_data
  }
  
  # Set up base plot
  p <- ggplot2::ggplot(posterior_long, ggplot2::aes(x = value, y = parameter))
  
  # Add layers based on plot type
  if (plot_type == "halfeye") {
    # Check for ggdist package
    if (!requireNamespace("ggdist", quietly = TRUE)) {
      stop("ggdist package is required for halfeye plots but not available", call. = FALSE)
    }
    
    # Create halfeye plot
    if (!is.null(color_by)) {
      p <- p + ggdist::stat_halfeye(ggplot2::aes(fill = .data[[color_by]]), 
                                  alpha = 0.8, normalize = "groups")
    } else {
      p <- p + ggdist::stat_halfeye(fill = "skyblue", alpha = 0.8)
    }
    
  } else if (plot_type == "interval") {
    # Create interval plot
    if (!is.null(color_by)) {
      p <- p + ggplot2::stat_summary(
        ggplot2::aes(color = .data[[color_by]]),
        fun.data = function(x) {
          quantiles <- stats::quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
          return(list(
            y = mean(x),
            ymin = quantiles[1],
            ymax = quantiles[5],
            lower = quantiles[2],
            upper = quantiles[4]
          ))
        },
        geom = "pointrange"
      )
    } else {
      p <- p + ggplot2::stat_summary(
        fun.data = function(x) {
          quantiles <- stats::quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
          return(list(
            y = mean(x),
            ymin = quantiles[1],
            ymax = quantiles[5],
            lower = quantiles[2],
            upper = quantiles[4]
          ))
        },
        geom = "pointrange",
        color = "skyblue"
      )
    }
    
  } else if (plot_type == "dots") {
    # Check for ggdist package
    if (!requireNamespace("ggdist", quietly = TRUE)) {
      stop("ggdist package is required for dots plots but not available", call. = FALSE)
    }
    
    # Create dots plot
    if (!is.null(color_by)) {
      p <- p + ggdist::stat_dotsinterval(ggplot2::aes(fill = .data[[color_by]]), 
                                       alpha = 0.8, normalize = "groups")
    } else {
      p <- p + ggdist::stat_dotsinterval(fill = "skyblue", alpha = 0.8)
    }
  }
  
  # Add reference line at 0
  p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")
  
  # Add faceting if requested
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_by]]), scales = "free")
  }
  
  # Apply theme
  if (theme == "publication") {
    p <- p + ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, color = "gray80"),
        strip.background = ggplot2::element_rect(fill = "gray95", color = NA),
        strip.text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )
  } else if (theme == "minimal") {
    p <- p + ggplot2::theme_minimal()
  } else if (theme == "classic") {
    p <- p + ggplot2::theme_classic()
  } else if (theme == "bw") {
    p <- p + ggplot2::theme_bw()
  }
  
  # Add labels
  p <- p + ggplot2::labs(
    title = "Posterior Distributions",
    x = "Parameter Value",
    y = NULL
  )
  
  return(p)
}

#' Create forest plots for clinical trial effect comparisons
#'
#' Generates forest plots for comparing treatment effects across multiple trials,
#' subgroups, or analyses. Essential for meta-analyses, subgroup analyses, and 
#' clinical publications. Shows individual study effects with confidence intervals
#' and optional overall summary effect.
#'
#' @param effect_data Data frame containing effect estimates and confidence intervals.
#'   Must include columns for trial/study identifier and effect estimates. 
#'   Flexible column naming - function will automatically detect:
#'   \itemize{
#'     \item \strong{Effect estimates:} \code{estimate}, \code{mean}, \code{effect}, \code{effect_size}
#'     \item \strong{Lower CI:} \code{ci_lower}, \code{lower}, \code{lower_ci}, \code{lb}, \code{q2.5}
#'     \item \strong{Upper CI:} \code{ci_upper}, \code{upper}, \code{upper_ci}, \code{ub}, \code{q97.5}
#'     \item \strong{Trial ID:} \code{trial}, \code{group}, \code{study}, \code{subgroup}
#'   }
#'   
#' @param parameter Character string naming the parameter being plotted.
#'   Used in plot title and for interpretation. Examples:
#'   \itemize{
#'     \item \code{"treatment_effect"} - Main treatment comparison
#'     \item \code{"hazard_ratio"} - Survival analysis results
#'     \item \code{"odds_ratio"} - Logistic regression effects
#'     \item \code{"mean_difference"} - Continuous outcome differences
#'   }
#'   
#' @param order_by Character string specifying how to order studies in the plot:
#'   \itemize{
#'     \item \code{"effect_size"} - Order by magnitude of effect (default)
#'     \item Any column name in \code{effect_data} - Order by that variable
#'     \item Studies are ordered from smallest to largest effect
#'   }
#'   
#' @param xlab Character string for x-axis label. Should describe the 
#'   parameter scale and direction. Examples:
#'   \itemize{
#'     \item \code{"Treatment Effect (Cohen's d)"}
#'     \item \code{"Hazard Ratio (log scale)"}
#'     \item \code{"Mean Difference in Primary Endpoint"}
#'   }
#'   
#' @param show_summary Logical indicating whether to calculate and display 
#'   an overall summary effect. Uses inverse-variance weighting:
#'   \itemize{
#'     \item \code{TRUE} - Shows "Overall" row with meta-analytic summary
#'     \item \code{FALSE} - Individual studies only
#'   }
#'   Summary calculation assumes normal approximation for confidence intervals.
#'
#' @return A ggplot object showing the forest plot. Includes:
#'   \itemize{
#'     \item Point estimates as squares (larger for summary effect)
#'     \item Horizontal lines for confidence intervals  
#'     \item Vertical reference line at zero (no effect)
#'     \item Numerical values displayed if ggtext package available
#'   }
#'
#' @details
#' \strong{Clinical Interpretation Guidelines:}
#' 
#' \strong{Effect Direction:}
#' \itemize{
#'   \item Points to the right of zero line indicate beneficial treatment effects
#'   \item Points to the left indicate harmful effects
#'   \item Confidence intervals crossing zero suggest non-significant effects
#' }
#' 
#' \strong{Heterogeneity Assessment:}
#' \itemize{
#'   \item Overlapping confidence intervals suggest consistent effects
#'   \item Non-overlapping intervals indicate significant heterogeneity
#'   \item Wide spread of point estimates suggests population differences
#' }
#' 
#' \strong{Summary Effect Interpretation:}
#' \itemize{
#'   \item Overall effect summarizes evidence across all studies
#'   \item Larger square indicates more precise (narrower CI) summary
#'   \item Use with caution if substantial heterogeneity observed
#' }
#' 
#' \strong{Best Practices:}
#' \itemize{
#'   \item Include sample sizes and study characteristics in caption
#'   \item Document fixed vs. random effects meta-analysis approach
#'   \item Consider forest plots for key subgroups and sensitivity analyses
#' }
#'
#' @examples
#' \dontrun{
#' # Prepare effect data for multiple trials
#' trial_effects <- data.frame(
#'   trial = c("TRIAL_A", "TRIAL_B", "TRIAL_C", "TRIAL_D"),
#'   estimate = c(0.3, 0.5, 0.2, 0.4),
#'   ci_lower = c(0.1, 0.2, -0.1, 0.1),
#'   ci_upper = c(0.5, 0.8, 0.5, 0.7),
#'   sample_size = c(100, 150, 80, 120)
#' )
#' 
#' # Basic forest plot
#' plot_forest(trial_effects, 
#'            parameter = "treatment_effect",
#'            xlab = "Standardized Mean Difference")
#' 
#' # Forest plot without summary effect
#' plot_forest(trial_effects,
#'            parameter = "treatment_effect", 
#'            xlab = "Effect Size (Cohen's d)",
#'            show_summary = FALSE)
#' 
#' # Order by sample size instead of effect size
#' plot_forest(trial_effects,
#'            parameter = "treatment_effect",
#'            order_by = "sample_size",
#'            xlab = "Treatment Effect")
#' 
#' # Subgroup analysis forest plot
#' subgroup_effects <- data.frame(
#'   subgroup = c("Age < 65", "Age >= 65", "Male", "Female"),
#'   effect_size = c(0.4, 0.2, 0.3, 0.5),
#'   lower = c(0.1, -0.1, 0.0, 0.2),
#'   upper = c(0.7, 0.5, 0.6, 0.8)
#' )
#' 
#' plot_forest(subgroup_effects,
#'            parameter = "subgroup_effect",
#'            xlab = "Treatment Effect by Subgroup") +
#'   ggplot2::labs(
#'     title = "Subgroup Analysis of Treatment Effects",
#'     subtitle = "Primary efficacy endpoint",
#'     caption = "Error bars represent 95% confidence intervals"
#'   )
#' 
#' # Meta-analysis of hazard ratios (log scale)
#' survival_data <- data.frame(
#'   study = paste("Study", 1:5),
#'   hr = c(0.8, 0.7, 0.9, 0.6, 0.8),
#'   hr_lower = c(0.6, 0.5, 0.7, 0.4, 0.6),
#'   hr_upper = c(1.0, 0.9, 1.1, 0.8, 1.0)
#' )
#' 
#' # Transform to log scale for plotting
#' survival_data$log_hr <- log(survival_data$hr)
#' survival_data$log_lower <- log(survival_data$hr_lower)
#' survival_data$log_upper <- log(survival_data$hr_upper)
#' 
#' plot_forest(survival_data,
#'            parameter = "hazard_ratio",
#'            xlab = "Log Hazard Ratio") +
#'   ggplot2::scale_x_continuous(
#'     labels = function(x) paste("HR =", round(exp(x), 2))
#'   )
#' }
#'
#' @seealso 
#' \code{\link{plot_posterior}} for posterior distribution plots,
#' \code{\link{plot_model_comparison}} for comparing different models,
#' \code{\link{extract_treatment_effect}} for extracting effects from models,
#' \code{\link{run_sensitivity_analyses}} for generating sensitivity data
#'
#' @export
plot_forest <- function(effect_data,
                      parameter = "treatment",
                      order_by = "effect_size",
                      xlab = "Effect Size",
                      show_summary = TRUE) {
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Check input data
  if (!is.data.frame(effect_data)) {
    stop("effect_data must be a data frame", call. = FALSE)
  }
  
  # Determine required columns
  required_cols <- c("trial", "estimate", "ci_lower", "ci_upper")
  
  # Check different common column names
  estimate_cols <- c("estimate", "mean", "effect", "effect_size")
  ci_lower_cols <- c("ci_lower", "lower", "lower_ci", "lb", "q2.5")
  ci_upper_cols <- c("ci_upper", "upper", "upper_ci", "ub", "q97.5")
  
  # Find best matching columns
  estimate_col <- intersect(estimate_cols, names(effect_data))[1]
  ci_lower_col <- intersect(ci_lower_cols, names(effect_data))[1]
  ci_upper_col <- intersect(ci_upper_cols, names(effect_data))[1]
  
  # Check if we have all required info
  missing_cols <- character(0)
  
  if (is.null(estimate_col)) missing_cols <- c(missing_cols, "estimate")
  if (is.null(ci_lower_col)) missing_cols <- c(missing_cols, "ci_lower")
  if (is.null(ci_upper_col)) missing_cols <- c(missing_cols, "ci_upper")
  
  if (length(missing_cols) > 0) {
    stop(paste("Required columns missing from effect_data:", paste(missing_cols, collapse = ", ")), 
         call. = FALSE)
  }
  
  # Get trial/group column
  if ("trial" %in% names(effect_data)) {
    trial_col <- "trial"
  } else if ("group" %in% names(effect_data)) {
    trial_col <- "group"
  } else if ("study" %in% names(effect_data)) {
    trial_col <- "study"
  } else if ("subgroup" %in% names(effect_data)) {
    trial_col <- "subgroup"
  } else {
    # Use first character column as trial identifier
    char_cols <- names(effect_data)[sapply(effect_data, is.character)]
    if (length(char_cols) > 0) {
      trial_col <- char_cols[1]
    } else {
      stop("Could not identify a column for trial/group names", call. = FALSE)
    }
  }
  
  # Prepare data for plotting
  plot_data <- effect_data
  names(plot_data)[names(plot_data) == estimate_col] <- "estimate"
  names(plot_data)[names(plot_data) == ci_lower_col] <- "ci_lower"
  names(plot_data)[names(plot_data) == ci_upper_col] <- "ci_upper"
  names(plot_data)[names(plot_data) == trial_col] <- "trial"
  
  # Order by requested variable
  if (order_by == "effect_size") {
    plot_data <- plot_data[order(plot_data$estimate), ]
  } else if (order_by %in% names(plot_data)) {
    plot_data <- plot_data[order(plot_data[[order_by]]), ]
  }
  
  # Convert trial to factor with ordered levels
  plot_data$trial <- factor(plot_data$trial, levels = unique(plot_data$trial))
  
  # Calculate overall summary if requested
  if (show_summary && nrow(plot_data) > 1) {
    # Simple inverse-variance weighted mean
    weights <- 1 / ((plot_data$ci_upper - plot_data$ci_lower) / 3.92)^2
    weighted_mean <- sum(plot_data$estimate * weights) / sum(weights)
    
    # Approximate CI for weighted mean
    weighted_var <- 1 / sum(weights)
    weighted_lower <- weighted_mean - 1.96 * sqrt(weighted_var)
    weighted_upper <- weighted_mean + 1.96 * sqrt(weighted_var)
    
    # Add to plot data
    summary_row <- data.frame(
      trial = "Overall",
      estimate = weighted_mean,
      ci_lower = weighted_lower,
      ci_upper = weighted_upper,
      stringsAsFactors = FALSE
    )
    
    # Add any other columns in plot_data with NA
    for (col in setdiff(names(plot_data), names(summary_row))) {
      summary_row[[col]] <- NA
    }
    
    # Combine with plot data
    plot_data <- rbind(plot_data, summary_row)
    
    # Update factor levels to put summary at bottom
    plot_data$trial <- factor(plot_data$trial, 
                            levels = c(setdiff(levels(plot_data$trial), "Overall"), "Overall"))
  }
  
  # Create forest plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = trial, x = estimate)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      height = 0.2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = ifelse(trial == "Overall", "Summary", "Study")),
      shape = 15
    ) +
    ggplot2::scale_size_manual(values = c("Summary" = 5, "Study" = 3)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste("Forest Plot for", parameter),
      x = xlab,
      y = NULL
    )
  
  # Add data table if possible
  if (requireNamespace("ggtext", quietly = TRUE)) {
    # Format estimates and CIs
    plot_data$formatted <- paste0(
      sprintf("%.2f", plot_data$estimate),
      " [",
      sprintf("%.2f", plot_data$ci_lower),
      ", ",
      sprintf("%.2f", plot_data$ci_upper),
      "]"
    )
    
    # Add text labels
    p <- p + ggplot2::geom_text(
      ggplot2::aes(x = max(ci_upper) * 1.1, label = formatted),
      hjust = 0
    ) +
      ggplot2::scale_x_continuous(
        limits = c(min(plot_data$ci_lower) * 1.1, max(plot_data$ci_upper) * 1.3)
      )
  }
  
  return(p)
}

#' Visualize treatment by covariate interactions in clinical trials
#'
#' Creates interaction plots to explore how treatment effects vary across 
#' patient subgroups or covariate values. Essential for personalized medicine
#' and clinical publications investigating effect modification. Supports both
#' continuous and categorical effect modifiers.
#'
#' @param model A BayesianModel object from \code{\link{fit_model}} or a 
#'   brmsfit object. Must contain interaction terms between treatment and 
#'   modifier variables. Model should include terms like:
#'   \itemize{
#'     \item \code{treatment * age} - Treatment by continuous covariate
#'     \item \code{treatment * sex} - Treatment by categorical covariate
#'     \item \code{treatment * baseline_severity} - Treatment by clinical measure
#'   }
#'   
#' @param treatment_var Character string naming the treatment variable in the model.
#'   Must match exactly the variable name used in model fitting:
#'   \itemize{
#'     \item \code{"treatment"} - Binary treatment indicator
#'     \item \code{"dose_group"} - Multi-level dose comparison
#'     \item \code{"intervention"} - Named intervention factor
#'   }
#'   
#' @param modifier_var Character string naming the effect modifier variable.
#'   Can be continuous or categorical:
#'   \itemize{
#'     \item \strong{Continuous:} age, weight, baseline_score, biomarker_level
#'     \item \strong{Categorical:} sex, race, disease_stage, comorbidity_group
#'   }
#'   
#' @param data Data frame containing the variables (optional). If not provided,
#'   function will attempt to extract from the model object. Should contain
#'   the same data used for model fitting.
#'   
#' @param type Character string specifying interaction plot type:
#'   \itemize{
#'     \item \code{"continuous"} - Line plots for continuous modifiers.
#'       Shows predicted outcomes across modifier range for each treatment group.
#'     \item \code{"categorical"} - Point plots for categorical modifiers.
#'       Shows predicted outcomes for each category, grouped by treatment.
#'   }
#'   Function auto-detects appropriate type based on modifier variable class.
#'
#' @return A ggplot object showing the interaction. Plot type depends on 
#'   modifier variable:
#'   \itemize{
#'     \item \strong{Continuous:} Lines for each treatment group with confidence ribbons
#'     \item \strong{Categorical:} Points with error bars for each treatment-category combination
#'   }
#'
#' @details
#' \strong{Clinical Interpretation Guidelines:}
#' 
#' \strong{Interaction Patterns:}
#' \itemize{
#'   \item \strong{Parallel lines/points:} No interaction - consistent treatment effect
#'   \item \strong{Converging lines:} Diminishing treatment effect in one subgroup
#'   \item \strong{Crossing lines:} Qualitative interaction - treatment harmful in some patients
#'   \item \strong{Diverging lines:} Enhanced treatment effect in certain subgroups
#' }
#' 
#' \strong{Statistical Considerations:}
#' \itemize{
#'   \item Wide confidence intervals suggest inadequate power for interaction detection
#'   \item Formal interaction tests should accompany visual inspection
#'   \item Multiple testing correction needed for multiple subgroups
#' }
#' 
#' \strong{Clinical Implications:}
#' \itemize{
#'   \item Significant interactions may require subgroup-specific dosing recommendations
#'   \item Pre-specification of clinically relevant subgroups strengthens interpretation
#'   \item Biological plausibility should support observed interactions
#' }
#' 
#' \strong{For Continuous Modifiers:}
#' \itemize{
#'   \item Look for threshold effects or smooth transitions
#'   \item Consider clinical cutpoints for continuous variables (e.g., age 65)
#'   \item Evaluate extrapolation beyond observed data ranges
#' }
#' 
#' \strong{For Categorical Modifiers:}
#' \itemize{
#'   \item Compare effect sizes and confidence intervals across categories
#'   \item Consider clinical meaningfulness of differences between subgroups
#'   \item Evaluate whether interaction justifies stratified analysis
#' }
#'
#' @examples
#' \dontrun{
#' # Load data and create model with interactions
#' data(synthetic_trial_data)
#' priors <- default_priors()
#' 
#' # Model with treatment by age interaction
#' model_spec <- create_model_spec(
#'   outcome ~ treatment * age + sex + baseline_score,
#'   family = "gaussian"
#' )
#' 
#' model <- fit_model(model_spec, synthetic_trial_data, priors)
#' 
#' # Continuous interaction plot (age)
#' plot_interaction(model, 
#'                 treatment_var = "treatment",
#'                 modifier_var = "age",
#'                 type = "continuous")
#' 
#' # Categorical interaction plot (automatically detected)
#' plot_interaction(model,
#'                 treatment_var = "treatment", 
#'                 modifier_var = "sex")
#' 
#' # Model with treatment by disease stage interaction
#' model_stage <- create_model_spec(
#'   survival_time ~ treatment * disease_stage + age + sex,
#'   family = "gamma"
#' )
#' 
#' fitted_model <- fit_model(model_stage, clinical_data, priors)
#' 
#' # Categorical interaction with custom styling
#' p <- plot_interaction(fitted_model,
#'                      treatment_var = "treatment",
#'                      modifier_var = "disease_stage",
#'                      type = "categorical")
#' 
#' p + ggplot2::labs(
#'     title = "Treatment Efficacy by Disease Stage",
#'     subtitle = "Significant interaction detected (p < 0.05)",
#'     caption = "Error bars: 95% credible intervals from Bayesian model"
#'   ) +
#'   ggplot2::theme_minimal() +
#'   ggplot2::scale_color_manual(
#'     values = c("Control" = "red", "Treatment" = "blue"),
#'     name = "Study Arm"
#'   )
#' 
#' # Biomarker interaction (continuous)
#' biomarker_model <- create_model_spec(
#'   response ~ treatment * biomarker_level + age + sex,
#'   family = "binomial"
#' )
#' 
#' fitted_bio <- fit_model(biomarker_model, biomarker_data, priors)
#' 
#' plot_interaction(fitted_bio,
#'                 treatment_var = "treatment",
#'                 modifier_var = "biomarker_level",
#'                 type = "continuous") +
#'   ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", 
#'                      color = "gray", alpha = 0.7) +
#'   ggplot2::labs(
#'     title = "Biomarker-Treatment Interaction",
#'     subtitle = "Response probability by biomarker level",
#'     x = "Biomarker Level (standardized)",
#'     y = "Predicted Response Probability"
#'   )
#' }
#'
#' @seealso 
#' \code{\link{create_model_spec}} for specifying interaction models,
#' \code{\link{fit_model}} for Bayesian model fitting,
#' \code{\link{identify_effect_modifiers}} for automated interaction detection,
#' \code{\link{assess_subgroup_credibility}} for subgroup analysis guidelines
#'
#' @export
plot_interaction <- function(model,
                           treatment_var,
                           modifier_var,
                           data = NULL,
                           type = c("continuous", "categorical")) {
  
  type <- match.arg(type)
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Check model
  if (inherits(model, "BayesianModel")) {
    if (is.null(model$model)) {
      stop("Model has not been fitted", call. = FALSE)
    }
    
    if (is.null(data)) {
      data <- model$model$data
    }
  } else if (inherits(model, "brmsfit")) {
    if (is.null(data)) {
      data <- model$data
    }
  } else {
    stop("model must be a BayesianModel or brmsfit object", call. = FALSE)
  }
  
  # Check variables
  if (!treatment_var %in% names(data)) {
    stop(paste("Treatment variable", treatment_var, "not found in data"), call. = FALSE)
  }
  
  if (!modifier_var %in% names(data)) {
    stop(paste("Modifier variable", modifier_var, "not found in data"), call. = FALSE)
  }
  
  # Create plot based on type and variable types
  if (type == "continuous" && is.numeric(data[[modifier_var]])) {
    # Continuous modifier
    
    # Create prediction grid
    grid_points <- 100
    mod_range <- range(data[[modifier_var]], na.rm = TRUE)
    mod_seq <- seq(mod_range[1], mod_range[2], length.out = grid_points)
    
    # Create prediction data
    newdata <- data.frame(mod_value = mod_seq)
    names(newdata)[1] <- modifier_var
    
    # Handle centered version if used in model
    if (paste0(modifier_var, "_c") %in% names(data)) {
      newdata[[paste0(modifier_var, "_c")]] <- scale(mod_seq, 
                                                   center = mean(data[[modifier_var]], na.rm = TRUE),
                                                   scale = FALSE)
    }
    
    # Get response variable
    if (inherits(model, "BayesianModel")) {
      response_var <- model$model_spec$outcome
    } else {
      # For brmsfit, extract from formula
      response_var <- model$formula$resp
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
      
      # Get level names for legend
      treated_name <- trt_levels[2]
      control_name <- trt_levels[1]
    } else if (is.numeric(data[[treatment_var]]) && length(unique(data[[treatment_var]])) <= 2) {
      trt_values <- sort(unique(data[[treatment_var]]))
      newdata_treated[[treatment_var]] <- trt_values[2]
      newdata_control[[treatment_var]] <- trt_values[1]
      
      # Get values for legend
      treated_name <- paste(treatment_var, "=", trt_values[2])
      control_name <- paste(treatment_var, "=", trt_values[1])
    } else {
      # For other cases, use first two unique values
      trt_values <- unique(data[[treatment_var]])
      newdata_treated[[treatment_var]] <- trt_values[min(2, length(trt_values))]
      newdata_control[[treatment_var]] <- trt_values[1]
      
      # Get values for legend
      treated_name <- paste(treatment_var, "=", trt_values[min(2, length(trt_values))])
      control_name <- paste(treatment_var, "=", trt_values[1])
    }
    
    # Generate predictions
    if (inherits(model, "BayesianModel")) {
      brms_model <- model$model
    } else {
      brms_model <- model
    }
    
    pred_treated <- brms::posterior_epred(brms_model, newdata = newdata_treated)
    pred_control <- brms::posterior_epred(brms_model, newdata = newdata_control)
    
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
      group = rep(c(treated_name, control_name), each = grid_points),
      mean = c(treated_mean, control_mean),
      lower = c(treated_lower, control_lower),
      upper = c(treated_upper, control_upper),
      stringsAsFactors = FALSE
    )
    
    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = modifier, y = mean, color = group, fill = group)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
      ggplot2::labs(
        title = paste("Interaction between", treatment_var, "and", modifier_var),
        x = modifier_var,
        y = response_var,
        color = NULL,
        fill = NULL
      ) +
      ggplot2::theme_minimal()
    
  } else if (type == "categorical" || !is.numeric(data[[modifier_var]])) {
    # Categorical interaction
    
    # If numeric, discretize into quantile groups
    if (is.numeric(data[[modifier_var]])) {
      # Create quantile groups
      n_groups <- 5
      cuts <- stats::quantile(data[[modifier_var]], probs = seq(0, 1, length.out = n_groups + 1), 
                            na.rm = TRUE)
      
      # Create labels
      labels <- character(n_groups)
      for (i in 1:n_groups) {
        labels[i] <- sprintf("[%.2f, %.2f]", cuts[i], cuts[i+1])
      }
      
      # Create grouping variable
      data$mod_group <- cut(data[[modifier_var]], breaks = cuts, labels = labels, include.lowest = TRUE)
      mod_group_var <- "mod_group"
    } else {
      # Use as is for categorical variables
      mod_group_var <- modifier_var
    }
    
    # Get unique groups
    if (is.factor(data[[mod_group_var]])) {
      groups <- levels(data[[mod_group_var]])
    } else {
      groups <- unique(data[[mod_group_var]])
    }
    
    # Get response variable
    if (inherits(model, "BayesianModel")) {
      response_var <- model$model_spec$outcome
    } else {
      # For brmsfit, extract from formula
      response_var <- model$formula$resp
    }
    
    # Generate fitted values
    if (inherits(model, "BayesianModel")) {
      brms_model <- model$model
    } else {
      brms_model <- model
    }
    
    # Create prediction data frames
    pred_data <- list()
    
    for (group in groups) {
      # Create data for this group
      group_data <- data.frame(
        group = group
      )
      names(group_data)[1] <- mod_group_var
      
      # Add modifier value if different from group variable
      if (mod_group_var != modifier_var) {
        # For numeric, use group midpoint
        if (is.numeric(data[[modifier_var]])) {
          group_idx <- which(labels == group)
          if (length(group_idx) > 0) {
            midpoint <- (cuts[group_idx] + cuts[group_idx + 1]) / 2
            group_data[[modifier_var]] <- midpoint
            
            # Add centered version if needed
            if (paste0(modifier_var, "_c") %in% names(data)) {
              group_data[[paste0(modifier_var, "_c")]] <- scale(midpoint, 
                                                             center = mean(data[[modifier_var]], na.rm = TRUE),
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
      
      # Get treatment values
      if (is.factor(data[[treatment_var]])) {
        trt_levels <- levels(data[[treatment_var]])
        trt_values <- trt_levels
      } else if (is.numeric(data[[treatment_var]]) && length(unique(data[[treatment_var]])) <= 2) {
        trt_values <- sort(unique(data[[treatment_var]]))
      } else {
        # For other cases, use first two unique values
        trt_values <- unique(data[[treatment_var]])
        trt_values <- trt_values[1:min(2, length(trt_values))]
      }
      
      # Create versions for each treatment value
      for (trt in trt_values) {
        trt_data <- group_data
        trt_data[[treatment_var]] <- trt
        
        # Generate predictions
        pred <- brms::posterior_epred(brms_model, newdata = trt_data)
        
        # Calculate summary statistics
        pred_mean <- mean(pred)
        pred_lower <- stats::quantile(pred, 0.025)
        pred_upper <- stats::quantile(pred, 0.975)
        
        # Store results
        pred_key <- paste(group, trt, sep = "_")
        pred_data[[pred_key]] <- list(
          group = group,
          treatment = trt,
          mean = pred_mean,
          lower = pred_lower,
          upper = pred_upper
        )
      }
    }
    
    # Combine results
    plot_data <- do.call(rbind, lapply(pred_data, function(x) {
      data.frame(
        group = x$group,
        treatment = x$treatment,
        mean = x$mean,
        lower = x$lower,
        upper = x$upper,
        stringsAsFactors = FALSE
      )
    }))
    
    # Create categorical treatment variable if numeric
    if (is.numeric(plot_data$treatment)) {
      plot_data$treatment_factor <- factor(plot_data$treatment)
    } else {
      plot_data$treatment_factor <- factor(plot_data$treatment)
    }
    
    # Create categorical group variable
    plot_data$group_factor <- factor(plot_data$group)
    
    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group_factor, y = mean, 
                                               color = treatment_factor, group = treatment_factor)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3), size = 3) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lower, ymax = upper),
        position = ggplot2::position_dodge(width = 0.3),
        width = 0.2
      ) +
      ggplot2::labs(
        title = paste("Interaction between", treatment_var, "and", modifier_var),
        x = modifier_var,
        y = response_var,
        color = treatment_var
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    stop("Invalid combination of plot type and variable type", call. = FALSE)
  }
  
  return(p)
}

#' Compare Bayesian models through visualization
#'
#' Creates comparative plots showing how different models estimate the same 
#' parameter. Essential for model selection, sensitivity analysis, and demonstrating
#' robustness of findings across different modeling assumptions. Supports multiple
#' visualization approaches for different audiences and purposes.
#'
#' @param models Named list of BayesianModel or brmsfit objects to compare.
#'   Names will be used as model labels in the plot. Examples:
#'   \itemize{
#'     \item \code{list("Neutral Prior" = model1, "Skeptical Prior" = model2)}
#'     \item \code{list("Linear" = linear_model, "Non-linear" = spline_model)}
#'     \item \code{list("Complete Case" = cc_model, "Multiple Imputation" = mi_model)}
#'   }
#'   
#' @param parameter Character string naming the parameter to compare across models.
#'   Must exist in all models. Common examples:
#'   \itemize{
#'     \item \code{"b_treatment"} - Main treatment effect
#'     \item \code{"b_Intercept"} - Model intercept/baseline
#'     \item \code{"sigma"} - Residual standard deviation
#'     \item \code{"b_treatment:age"} - Interaction effect
#'   }
#'   Use \code{extract_posterior()} to see available parameter names.
#'   
#' @param plot_type Character string specifying comparison visualization:
#'   \itemize{
#'     \item \code{"forest"} - Forest plot with point estimates and intervals.
#'       Best for formal documentation and clinical comparisons.
#'     \item \code{"ridgeline"} - Overlapping density plots for full distributions.
#'       Ideal for showing distributional differences. Requires ggridges package.
#'     \item \code{"intervals"} - Multi-level interval plots (50%, 95% CIs).
#'       Good compromise showing uncertainty without full distributions.
#'   }
#'
#' @return A ggplot object showing the model comparison. All plots include:
#'   \itemize{
#'     \item Reference line at zero for effect parameters
#'     \item Model names as labels
#'     \item Consistent color/styling for easy interpretation
#'   }
#'
#' @details
#' \strong{Model Comparison Guidelines:}
#' 
#' \strong{Prior Sensitivity:}
#' \itemize{
#'   \item Compare models with different prior specifications
#'   \item Similar posteriors indicate robustness to prior choice
#'   \item Large differences suggest data are sparse or priors are strongly informative
#'   \item Document prior sensitivity in formal reports
#' }
#' 
#' \strong{Model Structure:}
#' \itemize{
#'   \item Compare linear vs. non-linear specifications
#'   \item Evaluate different covariate adjustments
#'   \item Test alternative error distributions (normal, t, skewed)
#'   \item Assess impact of outlier treatment approaches
#' }
#' 
#' \strong{Missing Data Handling:}
#' \itemize{
#'   \item Compare complete case vs. imputation approaches
#'   \item Evaluate different imputation models
#'   \item Document sensitivity to missingness assumptions
#' }
#' 
#' \strong{Interpretation Guidelines:}
#' \itemize{
#'   \item \strong{Overlapping intervals:} Models agree on parameter estimate
#'   \item \strong{Non-overlapping intervals:} Substantial model disagreement
#'   \item \strong{Similar means, different widths:} Uncertainty differences
#'   \item \strong{Different means, similar widths:} Systematic bias differences
#' }
#' 
#' \strong{Analysis Considerations:}
#' \itemize{
#'   \item Pre-specify primary analysis and sensitivity analyses
#'   \item Show consistency of treatment effect across reasonable model variations
#'   \item Document any model-dependent conclusions
#'   \item Include model comparison in statistical analysis plan
#' }
#'
#' @examples
#' \dontrun{
#' # Fit models with different priors for comparison
#' data(synthetic_trial_data)
#' model_spec <- create_model_spec(outcome ~ treatment + age, family = "gaussian")
#' 
#' # Different prior specifications
#' neutral_priors <- default_priors("neutral")
#' skeptical_priors <- default_priors("skeptical")
#' optimistic_priors <- default_priors("optimistic")
#' 
#' # Fit models
#' model_neutral <- fit_model(model_spec, synthetic_trial_data, neutral_priors)
#' model_skeptical <- fit_model(model_spec, synthetic_trial_data, skeptical_priors)
#' model_optimistic <- fit_model(model_spec, synthetic_trial_data, optimistic_priors)
#' 
#' # Compare treatment effects across prior specifications
#' models_list <- list(
#'   "Neutral Prior" = model_neutral,
#'   "Skeptical Prior" = model_skeptical,
#'   "Optimistic Prior" = model_optimistic
#' )
#' 
#' # Forest plot comparison
#' plot_model_comparison(models_list, 
#'                      parameter = "b_treatment",
#'                      plot_type = "forest")
#' 
#' # Ridgeline plot for full distributions
#' plot_model_comparison(models_list,
#'                      parameter = "b_treatment",
#'                      plot_type = "ridgeline") +
#'   ggplot2::labs(
#'     title = "Prior Sensitivity Analysis",
#'     subtitle = "Treatment effect estimates across prior specifications"
#'   )
#' 
#' # Compare model structures
#' linear_spec <- create_model_spec(outcome ~ treatment + age, family = "gaussian")
#' nonlinear_spec <- create_model_spec(outcome ~ treatment + splines::bs(age, 3), 
#'                                    family = "gaussian")
#' 
#' linear_model <- fit_model(linear_spec, data, priors)
#' nonlinear_model <- fit_model(nonlinear_spec, data, priors)
#' 
#' structure_models <- list(
#'   "Linear Age" = linear_model,
#'   "Spline Age" = nonlinear_model
#' )
#' 
#' plot_model_comparison(structure_models,
#'                      parameter = "b_treatment",
#'                      plot_type = "intervals")
#' 
#' # Missing data sensitivity
#' complete_case_model <- fit_model(model_spec, 
#'                                 data[complete.cases(data), ], 
#'                                 priors)
#' imputed_model <- fit_model(model_spec, imputed_data, priors)
#' 
#' missing_models <- list(
#'   "Complete Cases" = complete_case_model,
#'   "Multiple Imputation" = imputed_model
#' )
#' 
#' # Create comprehensive comparison plot
#' p <- plot_model_comparison(missing_models,
#'                           parameter = "b_treatment", 
#'                           plot_type = "forest")
#' 
#' final_plot <- p + 
#'   ggplot2::labs(
#'     title = "Missing Data Sensitivity Analysis",
#'     subtitle = "Treatment effect: Complete cases vs. Multiple imputation",
#'     caption = "95% credible intervals from Bayesian models",
#'     x = "Treatment Effect (standardized)"
#'   ) +
#'   ggplot2::theme_minimal()
#' }
#'
#' @seealso 
#' \code{\link{fit_model}} for fitting Bayesian models,
#' \code{\link{default_priors}} for different prior specifications,
#' \code{\link{run_sensitivity_analyses}} for systematic sensitivity analysis,
#' \code{\link{extract_posterior}} for extracting posterior samples
#'
#' @export
plot_model_comparison <- function(models,
                                 parameter,
                                 plot_type = c("forest", "ridgeline", "intervals")) {
  
  plot_type <- match.arg(plot_type)
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Check input
  if (!is.list(models)) {
    stop("models must be a list", call. = FALSE)
  }
  
  # Extract posteriors for the parameter from each model
  posteriors <- list()
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    if (inherits(model, "BayesianModel")) {
      if (!is.null(model$model)) {
        post_samples <- extract_posterior(model, parameter)
        
        if (ncol(post_samples) > 0) {
          # Use first matching parameter
          param_col <- grep(parameter, colnames(post_samples), value = TRUE)[1]
          posteriors[[model_name]] <- post_samples[[param_col]]
        } else {
          warning(paste("Parameter", parameter, "not found in model", model_name), call. = FALSE)
        }
      } else {
        warning(paste("Model", model_name, "has not been fitted"), call. = FALSE)
      }
    } else if (inherits(model, "brmsfit")) {
      post_samples <- brms::posterior_samples(model)
      
      if (any(grepl(parameter, colnames(post_samples)))) {
        # Use first matching parameter
        param_col <- grep(parameter, colnames(post_samples), value = TRUE)[1]
        posteriors[[model_name]] <- post_samples[[param_col]]
      } else {
        warning(paste("Parameter", parameter, "not found in model", model_name), call. = FALSE)
      }
    } else {
      warning(paste("Model", model_name, "is not a BayesianModel or brmsfit object"), call. = FALSE)
    }
  }
  
  # Check if we have any posteriors
  if (length(posteriors) == 0) {
    stop("No posteriors could be extracted from the models", call. = FALSE)
  }
  
  # Create plot based on plot type
  if (plot_type == "forest") {
    # Create summary data frame
    summary_data <- data.frame(
      model = character(0),
      estimate = numeric(0),
      ci_lower = numeric(0),
      ci_upper = numeric(0),
      stringsAsFactors = FALSE
    )
    
    for (model_name in names(posteriors)) {
      post <- posteriors[[model_name]]
      
      # Calculate summary statistics
      estimate <- mean(post)
      ci_lower <- stats::quantile(post, 0.025)
      ci_upper <- stats::quantile(post, 0.975)
      
      # Add to summary data
      summary_data <- rbind(summary_data, data.frame(
        model = model_name,
        estimate = estimate,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        stringsAsFactors = FALSE
      ))
    }
    
    # Create forest plot
    p <- plot_forest(
      summary_data, 
      parameter = parameter, 
      order_by = "estimate", 
      xlab = "Parameter Estimate",
      show_summary = FALSE
    )
    
    p <- p + ggplot2::labs(
      title = paste("Model Comparison for", parameter),
      y = "Model"
    )
    
  } else if (plot_type == "ridgeline") {
    # Check for ggridges package
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("ggridges package is required for ridgeline plots but not available", call. = FALSE)
    }
    
    # Combine posteriors into a data frame
    combined_data <- data.frame()
    
    for (model_name in names(posteriors)) {
      post <- posteriors[[model_name]]
      
      # Add to combined data
      model_data <- data.frame(
        model = model_name,
        value = post,
        stringsAsFactors = FALSE
      )
      
      combined_data <- rbind(combined_data, model_data)
    }
    
    # Create ridgeline plot
    p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = value, y = model, fill = model)) +
      ggridges::geom_density_ridges(alpha = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(
        title = paste("Posterior Distributions for", parameter),
        subtitle = "Comparison across models",
        x = "Parameter Value",
        y = NULL
      )
    
  } else if (plot_type == "intervals") {
    # Create summary data frame
    summary_data <- data.frame(
      model = character(0),
      mean = numeric(0),
      median = numeric(0),
      sd = numeric(0),
      q2.5 = numeric(0),
      q25 = numeric(0),
      q75 = numeric(0),
      q97.5 = numeric(0),
      stringsAsFactors = FALSE
    )
    
    for (model_name in names(posteriors)) {
      post <- posteriors[[model_name]]
      
      # Calculate summary statistics
      summary_data <- rbind(summary_data, data.frame(
        model = model_name,
        mean = mean(post),
        median = stats::median(post),
        sd = stats::sd(post),
        q2.5 = stats::quantile(post, 0.025),
        q25 = stats::quantile(post, 0.25),
        q75 = stats::quantile(post, 0.75),
        q97.5 = stats::quantile(post, 0.975),
        stringsAsFactors = FALSE
      ))
    }
    
    # Create interval plot
    p <- ggplot2::ggplot(summary_data, ggplot2::aes(y = model)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_point(ggplot2::aes(x = mean), size = 3, color = "blue") +
      ggplot2::geom_linerange(ggplot2::aes(xmin = q25, xmax = q75), 
                            size = 2, color = "skyblue") +
      ggplot2::geom_linerange(ggplot2::aes(xmin = q2.5, xmax = q97.5), 
                            size = 0.8, color = "skyblue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Parameter Estimates for", parameter),
        subtitle = "Comparison across models",
        x = "Parameter Value",
        y = NULL
      )
  }
  
  return(p)
}

#' Visualize predicted outcomes by treatment groups
#'
#' Creates plots showing predicted outcomes across treatment groups with uncertainty
#' quantification. Essential for communicating treatment effects to clinical audiences
#' and clinical audiences. Automatically adapts to continuous vs. binary outcomes
#' and supports stratification by additional variables.
#'
#' @param model A BayesianModel object from \code{\link{fit_model}} or brmsfit object.
#'   Model must be fitted and contain treatment effects. Works with various 
#'   outcome types:
#'   \itemize{
#'     \item Continuous outcomes (Gaussian, t-distributed)
#'     \item Binary outcomes (Binomial, Bernoulli)
#'     \item Count outcomes (Poisson, Negative Binomial)
#'   }
#'   
#' @param treatment_var Character string naming the treatment variable in the model.
#'   Must match the variable name used in model fitting:
#'   \itemize{
#'     \item \code{"treatment"} - Binary indicator (0/1, Control/Treatment)
#'     \item \code{"dose_group"} - Multi-level factor (Low/Medium/High)
#'     \item \code{"intervention_arm"} - Named treatment arms
#'   }
#'   
#' @param data Data frame containing the analysis dataset (optional). If not 
#'   provided, function attempts to extract from model object. Should be the 
#'   same data used for model fitting.
#'   
#' @param outcome_var Character string naming the outcome variable (optional).
#'   If not provided, function attempts to extract from model specification.
#'   Used for axis labels and clinical interpretation:
#'   \itemize{
#'     \item \code{"primary_endpoint"} - Main efficacy measure
#'     \item \code{"response_rate"} - Binary response indicator
#'     \item \code{"adverse_events"} - Safety outcome count
#'   }
#'   
#' @param color_var Character string naming a stratification variable (optional).
#'   Creates separate colors/groups within each treatment:
#'   \itemize{
#'     \item \code{"sex"} - Male/Female stratification
#'     \item \code{"age_group"} - Age-based subgroups
#'     \item \code{"disease_stage"} - Clinical severity strata
#'   }
#'   Variable must be present in the dataset.
#'
#' @return A ggplot object showing outcomes by treatment group:
#'   \itemize{
#'     \item \strong{Continuous outcomes:} Box plots with individual points and model-based confidence intervals
#'     \item \strong{Binary outcomes:} Bar plots showing proportions with credible intervals
#'   }
#'   Red error bars show model-based uncertainty (95% credible intervals).
#'
#' @details
#' \strong{Plot Interpretation Guidelines:}
#' 
#' \strong{For Continuous Outcomes:}
#' \itemize{
#'   \item Box plots show data distribution within each treatment group
#'   \item Individual points (jittered) show all observations
#'   \item Red error bars show model-based 95% credible intervals for group means
#'   \item Non-overlapping error bars suggest significant treatment differences
#' }
#' 
#' \strong{For Binary Outcomes:}
#' \itemize{
#'   \item Bar heights show response rates/event proportions
#'   \item Error bars show uncertainty in population response rates
#'   \item Y-axis scaled 0-1 for proportions
#'   \item Compare bar heights and interval overlap for treatment effect assessment
#' }
#' 
#' \strong{Clinical Communication:}
#' \itemize{
#'   \item Use for investigator meetings and clinical presentations
#'   \item Include sample sizes in plot caption
#'   \item Document any important baseline differences between groups
#'   \item Consider clinical significance thresholds when interpreting differences
#' }
#' 
#' \strong{Statistical Notes:}
#' \itemize{
#'   \item Model-based intervals account for covariate adjustment
#'   \item Bayesian credible intervals (not frequentist confidence intervals)
#'   \item Intervals reflect posterior uncertainty given the data and priors
#'   \item Stratified plots help identify effect modification
#' }
#'
#' @examples
#' \dontrun{
#' # Load data and fit model
#' data(synthetic_trial_data)
#' priors <- default_priors()
#' model_spec <- create_model_spec(
#'   outcome ~ treatment + age + sex + baseline_score,
#'   family = "gaussian"
#' )
#' 
#' model <- fit_model(model_spec, synthetic_trial_data, priors)
#' 
#' # Basic outcomes plot
#' plot_outcomes(model, treatment_var = "treatment")
#' 
#' # Stratified by sex
#' plot_outcomes(model, 
#'              treatment_var = "treatment",
#'              color_var = "sex") +
#'   ggplot2::labs(
#'     title = "Treatment Effects by Sex",
#'     subtitle = "Primary efficacy endpoint",
#'     caption = "N=500. Error bars: 95% credible intervals."
#'   )
#' 
#' # Binary outcome example
#' response_spec <- create_model_spec(
#'   response ~ treatment + age + baseline_severity,
#'   family = "binomial"
#' )
#' 
#' response_model <- fit_model(response_spec, trial_data, priors)
#' 
#' plot_outcomes(response_model,
#'              treatment_var = "treatment",
#'              outcome_var = "response") +
#'   ggplot2::labs(
#'     title = "Response Rates by Treatment Group",
#'     subtitle = "Adjusted for age and baseline severity",
#'     y = "Response Rate",
#'     caption = "Bayesian logistic regression model"
#'   )
#' 
#' # Multi-arm trial visualization
#' dose_spec <- create_model_spec(
#'   efficacy_score ~ dose_group + age + weight,
#'   family = "gaussian"
#' )
#' 
#' dose_model <- fit_model(dose_spec, dose_data, priors)
#' 
#' plot_outcomes(dose_model,
#'              treatment_var = "dose_group",
#'              outcome_var = "efficacy_score",
#'              color_var = "age_group") +
#'   ggplot2::labs(
#'     title = "Dose-Response Relationship",
#'     subtitle = "Efficacy by dose group and age",
#'     x = "Dose Group",
#'     y = "Efficacy Score (0-100)"
#'   ) +
#'   ggplot2::theme_minimal()
#' 
#' # Safety outcome (count data)
#' safety_spec <- create_model_spec(
#'   adverse_events ~ treatment + age + comorbidity_count,
#'   family = "poisson"
#' )
#' 
#' safety_model <- fit_model(safety_spec, safety_data, priors)
#' 
#' plot_outcomes(safety_model,
#'              treatment_var = "treatment",
#'              outcome_var = "adverse_events") +
#'   ggplot2::labs(
#'     title = "Adverse Event Rates",
#'     subtitle = "Count per patient by treatment group",
#'     y = "Mean Adverse Events per Patient",
#'     caption = "Poisson regression adjusted for age and comorbidities"
#'   )
#' }
#'
#' @seealso 
#' \code{\link{fit_model}} for Bayesian model fitting,
#' \code{\link{plot_posterior}} for parameter-focused visualizations,
#' \code{\link{plot_interaction}} for exploring effect modification,
#' \code{\link{extract_treatment_effect}} for numerical effect estimates
#'
#' @export
plot_outcomes <- function(model,
                         treatment_var,
                         data = NULL,
                         outcome_var = NULL,
                         color_var = NULL) {
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Check model
  if (inherits(model, "BayesianModel")) {
    if (is.null(model$model)) {
      stop("Model has not been fitted", call. = FALSE)
    }
    
    if (is.null(data)) {
      data <- model$model$data
    }
    
    if (is.null(outcome_var)) {
      outcome_var <- model$model_spec$outcome
    }
  } else if (inherits(model, "brmsfit")) {
    if (is.null(data)) {
      data <- model$data
    }
    
    if (is.null(outcome_var)) {
      outcome_var <- model$formula$resp
    }
  } else {
    stop("model must be a BayesianModel or brmsfit object", call. = FALSE)
  }
  
  # Check variables
  if (!treatment_var %in% names(data)) {
    stop(paste("Treatment variable", treatment_var, "not found in data"), call. = FALSE)
  }
  
  if (!outcome_var %in% names(data)) {
    stop(paste("Outcome variable", outcome_var, "not found in data"), call. = FALSE)
  }
  
  # Check color variable if provided
  if (!is.null(color_var) && !color_var %in% names(data)) {
    stop(paste("Color variable", color_var, "not found in data"), call. = FALSE)
  }
  
  # Get treatment values
  if (is.factor(data[[treatment_var]])) {
    trt_levels <- levels(data[[treatment_var]])
  } else {
    trt_levels <- sort(unique(data[[treatment_var]]))
  }
  
  # Generate fitted values
  if (inherits(model, "BayesianModel")) {
    brms_model <- model$model
  } else {
    brms_model <- model
  }
  
  # Get family
  family <- brms_model$family$family
  
  # Generate posterior predictions
  fitted_values <- brms::posterior_epred(brms_model)
  
  # Calculate mean and credible intervals
  fitted_mean <- colMeans(fitted_values)
  fitted_lower <- apply(fitted_values, 2, function(x) stats::quantile(x, 0.025))
  fitted_upper <- apply(fitted_values, 2, function(x) stats::quantile(x, 0.975))
  
  # Create plotting data frame
  plot_data <- data.frame(
    outcome = data[[outcome_var]],
    treatment = data[[treatment_var]],
    fitted = fitted_mean,
    lower = fitted_lower,
    upper = fitted_upper,
    stringsAsFactors = FALSE
  )
  
  # Add color variable if provided
  if (!is.null(color_var)) {
    plot_data$color <- data[[color_var]]
  }
  
  # Create plot based on outcome type
  if (is.numeric(data[[outcome_var]])) {
    # Continuous outcome
    if (is.null(color_var)) {
      # Basic plot without color variable
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = treatment, y = outcome)) +
        ggplot2::geom_boxplot(alpha = 0.5, fill = "skyblue") +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("Outcomes by", treatment_var),
          subtitle = paste("Based on", family, "model"),
          x = treatment_var,
          y = outcome_var
        )
    } else {
      # Plot with color variable
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = treatment, y = outcome, color = color)) +
        ggplot2::geom_boxplot(alpha = 0.5) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("Outcomes by", treatment_var),
          subtitle = paste("Based on", family, "model"),
          x = treatment_var,
          y = outcome_var,
          color = color_var
        )
    }
    
    # Add means and CIs
    trt_means <- tapply(plot_data$fitted, plot_data$treatment, mean)
    trt_lower <- tapply(plot_data$lower, plot_data$treatment, mean)
    trt_upper <- tapply(plot_data$upper, plot_data$treatment, mean)
    
    trt_summary <- data.frame(
      treatment = as.numeric(names(trt_means)),
      mean = as.numeric(trt_means),
      lower = as.numeric(trt_lower),
      upper = as.numeric(trt_upper),
      stringsAsFactors = FALSE
    )
    
    # Add to plot
    p <- p + ggplot2::geom_pointrange(
      data = trt_summary,
      ggplot2::aes(x = treatment, y = mean, ymin = lower, ymax = upper),
      color = "red",
      size = 1,
      position = ggplot2::position_dodge(width = 0.5)
    )
    
  } else {
    # Binary outcome
    
    # Create proportion data
    if (is.null(color_var)) {
      # Aggregate without color variable
      prop_data <- aggregate(
        cbind(outcome, fitted, lower, upper) ~ treatment, 
        data = plot_data, 
        FUN = mean
      )
      
      # Create plot
      p <- ggplot2::ggplot(prop_data, ggplot2::aes(x = treatment, y = outcome)) +
        ggplot2::geom_bar(stat = "identity", fill = "skyblue", alpha = 0.7) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = lower, ymax = upper),
          width = 0.2,
          color = "red"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("Outcomes by", treatment_var),
          subtitle = paste("Based on", family, "model"),
          x = treatment_var,
          y = paste("Proportion of", outcome_var)
        ) +
        ggplot2::ylim(0, 1)
    } else {
      # Aggregate with color variable
      prop_data <- aggregate(
        cbind(outcome, fitted, lower, upper) ~ treatment + color, 
        data = plot_data, 
        FUN = mean
      )
      
      # Create plot
      p <- ggplot2::ggplot(prop_data, ggplot2::aes(x = treatment, y = outcome, fill = color)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = lower, ymax = upper),
          width = 0.2,
          color = "red",
          position = ggplot2::position_dodge(width = 0.9)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("Outcomes by", treatment_var),
          subtitle = paste("Based on", family, "model"),
          x = treatment_var,
          y = paste("Proportion of", outcome_var),
          fill = color_var
        ) +
        ggplot2::ylim(0, 1)
    }
  }
  
  return(p)
}

#' Assess model fit through diagnostic visualizations
#'
#' Creates diagnostic plots to evaluate how well a Bayesian model fits the observed
#' data. Essential for model validation, identifying outliers, and ensuring model
#' assumptions are met before drawing clinical conclusions. Multiple plot types
#' address different aspects of model adequacy.
#'
#' @param model A BayesianModel object from \code{\link{fit_model}} or brmsfit object.
#'   Model must be fitted with posterior samples available. Works with various
#'   model families and outcome types.
#'   
#' @param type Character string specifying the diagnostic plot type:
#'   \itemize{
#'     \item \code{"predicted_vs_observed"} - Scatter plot of model predictions vs. actual data.
#'       Points should cluster around diagonal line. Systematic deviations indicate model bias.
#'     \item \code{"residuals"} - Residual plots showing unexplained variation.
#'       Includes fitted vs. residuals plot and residual histogram. Random scatter indicates good fit.
#'     \item \code{"pp_check"} - Posterior predictive checks comparing simulated vs. observed data.
#'       Uses brms built-in functionality. Overlapping distributions suggest adequate model.
#'   }
#'   
#' @param n_draws Integer specifying number of posterior draws to use for 
#'   posterior predictive checks. More draws give smoother distributions but 
#'   slower computation. Recommended range: 50-200 for exploratory analysis,
#'   500+ for final model assessment.
#'
#' @return A ggplot object (or combination via patchwork) showing model diagnostics:
#'   \itemize{
#'     \item \strong{Predicted vs. Observed:} Scatter plot with diagonal reference line and loess smoother
#'     \item \strong{Residuals:} Scatter plot with horizontal reference line, plus histogram
#'     \item \strong{PP Check:} Density overlay plot from brms::pp_check()
#'   }
#'
#' @details
#' \strong{Diagnostic Interpretation Guidelines:}
#' 
#' \strong{Predicted vs. Observed Plot:}
#' \itemize{
#'   \item \strong{Good fit:} Points cluster tightly around diagonal (y = x) line
#'   \item \strong{Systematic bias:} Curved pattern or consistent deviation from diagonal
#'   \item \strong{Heteroscedasticity:} Fan-shaped pattern indicating non-constant variance
#'   \item \strong{Outliers:} Points far from diagonal may be influential observations
#' }
#' 
#' \strong{Residual Plots:}
#' \itemize{
#'   \item \strong{Good fit:} Random scatter around horizontal line (residuals = 0)
#'   \item \strong{Non-linearity:} Curved patterns suggest missing non-linear terms
#'   \item \strong{Heteroscedasticity:} Increasing/decreasing variance with fitted values
#'   \item \strong{Normal residuals:} Histogram should be approximately bell-shaped
#' }
#' 
#' \strong{Posterior Predictive Checks:}
#' \itemize{
#'   \item \strong{Good fit:} Observed data (dark line) falls within simulated data envelope
#'   \item \strong{Poor fit:} Observed data systematically different from simulations
#'   \item \strong{Overdispersion:} Observed data more variable than model predicts
#'   \item \strong{Underdispersion:} Observed data less variable than model predicts
#' }
#' 
#' \strong{Clinical Implications:}
#' \itemize{
#'   \item Poor model fit may invalidate treatment effect estimates
#'   \item Document model diagnostic results in formal reports
#'   \item Consider alternative model specifications if diagnostics indicate problems
#'   \item Outliers may represent important clinical subgroups or data errors
#' }
#' 
#' \strong{Model Improvement Strategies:}
#' \itemize{
#'   \item \strong{Non-linearity:} Add polynomial terms, splines, or transformations
#'   \item \strong{Heteroscedasticity:} Consider robust error distributions (t-distribution)
#'   \item \strong{Overdispersion:} Use negative binomial instead of Poisson for counts
#'   \item \strong{Outliers:} Investigate data quality or consider robust modeling approaches
#' }
#'
#' @examples
#' \dontrun{
#' # Fit a model for diagnostic assessment
#' data(synthetic_trial_data)
#' priors <- default_priors()
#' model_spec <- create_model_spec(
#'   outcome ~ treatment + age + sex + baseline_score,
#'   family = "gaussian"
#' )
#' 
#' model <- fit_model(model_spec, synthetic_trial_data, priors)
#' 
#' # Basic predicted vs observed plot
#' plot_model_fit(model, type = "predicted_vs_observed")
#' 
#' # Comprehensive residual analysis
#' resid_plot <- plot_model_fit(model, type = "residuals")
#' print(resid_plot)
#' 
#' # Posterior predictive checks
#' pp_plot <- plot_model_fit(model, type = "pp_check", n_draws = 100)
#' pp_plot + 
#'   ggplot2::labs(
#'     title = "Posterior Predictive Check",
#'     subtitle = "Model adequacy assessment",
#'     caption = "Dark line: observed data. Light lines: simulated data from model."
#'   )
#' 
#' # Create comprehensive diagnostic panel
#' pred_obs <- plot_model_fit(model, type = "predicted_vs_observed")
#' residuals <- plot_model_fit(model, type = "residuals")
#' pp_check <- plot_model_fit(model, type = "pp_check")
#' 
#' # Combine plots if patchwork available
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   diagnostic_panel <- (pred_obs | pp_check) / residuals
#'   diagnostic_panel + 
#'     patchwork::plot_annotation(
#'       title = "Model Diagnostic Assessment",
#'       subtitle = "Evaluating fit adequacy for clinical trial model"
#'     )
#' }
#' 
#' # Model comparison through diagnostics
#' # Compare linear vs. non-linear age effect
#' linear_spec <- create_model_spec(
#'   outcome ~ treatment + age + sex,
#'   family = "gaussian"
#' )
#' 
#' nonlinear_spec <- create_model_spec(
#'   outcome ~ treatment + splines::bs(age, df = 3) + sex,
#'   family = "gaussian"
#' )
#' 
#' linear_model <- fit_model(linear_spec, synthetic_trial_data, priors)
#' nonlinear_model <- fit_model(nonlinear_spec, synthetic_trial_data, priors)
#' 
#' # Compare residual patterns
#' linear_resid <- plot_model_fit(linear_model, type = "residuals") +
#'   ggplot2::labs(title = "Linear Age Model")
#' 
#' nonlinear_resid <- plot_model_fit(nonlinear_model, type = "residuals") +
#'   ggplot2::labs(title = "Non-linear Age Model")
#' 
#' # Side-by-side comparison
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   linear_resid | nonlinear_resid
#' }
#' }
#'
#' @seealso 
#' \code{\link{fit_model}} for Bayesian model fitting,
#' \code{\link{check_diagnostics}} for numerical diagnostic summaries,
#' \code{\link{model_fit_stats}} for quantitative fit statistics,
#' \code{\link{posterior_predictive_pvalue}} for formal adequacy testing
#'
#' @export
plot_model_fit <- function(model,
                          type = c("predicted_vs_observed", "residuals", "pp_check"),
                          n_draws = 50) {
  
  type <- match.arg(type)
  
  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available", call. = FALSE)
  }
  
  # Check model
  if (inherits(model, "BayesianModel")) {
    if (is.null(model$model)) {
      stop("Model has not been fitted", call. = FALSE)
    }
    
    brms_model <- model$model
  } else if (inherits(model, "brmsfit")) {
    brms_model <- model
  } else {
    stop("model must be a BayesianModel or brmsfit object", call. = FALSE)
  }
  
  # Get data and response variable
  data <- brms_model$data
  
  # Extract response variable
  response_var <- brms_model$formula$resp
  
  # Generate plot based on type
  if (type == "predicted_vs_observed") {
    # Generate fitted values
    fitted_values <- brms::posterior_epred(brms_model)
    
    # Calculate mean fitted values
    fitted_mean <- colMeans(fitted_values)
    
    # Create data frame for plotting
    plot_data <- data.frame(
      observed = data[[response_var]],
      predicted = fitted_mean,
      stringsAsFactors = FALSE
    )
    
    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = predicted, y = observed)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "blue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Predicted vs Observed Values",
        x = "Predicted",
        y = "Observed"
      )
    
  } else if (type == "residuals") {
    # Generate fitted values
    fitted_values <- brms::posterior_epred(brms_model)
    
    # Calculate mean fitted values
    fitted_mean <- colMeans(fitted_values)
    
    # Calculate residuals
    residuals <- data[[response_var]] - fitted_mean
    
    # Create data frame for plotting
    plot_data <- data.frame(
      fitted = fitted_mean,
      residuals = residuals,
      stringsAsFactors = FALSE
    )
    
    # Create residual plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fitted, y = residuals)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "blue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Residual Plot",
        x = "Fitted Values",
        y = "Residuals"
      )
    
    # Add residual histogram
    p_hist <- ggplot2::ggplot(plot_data, ggplot2::aes(x = residuals)) +
      ggplot2::geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Residual Distribution",
        x = "Residuals",
        y = "Count"
      )
    
    # Combine plots if patchwork is available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- p / p_hist
    }
    
  } else if (type == "pp_check") {
    # Use brms pp_check function
    p <- brms::pp_check(brms_model, ndraws = n_draws)
  }
  
  return(p)
}