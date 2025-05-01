#' Plot posterior distributions
#' 
#' @param posterior_data Data frame or matrix of posterior samples
#' @param plot_type Type of posterior plot to create
#' @param facet_by Variable to facet plots by
#' @param color_by Variable to color plots by
#' @param theme Plot theme to use
#' @return A ggplot object
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

#' Plot forest plot for comparison
#' 
#' @param effect_data Data frame with effect estimates
#' @param parameter Name of the parameter to plot
#' @param order_by Variable to order points by
#' @param xlab Label for x-axis
#' @param show_summary Logical indicating whether to show summary effect
#' @return A ggplot object
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

#' Plot interaction effects
#' 
#' @param model A BayesianModel object
#' @param treatment_var Name of the treatment variable
#' @param modifier_var Name of the effect modifier variable
#' @param data Data frame (optional)
#' @param type Type of interaction plot
#' @return A ggplot object
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

#' Plot model comparison
#' 
#' @param models List of BayesianModel objects
#' @param parameter Name of the parameter to compare
#' @param plot_type Type of comparison plot
#' @return A ggplot object
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

#' Plot predicted outcomes across treatment groups
#' 
#' @param model A BayesianModel object
#' @param treatment_var Name of the treatment variable
#' @param data Data frame (optional)
#' @param outcome_var Name of the outcome variable (optional)
#' @param color_var Name of a variable to color by (optional)
#' @return A ggplot object
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

#' Plot model fit
#' 
#' @param model A BayesianModel object
#' @param type Type of fit plot
#' @param n_draws Number of posterior draws to plot
#' @return A ggplot object
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