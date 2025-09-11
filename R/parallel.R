#' Setup Parallel Processing for Bayesian Models
#'
#' Configures parallel processing backend for fitting multiple Bayesian models.
#'
#' @param strategy Parallel processing strategy ("multisession", "multicore", "sequential")
#' @param workers Number of workers (NULL for automatic detection)
#' @param verbose Whether to print setup information
#'
#' @return Invisibly returns the strategy used
#' @export
#'
#' @examples
#' \dontrun{
#' # Setup parallel processing
#' setup_parallel("multisession", workers = 4)
#' 
#' # Use sequential processing (no parallelization)
#' setup_parallel("sequential")
#' }
setup_parallel <- function(strategy = c("multisession", "multicore", "sequential"),
                          workers = NULL,
                          verbose = TRUE) {
  
  strategy <- match.arg(strategy)
  
  # Check for required packages
  if (!requireNamespace("future", quietly = TRUE)) {
    cli::cli_warn("Package 'future' not available, using sequential processing")
    strategy <- "sequential"
  }
  
  if (strategy != "sequential" && !requireNamespace("furrr", quietly = TRUE)) {
    cli::cli_warn("Package 'furrr' not available, using sequential processing") 
    strategy <- "sequential"
  }
  
  # Set up parallel backend
  if (strategy == "sequential") {
    if (requireNamespace("future", quietly = TRUE)) {
      future::plan(future::sequential)
    }
    if (verbose) {
      cli::cli_inform("Using sequential processing")
    }
    
  } else {
    # Determine number of workers
    if (is.null(workers)) {
      workers <- min(parallel::detectCores() - 1, 8)  # Leave one core free, max 8
      workers <- max(workers, 1)  # At least 1 worker
    }
    
    # Validate workers
    max_cores <- parallel::detectCores()
    if (workers > max_cores) {
      cli::cli_warn("Requested {workers} workers but only {max_cores} cores available")
      workers <- max_cores
    }
    
    # Set up future plan
    if (strategy == "multisession") {
      future::plan(future::multisession, workers = workers)
      if (verbose) {
        cli::cli_inform("Using multisession parallel processing with {workers} workers")
      }
      
    } else if (strategy == "multicore") {
      if (.Platform$OS.type == "windows") {
        cli::cli_warn("Multicore not supported on Windows, using multisession")
        future::plan(future::multisession, workers = workers)
        if (verbose) {
          cli::cli_inform("Using multisession parallel processing with {workers} workers")
        }
      } else {
        future::plan(future::multicore, workers = workers)
        if (verbose) {
          cli::cli_inform("Using multicore parallel processing with {workers} workers")
        }
      }
    }
    
    # Configure progressr if available
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::handlers(global = TRUE)
      if (verbose) {
        cli::cli_inform("Progress reporting enabled")
      }
    }
  }
  
  # Store configuration in options
  options(
    bayestrials.parallel.strategy = strategy,
    bayestrials.parallel.workers = workers
  )
  
  return(invisible(strategy))
}

#' Get Current Parallel Configuration
#'
#' Returns information about the current parallel processing configuration.
#'
#' @return A list with parallel configuration details
#' @export
get_parallel_config <- function() {
  
  # Get future plan info
  if (requireNamespace("future", quietly = TRUE)) {
    plan_info <- future::plan()
    strategy <- class(plan_info)[1]
    
    # Get number of workers
    if (strategy == "sequential") {
      workers <- 1
    } else {
      workers <- future::nbrOfWorkers()
    }
  } else {
    strategy <- "sequential"
    workers <- 1
  }
  
  # Get stored options
  stored_strategy <- getOption("bayestrials.parallel.strategy", "unknown")
  stored_workers <- getOption("bayestrials.parallel.workers", NA)
  
  config <- list(
    current_strategy = strategy,
    current_workers = workers,
    stored_strategy = stored_strategy,
    stored_workers = stored_workers,
    available_cores = parallel::detectCores(),
    future_workers = workers,
    progressr_available = requireNamespace("progressr", quietly = TRUE)
  )
  
  class(config) <- c("bayestrials_parallel_config", "list")
  return(config)
}

#' Print Parallel Configuration
#'
#' @param x A bayestrials_parallel_config object
#' @param ... Additional arguments (ignored)
#'
#' @return x invisibly
#' @export
print.bayestrials_parallel_config <- function(x, ...) {
  cli::cli_h2("Parallel Processing Configuration")
  
  cli::cli_ul(c(
    "Current strategy: {x$current_strategy}",
    "Current workers: {x$current_workers}",
    "Available cores: {x$available_cores}",
    "Progress reporting: {if (x$progressr_available) 'Available' else 'Not available'}"
  ))
  
  if (x$current_strategy != "sequential") {
    if (x$current_workers == x$available_cores) {
      cli::cli_alert_warning("Using all available cores - consider leaving one free")
    } else if (x$current_workers > x$available_cores * 0.8) {
      cli::cli_alert_info("Using most available cores - good for compute-intensive tasks")
    }
  }
  
  invisible(x)
}

#' Reset Parallel Processing
#'
#' Resets parallel processing to sequential mode and clears configuration.
#'
#' @param verbose Whether to print reset information
#'
#' @return Invisibly returns NULL
#' @export
reset_parallel <- function(verbose = TRUE) {
  
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  
  # Clear options
  options(
    bayestrials.parallel.strategy = NULL,
    bayestrials.parallel.workers = NULL
  )
  
  if (verbose) {
    cli::cli_inform("Parallel processing reset to sequential mode")
  }
  
  return(invisible(NULL))
}

#' Fit Models in Parallel
#'
#' Internal function to fit multiple models in parallel using the configured backend.
#'
#' @param model_specs List of model specifications
#' @param data Data frame for fitting
#' @param priors_list List of prior specifications
#' @param fit_args Additional arguments for model fitting
#' @param progress Whether to show progress
#'
#' @return List of fitted models
#' @keywords internal
fit_models_parallel <- function(model_specs, data, priors_list = NULL, fit_args = list(), progress = TRUE) {
  
  # Check if parallel processing is available
  config <- get_parallel_config()
  use_parallel <- config$current_strategy != "sequential" && 
                 requireNamespace("furrr", quietly = TRUE)
  
  # Create combinations of models and priors
  if (is.null(priors_list)) {
    priors_list <- list(default = NULL)
  }
  
  combinations <- tidyr::crossing(
    model_idx = seq_along(model_specs),
    prior_idx = seq_along(priors_list)
  )
  
  # Function to fit a single model
  fit_single_model <- function(combo_idx) {
    combo <- combinations[combo_idx, ]
    model_spec <- model_specs[[combo$model_idx]]
    priors <- priors_list[[combo$prior_idx]]
    
    # Set up model name
    model_name <- model_spec$model_name %||% paste0("model_", combo$model_idx)
    prior_name <- names(priors_list)[combo$prior_idx] %||% "default"
    
    # Attach priors to model spec
    if (!is.null(priors)) {
      attr(model_spec, "priors") <- priors
    }
    
    # Fit model
    result <- tryCatch({
      do.call(fit_model, c(list(model_spec = model_spec, data = data), fit_args))
    }, error = function(e) {
      cli::cli_warn("Error fitting {model_name} with {prior_name} priors: {e$message}")
      NULL
    })
    
    # Add metadata
    if (!is.null(result)) {
      result$model_name <- model_name
      result$prior_name <- prior_name
    }
    
    return(result)
  }
  
  # Fit models
  if (use_parallel) {
    if (progress && requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(along = seq_len(nrow(combinations)))
        
        results <- furrr::future_map(seq_len(nrow(combinations)), function(i) {
          p(sprintf("Fitting model %d/%d", i, nrow(combinations)))
          fit_single_model(i)
        }, .options = furrr::furrr_options(seed = TRUE))
      })
    } else {
      results <- furrr::future_map(seq_len(nrow(combinations)), fit_single_model,
                                   .options = furrr::furrr_options(seed = TRUE))
    }
  } else {
    # Sequential processing
    results <- list()
    if (progress) {
      cli::cli_progress_bar("Fitting models", total = nrow(combinations))
    }
    
    for (i in seq_len(nrow(combinations))) {
      if (progress) {
        cli::cli_progress_update()
      }
      results[[i]] <- fit_single_model(i)
    }
    
    if (progress) {
      cli::cli_progress_done()
    }
  }
  
  # Remove NULL results and name appropriately
  results <- results[!sapply(results, is.null)]
  
  # Create names for results
  names(results) <- sapply(results, function(r) {
    if (!is.null(r)) {
      paste(r$model_name, r$prior_name, sep = "_")
    } else {
      NA
    }
  })
  
  names(results) <- names(results)[!is.na(names(results))]
  
  return(results)
}