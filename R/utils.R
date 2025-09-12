#' Check if a package is installed and load it
#' 
#' @param pkg Character string of package name
#' @param min_version Minimum version required (optional)
#' @return Logical indicating if package was successfully loaded
#' @keywords internal
check_and_load_package <- function(pkg, min_version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg, "' is required but not installed. ",
               "Please install it with install.packages('", pkg, "')"), 
         call. = FALSE)
  }
  
  if (!is.null(min_version)) {
    pkg_version <- utils::packageVersion(pkg)
    if (pkg_version < min_version) {
      stop(paste0("Package '", pkg, "' version ", pkg_version, " is installed, but version >= ", 
                 min_version, " is required. Please update it with install.packages('", pkg, "')"),
           call. = FALSE)
    }
  }
  
  return(TRUE)
}

#' Check if a file exists and is readable
#' 
#' @param file_path Character string of file path
#' @param extension Expected file extension (optional)
#' @return Logical indicating if file exists and is readable
#' @keywords internal
check_file <- function(file_path, extension = NULL) {
  if (!file.exists(file_path)) {
    stop(paste0("File does not exist: ", file_path), call. = FALSE)
  }
  
  if (!is.null(extension)) {
    actual_ext <- tools::file_ext(file_path)
    if (actual_ext != extension) {
      stop(paste0("File '", file_path, "' does not have the expected extension '", 
                 extension, "'. Found: '", actual_ext, "'"), call. = FALSE)
    }
  }
  
  if (!file.access(file_path, mode = 4) == 0) {
    stop(paste0("File exists but is not readable: ", file_path), call. = FALSE)
  }
  
  return(TRUE)
}

#' Set a global random seed in a reproducible way
#' 
#' @param seed Integer seed for random number generation
#' @return Invisibly returns the seed used
#' @export
set_seed <- function(seed = NULL) {
  if (is.null(seed)) {
    seed <- as.integer(Sys.time())
  }
  
  # Ensure seed is a valid integer
  seed <- as.integer(seed)
  
  set.seed(seed)
  
  # For brms/Stan, seed is passed to the sampling function directly
  
  invisible(seed)
}

#' Create a standardized log message
#' 
#' @param msg The message to log
#' @param level Character string of log level (info, warning, error, debug)
#' @param timestamp Logical indicating whether to include timestamp
#' @return Character string with formatted log message
#' @keywords internal
create_log_message <- function(msg, level = "info", timestamp = TRUE) {
  level <- match.arg(level, c("info", "warning", "error", "debug"))
  
  level_prefix <- switch(level,
                         info = "[INFO] ",
                         warning = "[WARNING] ",
                         error = "[ERROR] ",
                         debug = "[DEBUG] ")
  
  time_prefix <- if (timestamp) paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ") else ""
  
  paste0(time_prefix, level_prefix, msg)
}

#' Log a message to console or file
#' 
#' @param msg The message to log
#' @param level Character string of log level (info, warning, error, debug)
#' @param file Optional file path to write log to
#' @param console Logical indicating whether to print to console
#' @return Invisibly returns the formatted log message
#' @export
log_message <- function(msg, level = "info", file = NULL, console = TRUE) {
  level <- match.arg(level, c("info", "warning", "error", "debug"))
  formatted_msg <- create_log_message(msg, level)
  
  if (console) {
    switch(level,
           info = cat(formatted_msg, "\n"),
           warning = warning(formatted_msg, call. = FALSE),
           error = stop(formatted_msg, call. = FALSE),
           debug = if (getOption("bayestrials.debug", FALSE)) cat(formatted_msg, "\n"))
  }
  
  if (!is.null(file)) {
    cat(formatted_msg, "\n", file = file, append = TRUE)
  }
  
  invisible(formatted_msg)
}

#' Get package option with fallback to default
#' 
#' @param option_name Name of the option to retrieve
#' @param default Default value if option is not set
#' @return The option value or default
#' @keywords internal
get_option <- function(option_name, default = NULL) {
  full_name <- paste0("bayestrials.", option_name)
  getOption(full_name, default)
}

#' Set package option
#' 
#' @param option_name Name of the option to set
#' @param value Value to set for the option
#' @return Invisibly returns the previous option value
#' @export
set_option <- function(option_name, value) {
  full_name <- paste0("bayestrials.", option_name)
  old_value <- getOption(full_name)
  options(setNames(list(value), full_name))
  invisible(old_value)
}

#' Format a value for printing based on its type
#' 
#' @param x The value to format
#' @return Character string representation of the value
#' @keywords internal
format_value <- function(x) {
  if (is.null(x)) {
    return("NULL")
  } else if (is.character(x)) {
    return(paste0("\"", x, "\""))
  } else if (is.logical(x)) {
    return(ifelse(x, "TRUE", "FALSE"))
  } else if (is.factor(x)) {
    return(paste0("factor(\"", as.character(x), "\")"))
  } else if (is.formula(x)) {
    return(paste0("formula(", deparse(x), ")"))
  } else {
    return(as.character(x))
  }
}

#' Check if a variable is a formula
#' 
#' @param x Object to check
#' @return Logical indicating if x is a formula
#' @keywords internal
is.formula <- function(x) {
  inherits(x, "formula")
}

#' Convert a string to a formula
#' 
#' @param formula_str String representing a formula
#' @return A formula object
#' @keywords internal
str_to_formula <- function(formula_str) {
  if (is.formula(formula_str)) {
    return(formula_str)
  }
  
  as.formula(formula_str)
}

#' Check if object is of a specific class
#' 
#' @param x Object to check
#' @param class_name Name of the class to check for
#' @return Logical indicating if object inherits from class_name
#' @keywords internal
is_class <- function(x, class_name) {
  inherits(x, class_name)
}

#' Get current working directory and format as absolute path
#' 
#' @return Character string of absolute path to current working directory
#' @keywords internal
get_working_dir <- function() {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

#' Create a timestamped directory name
#' 
#' @param prefix Prefix for the directory name
#' @return Character string with prefix and timestamp
#' @keywords internal
create_timestamped_dir <- function(prefix = "analysis") {
  paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
}

#' Get package version
#' 
#' @return Character string with package version
#' @export
get_version <- function() {
  as.character(utils::packageVersion("bayestrials"))
}

#' Get session info including package versions
#' 
#' @return List containing session information
#' @export
get_session_info <- function() {
  session_info <- list(
    bayestrials_version = get_version(),
    r_version = getRversion(),
    platform = Sys.info()["sysname"],
    date = Sys.Date(),
    locale = Sys.getlocale()
  )
  
  if (requireNamespace("brms", quietly = TRUE)) {
    session_info$brms_version <- as.character(utils::packageVersion("brms"))
  }
  
  if (requireNamespace("rstan", quietly = TRUE)) {
    session_info$rstan_version <- as.character(utils::packageVersion("rstan"))
  }
  
  return(session_info)
}