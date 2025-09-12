test_that("Data generation functions work correctly", {
  # Test generate_synthetic_data
  data_continuous <- generate_synthetic_data(
    n_observations = 100,
    n_variables = 5,
    outcome_type = "continuous",
    seed = 123
  )
  
  expect_s3_class(data_continuous, "data.frame")
  expect_equal(nrow(data_continuous), 100)
  expect_true("outcome_continuous" %in% names(data_continuous))
  expect_true("treatment" %in% names(data_continuous))
  expect_true("age_centered" %in% names(data_continuous))
  
  # Test get_variable_descriptions
  var_desc <- get_variable_descriptions(data_continuous)
  expect_s3_class(var_desc, "data.frame")
  expect_true(nrow(var_desc) > 0)
  expect_true(all(c("variable", "description", "type", "role") %in% names(var_desc)))
  
  # Test get_true_effects
  true_effects <- get_true_effects(data_continuous)
  expect_type(true_effects, "list")
  expect_true("treatment_main" %in% names(true_effects))
})

test_that("Binary outcome data generation works", {
  data_binary <- generate_synthetic_data(
    n_observations = 50,
    outcome_type = "binary",
    seed = 456
  )
  
  expect_true("outcome_binary" %in% names(data_binary))
  expect_true(all(data_binary$outcome_binary %in% c(0, 1)))
})

test_that("Count outcome data generation works", {
  data_count <- generate_synthetic_data(
    n_observations = 50,
    outcome_type = "count",
    seed = 789
  )
  
  expect_true("outcome_count" %in% names(data_count))
  expect_true(all(data_count$outcome_count >= 0))
  expect_true(all(data_count$outcome_count == round(data_count$outcome_count)))
})

test_that("Model specification functions work correctly", {
  # Test create_model_spec
  model_spec <- create_model_spec(
    "outcome_continuous",
    "treatment + age_centered",
    family = "gaussian",
    model_name = "test_model"
  )
  
  expect_s3_class(model_spec, "ModelSpecification")
  expect_equal(model_spec$outcome, "outcome_continuous")
  expect_equal(model_spec$predictors, "treatment + age_centered")
  expect_equal(model_spec$family, "gaussian")
  expect_equal(model_spec$model_name, "test_model")
  
  # Test validate_model_spec
  expect_true(validate_model_spec(model_spec))
  
  # Test with invalid family
  expect_error(
    create_model_spec("y", "x", family = "invalid_family"),
    "Invalid family"
  )
})

test_that("Prior specification functions work correctly", {
  # Test PriorSpecification creation
  priors <- PriorSpecification()
  expect_s3_class(priors, "PriorSpecification")
  expect_equal(nrow(priors$priors), 0)
  
  # Test add_prior
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
  expect_equal(nrow(priors$priors), 1)
  expect_equal(priors$priors$parameter[1], "b_treatment")
  expect_equal(priors$priors$distribution[1], "normal")
  expect_equal(priors$priors$location[1], 0)
  expect_equal(priors$priors$scale[1], 2.5)
  
  # Test add_prior with bounds
  priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)
  expect_equal(nrow(priors$priors), 2)
  expect_equal(priors$priors$lb[2], 0)
  
  # Test remove_prior
  priors <- remove_prior(priors, "b_treatment")
  expect_equal(nrow(priors$priors), 1)
  expect_equal(priors$priors$parameter[1], "sigma")
})

test_that("Standard priors creation works", {
  variables <- c("treatment", "age_centered", "baseline_severity_centered")
  
  # Test continuous outcome priors
  std_priors <- create_standard_priors(
    variables = variables,
    outcome_type = "continuous",
    prior_types = c("neutral", "skeptical")
  )
  
  expect_s3_class(std_priors, "PriorSpecification")
  expect_true(nrow(std_priors$priors) > 0)
  
  # Should have priors for treatment and other variables
  param_names <- std_priors$priors$parameter
  expect_true(any(grepl("treatment", param_names)))
  expect_true(any(grepl("age_centered", param_names)))
})

test_that("Prior conversion to brms format works", {
  # Create simple priors
  priors <- PriorSpecification()
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
  priors <- add_prior(priors, "Intercept", "normal", 0, 10)
  priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)
  
  # Test conversion
  brms_priors <- to_brms_prior(priors)
  
  expect_s3_class(brms_priors, "brmsprior")
  expect_s3_class(brms_priors, "data.frame")
  expect_true(nrow(brms_priors) > 0)
  
  # Check specific priors are present
  expect_true(any(grepl("normal\\(0, 2.5\\)", brms_priors$prior)))
  expect_true(any(grepl("normal\\(0, 10\\)", brms_priors$prior)))
})

test_that("Parallel processing setup works", {
  # Test setup_parallel
  original_strategy <- setup_parallel("sequential", verbose = FALSE)
  expect_equal(original_strategy, "sequential")
  
  # Test get_parallel_config
  config <- get_parallel_config()
  expect_type(config, "list")
  expect_true("current_strategy" %in% names(config))
  expect_true("available_cores" %in% names(config))
  
  # Test reset_parallel
  reset_parallel(verbose = FALSE)
  
  # Verify reset worked
  config_after <- get_parallel_config()
  expect_equal(config_after$stored_strategy, "unknown")
})

test_that("Utility functions work correctly", {
  # Test set_seed
  set_seed(123)
  r1 <- rnorm(5)
  set_seed(123)
  r2 <- rnorm(5)
  expect_equal(r1, r2)
  
  # Test get_version
  version <- get_version()
  expect_type(version, "character")
  expect_true(grepl("\\d+\\.\\d+\\.\\d+", version))
  
  # Test get_session_info
  session_info <- get_session_info()
  expect_type(session_info, "list")
  expect_true("bayestrials_version" %in% names(session_info))
  expect_true("r_version" %in% names(session_info))
})