# Comprehensive Synthetic Data Workflow Tests
# Tests the complete bayestrials workflow using synthetic data

context("Synthetic Data Workflow")

# Test configuration for faster testing
TEST_CONFIG <- list(
  n_obs_quick = 50,
  n_obs_full = 200,
  n_iter = 500,
  n_chains = 1,
  seed = 42
)

#-----------------------------------------------------------------------------
# Data Generation Tests
#-----------------------------------------------------------------------------

test_that("Synthetic data generation works for all outcome types", {
  # Continuous outcome
  data_cont <- generate_synthetic_data(
    n_observations = TEST_CONFIG$n_obs_quick,
    outcome_type = "continuous",
    seed = TEST_CONFIG$seed
  )

  expect_s3_class(data_cont, "data.frame")
  expect_equal(nrow(data_cont), TEST_CONFIG$n_obs_quick)
  expect_true("outcome_continuous" %in% names(data_cont))
  expect_true("treatment" %in% names(data_cont))
  expect_true(all(data_cont$treatment %in% c(0, 1)))

  # Binary outcome
  data_bin <- generate_synthetic_data(
    n_observations = TEST_CONFIG$n_obs_quick,
    outcome_type = "binary",
    seed = TEST_CONFIG$seed + 1
  )

  expect_true("outcome_binary" %in% names(data_bin))
  expect_true(all(data_bin$outcome_binary %in% c(0, 1)))

  # Count outcome
  data_count <- generate_synthetic_data(
    n_observations = TEST_CONFIG$n_obs_quick,
    outcome_type = "count",
    seed = TEST_CONFIG$seed + 2
  )

  expect_true("outcome_count" %in% names(data_count))
  expect_true(all(data_count$outcome_count >= 0))
  expect_true(all(data_count$outcome_count == floor(data_count$outcome_count)))
})

test_that("True effects are stored and retrievable", {
  data <- generate_synthetic_data(
    n_observations = 100,
    treatment_effect = 0.7,
    seed = TEST_CONFIG$seed
  )

  true_effects <- get_true_effects(data)

  expect_type(true_effects, "list")
  expect_true("treatment_main" %in% names(true_effects))
  expect_equal(true_effects$treatment_main, 0.7)
  expect_true("intercept" %in% names(true_effects))
  expect_true("age_effect" %in% names(true_effects))
  expect_true("baseline_severity_effect" %in% names(true_effects))
})

test_that("Variable descriptions are accurate", {
  data <- generate_synthetic_data(
    n_observations = 50,
    outcome_type = "continuous",
    seed = TEST_CONFIG$seed
  )

  var_desc <- get_variable_descriptions(data)

  expect_s3_class(var_desc, "data.frame")
  expect_true(all(c("variable", "description", "type", "role") %in% names(var_desc)))

  # Check for expected variables
  expect_true("treatment" %in% var_desc$variable)
  expect_true("age_centered" %in% var_desc$variable)
  expect_true("outcome_continuous" %in% var_desc$variable)

  # Check roles
  expect_true("outcome" %in% var_desc$role)
  expect_true("treatment" %in% var_desc$role)
  expect_true("covariate" %in% var_desc$role)
})

#-----------------------------------------------------------------------------
# Model Specification Tests
#-----------------------------------------------------------------------------

test_that("Model specifications are created correctly for all families", {
  # Gaussian family
  spec_gauss <- create_model_spec(
    "outcome_continuous",
    "treatment + age_centered",
    family = "gaussian",
    model_name = "gauss_model"
  )

  expect_s3_class(spec_gauss, "ModelSpecification")
  expect_equal(spec_gauss$family, "gaussian")
  expect_equal(spec_gauss$outcome, "outcome_continuous")
  expect_true(validate_model_spec(spec_gauss))

  # Binomial family
  spec_bin <- create_model_spec(
    "outcome_binary",
    "treatment + age_centered",
    family = "binomial",
    link = "logit",
    model_name = "binary_model"
  )

  expect_equal(spec_bin$family, "binomial")
  expect_equal(spec_bin$link, "logit")

  # Poisson family
  spec_pois <- create_model_spec(
    "outcome_count",
    "treatment + age_centered",
    family = "poisson",
    model_name = "count_model"
  )

  expect_equal(spec_pois$family, "poisson")
})

#-----------------------------------------------------------------------------
# Prior Specification Tests
#-----------------------------------------------------------------------------

test_that("Prior specifications work correctly", {
  priors <- PriorSpecification()
  expect_s3_class(priors, "PriorSpecification")
  expect_equal(nrow(priors$priors), 0)

  # Add various priors
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
  priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
  priors <- add_prior(priors, "Intercept", "normal", 0, 10)
  priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)

  expect_equal(nrow(priors$priors), 4)

  # Test conversion to brms format
  brms_priors <- to_brms_prior(priors)
  expect_s3_class(brms_priors, "brmsprior")
  expect_equal(nrow(brms_priors), 4)
})

test_that("Standard priors can be created", {
  variables <- c("treatment", "age_centered", "baseline_severity_centered")

  std_priors <- create_standard_priors(
    variables = variables,
    outcome_type = "continuous",
    prior_types = c("neutral", "skeptical", "optimistic")
  )

  expect_s3_class(std_priors, "PriorSpecification")
  expect_true(nrow(std_priors$priors) > 0)
})

#-----------------------------------------------------------------------------
# Full Workflow Integration Tests (Slow - Skip on CRAN)
#-----------------------------------------------------------------------------

test_that("Complete workflow works with continuous outcome", {
  skip_if_not_installed("brms")
  skip_on_cran()

  # Generate data
  data <- generate_synthetic_data(
    n_observations = TEST_CONFIG$n_obs_full,
    outcome_type = "continuous",
    treatment_effect = 0.5,
    seed = TEST_CONFIG$seed
  )

  true_effects <- get_true_effects(data)

  # Create model spec
  model_spec <- create_model_spec(
    "outcome_continuous",
    "treatment + age_centered + baseline_severity_centered",
    family = "gaussian",
    model_name = "workflow_test"
  )

  # Create priors
  priors <- PriorSpecification()
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
  priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
  priors <- add_prior(priors, "b_baseline_severity_centered", "normal", 0, 2.0)
  priors <- add_prior(priors, "Intercept", "normal", 10, 10)
  priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)

  attr(model_spec, "priors") <- priors

  # Fit model
  result <- fit_model(
    model_spec = model_spec,
    data = data,
    chains = TEST_CONFIG$n_chains,
    iter = TEST_CONFIG$n_iter,
    cores = 1,
    seed = TEST_CONFIG$seed
  )

  expect_s3_class(result, "BayesianModel")
  expect_true(!is.null(result$model))

  # Check diagnostics
  diagnostics <- check_diagnostics(result)
  expect_type(diagnostics, "list")

  # Extract estimates
  estimates <- extract_estimates(result)
  expect_s3_class(estimates, "data.frame")
  expect_true("b_treatment" %in% estimates$parameter)

  # Check parameter recovery (treatment effect should be close to true value)
  treatment_est <- estimates$mean[estimates$parameter == "b_treatment"]
  treatment_true <- true_effects$treatment_main
  treatment_bias <- abs(treatment_est - treatment_true)

  # Treatment effect should be within reasonable range
  expect_lt(treatment_bias, 0.5)  # Allow some estimation error
})

test_that("Complete workflow works with binary outcome", {
  skip_if_not_installed("brms")
  skip_on_cran()

  # Generate data
  data <- generate_synthetic_data(
    n_observations = TEST_CONFIG$n_obs_full,
    outcome_type = "binary",
    treatment_effect = 0.5,
    seed = TEST_CONFIG$seed
  )

  # Create model spec
  model_spec <- create_model_spec(
    "outcome_binary",
    "treatment + age_centered",
    family = "binomial",
    link = "logit",
    model_name = "binary_workflow"
  )

  # Create priors
  priors <- PriorSpecification()
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.0)
  priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
  priors <- add_prior(priors, "Intercept", "normal", 0, 5)

  attr(model_spec, "priors") <- priors

  # Fit model
  result <- fit_model(
    model_spec = model_spec,
    data = data,
    chains = TEST_CONFIG$n_chains,
    iter = TEST_CONFIG$n_iter,
    cores = 1,
    seed = TEST_CONFIG$seed
  )

  expect_s3_class(result, "BayesianModel")

  # Check estimates
  estimates <- extract_estimates(result)
  expect_true("b_treatment" %in% estimates$parameter)
})

#-----------------------------------------------------------------------------
# Sensitivity Analysis Tests
#-----------------------------------------------------------------------------

test_that("Sensitivity analysis with different priors works", {
  skip_if_not_installed("brms")
  skip_on_cran()

  # Generate data
  data <- generate_synthetic_data(
    n_observations = TEST_CONFIG$n_obs_quick,
    outcome_type = "continuous",
    seed = TEST_CONFIG$seed
  )

  # Create neutral prior spec
  priors_neutral <- PriorSpecification()
  priors_neutral <- add_prior(priors_neutral, "b_treatment", "normal", 0, 2.5)
  priors_neutral <- add_prior(priors_neutral, "Intercept", "normal", 10, 10)
  priors_neutral <- add_prior(priors_neutral, "sigma", "half_cauchy", NA, 5, 0)

  # Create skeptical prior spec
  priors_skeptical <- PriorSpecification()
  priors_skeptical <- add_prior(priors_skeptical, "b_treatment", "normal", 0, 0.5)
  priors_skeptical <- add_prior(priors_skeptical, "Intercept", "normal", 10, 10)
  priors_skeptical <- add_prior(priors_skeptical, "sigma", "half_cauchy", NA, 5, 0)

  # Fit with neutral priors
  model_spec_neutral <- create_model_spec(
    "outcome_continuous", "treatment", family = "gaussian"
  )
  attr(model_spec_neutral, "priors") <- priors_neutral

  result_neutral <- fit_model(
    model_spec = model_spec_neutral,
    data = data,
    chains = 1, iter = TEST_CONFIG$n_iter, cores = 1, seed = TEST_CONFIG$seed
  )

  # Fit with skeptical priors
  model_spec_skeptical <- create_model_spec(
    "outcome_continuous", "treatment", family = "gaussian"
  )
  attr(model_spec_skeptical, "priors") <- priors_skeptical

  result_skeptical <- fit_model(
    model_spec = model_spec_skeptical,
    data = data,
    chains = 1, iter = TEST_CONFIG$n_iter, cores = 1, seed = TEST_CONFIG$seed
  )

  # Extract and compare estimates
  est_neutral <- extract_estimates(result_neutral)
  est_skeptical <- extract_estimates(result_skeptical)

  treat_neutral <- est_neutral$mean[est_neutral$parameter == "b_treatment"]
  treat_skeptical <- est_skeptical$mean[est_skeptical$parameter == "b_treatment"]

  # Both should produce valid estimates
  expect_false(is.na(treat_neutral))
  expect_false(is.na(treat_skeptical))

  # Skeptical prior should shrink estimate toward zero
  # (unless data strongly supports effect)
  expect_true(abs(treat_skeptical) <= abs(treat_neutral) + 0.5)
})

#-----------------------------------------------------------------------------
# Visualization Tests
#-----------------------------------------------------------------------------

test_that("Visualization functions work with synthetic data", {
  skip_if_not_installed("brms")
  skip_if_not_installed("ggplot2")
  skip_on_cran()

  # Generate and fit a simple model
  data <- generate_synthetic_data(
    n_observations = TEST_CONFIG$n_obs_quick,
    outcome_type = "continuous",
    seed = TEST_CONFIG$seed
  )

  model_spec <- create_model_spec("outcome_continuous", "treatment", family = "gaussian")

  priors <- PriorSpecification()
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
  priors <- add_prior(priors, "Intercept", "normal", 10, 10)
  priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)
  attr(model_spec, "priors") <- priors

  result <- fit_model(
    model_spec = model_spec,
    data = data,
    chains = 1, iter = TEST_CONFIG$n_iter, cores = 1, seed = TEST_CONFIG$seed
  )

  # Test posterior plot
  post_plot <- plot_posterior(result)
  expect_true(inherits(post_plot, "gg") || inherits(post_plot, "ggplot"))

  # Test forest plot
  estimates <- extract_estimates(result)
  forest_data <- tibble::tibble(
    trial = "Test",
    estimate = estimates$mean[estimates$parameter == "b_treatment"],
    ci_lower = estimates$`2.5%`[estimates$parameter == "b_treatment"],
    ci_upper = estimates$`97.5%`[estimates$parameter == "b_treatment"]
  )

  forest_plot <- plot_forest(forest_data, parameter = "treatment")
  expect_true(inherits(forest_plot, "gg") || inherits(forest_plot, "ggplot"))
})

#-----------------------------------------------------------------------------
# Export/Import Tests
#-----------------------------------------------------------------------------

test_that("Prior export and import works correctly", {
  priors <- PriorSpecification()
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
  priors <- add_prior(priors, "b_age", "normal", 0, 1.0)
  priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)

  temp_file <- tempfile(fileext = ".csv")

  # Export
  export_priors(priors, temp_file)
  expect_true(file.exists(temp_file))

  # Import
  priors_imported <- import_priors(temp_file)
  expect_s3_class(priors_imported, "PriorSpecification")
  expect_equal(nrow(priors_imported$priors), nrow(priors$priors))

  # Clean up
  unlink(temp_file)
})

test_that("Reproducibility with seed setting works", {
  # Generate data twice with same seed
  data1 <- generate_synthetic_data(n_observations = 50, seed = 123)
  data2 <- generate_synthetic_data(n_observations = 50, seed = 123)

  expect_equal(data1$outcome_continuous, data2$outcome_continuous)
  expect_equal(data1$treatment, data2$treatment)
  expect_equal(data1$age, data2$age)

  # Different seed should produce different data
  data3 <- generate_synthetic_data(n_observations = 50, seed = 456)
  expect_false(all(data1$outcome_continuous == data3$outcome_continuous))
})
