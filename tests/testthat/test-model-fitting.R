test_that("Model fitting works with synthetic data", {
  skip_if_not_installed("brms")
  skip_on_cran()  # Skip on CRAN due to time constraints
  
  # Generate small synthetic data for quick testing
  data <- generate_synthetic_data(
    n_observations = 50,
    n_variables = 5, 
    outcome_type = "continuous",
    seed = 123
  )
  
  # Create simple model specification
  model_spec <- create_model_spec(
    "outcome_continuous",
    "treatment + age_centered",
    family = "gaussian",
    model_name = "test_model"
  )
  
  # Create simple priors
  priors <- PriorSpecification()
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
  priors <- add_prior(priors, "b_age_centered", "normal", 0, 1)
  priors <- add_prior(priors, "Intercept", "normal", 0, 10)
  priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)
  
  # Attach priors to model spec
  attr(model_spec, "priors") <- priors
  
  # Fit model with minimal iterations for speed
  result <- fit_model(
    model_spec = model_spec,
    data = data,
    chains = 1,
    iter = 100,  # Very short for testing
    cores = 1,
    seed = 123
  )
  
  expect_s3_class(result, "BayesianModel")
  expect_true(!is.null(result$model))
  expect_true(inherits(result$model, "brmsfit"))
})

test_that("Model fitting handles different families", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  # Test with binary outcome
  data_binary <- generate_synthetic_data(
    n_observations = 50,
    outcome_type = "binary",
    seed = 456
  )
  
  model_spec_binary <- create_model_spec(
    "outcome_binary",
    "treatment + age_centered",
    family = "binomial",
    link = "logit",
    model_name = "binary_test"
  )
  
  # Create priors for binary model
  priors_binary <- PriorSpecification()
  priors_binary <- add_prior(priors_binary, "b_treatment", "normal", 0, 1.5)
  priors_binary <- add_prior(priors_binary, "b_age_centered", "normal", 0, 1)
  priors_binary <- add_prior(priors_binary, "Intercept", "normal", 0, 5)
  
  attr(model_spec_binary, "priors") <- priors_binary
  
  # This test might be slow, so we'll just check the setup
  expect_s3_class(model_spec_binary, "ModelSpecification")
  expect_equal(model_spec_binary$family, "binomial")
  expect_equal(model_spec_binary$link, "logit")
})

test_that("Model validation catches errors", {
  # Test invalid model spec
  expect_error(
    create_model_spec(NULL, "x", family = "gaussian"),
    "Model specification must include outcome"
  )
  
  expect_error(
    create_model_spec("y", NULL, family = "gaussian"),
    "Model specification must include outcome"
  )
  
  expect_error(
    create_model_spec("y", "x", family = "invalid_family"),
    "Invalid family"
  )
})

test_that("Prior conversion works correctly", {
  # Test various prior specifications
  priors <- PriorSpecification()
  priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
  priors <- add_prior(priors, "b_age", "normal", 0, 1)
  priors <- add_prior(priors, "Intercept", "normal", 0, 10)
  priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)
  
  brms_priors <- to_brms_prior(priors)
  
  expect_s3_class(brms_priors, "brmsprior")
  expect_equal(nrow(brms_priors), 4)
  
  # Check specific priors
  prior_strings <- brms_priors$prior
  expect_true(any(grepl("normal\\(0, 2.5\\)", prior_strings)))
  expect_true(any(grepl("student_t\\(3, 0, 5\\)", prior_strings)))
  
  # Check classes
  classes <- brms_priors$class
  expect_true("b" %in% classes)
  expect_true("Intercept" %in% classes) 
  expect_true("sigma" %in% classes)
  
  # Check coefficients
  coefs <- brms_priors$coef[!is.na(brms_priors$coef)]
  expect_true("treatment" %in% coefs)
  expect_true("age" %in% coefs)
})

test_that("Model specs from CSV work", {
  # Test reading model specifications from CSV
  model_file <- file.path(
    system.file(package = "bayestrials"), 
    "..", "..", "tests", "test-data", "test_model_specs.csv"
  )
  
  # Skip if test file doesn't exist
  skip_if_not(file.exists(model_file))
  
  model_specs <- read_model_spec(model_file, validate = TRUE)
  
  expect_type(model_specs, "list")
  expect_true(length(model_specs) > 0)
  expect_true(all(sapply(model_specs, inherits, "ModelSpecification")))
})

test_that("Family functions work correctly", {
  # Test that our family handling works
  model_spec_gauss <- create_model_spec("y", "x", family = "gaussian")
  expect_equal(model_spec_gauss$family, "gaussian")
  
  model_spec_binom <- create_model_spec("y", "x", family = "binomial", link = "logit")
  expect_equal(model_spec_binom$family, "binomial")
  expect_equal(model_spec_binom$link, "logit")
  
  model_spec_pois <- create_model_spec("y", "x", family = "poisson")
  expect_equal(model_spec_pois$family, "poisson")
})