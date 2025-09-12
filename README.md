# BayesTrials: A Comprehensive Framework for Reproducible Bayesian Clinical Trial Analysis

[![R-CMD-check](https://github.com/AartikSarma/bayestrials/workflows/R-CMD-check/badge.svg)](https://github.com/AartikSarma/bayestrials/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> A comprehensive R package for reproducible Bayesian analyses of clinical datasets using brms

BayesTrials provides a standardized framework for conducting Bayesian clinical trial analyses with systematic prior specification, automated model fitting, and comprehensive reporting capabilities.

## Features

- **Standardized Prior Specification**: Tibble-based format for defining neutral, optimistic, and skeptical priors
- **Multiple Model Support**: Gaussian, Binomial, Poisson, and other distributions via brms
- **Parallel Processing**: Automated fitting of all prior × model permutations using `future`/`furrr`
- **Convergence Diagnostics**: Built-in Rhat, ESS, and divergence monitoring
- **Rich Visualizations**: Prior vs posterior plots, forest plots, and diagnostic dashboards
- **Automated Reporting**: Parameterized RMarkdown templates with HTML/PDF export
- **Model Comparison**: WAIC and LOO-IC for Bayesian model selection
- **CSV-based Workflow**: Define models and priors through simple CSV files
- **Data Harmonization**: Work with multiple trials with different variable names
- **Effect Modifier Identification**: Discover and validate treatment effect heterogeneity

## Installation

### System Requirements

Before installing bayestrials, ensure you have:

- **R** (≥ 4.0.0)
- **C++ compiler** for Stan compilation:
  - **macOS**: Install Xcode Command Line Tools: `xcode-select --install`
  - **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  - **Linux**: Install `g++` and related development tools

### Package Installation

```r
# Install development version from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("AartikSarma/bayestrials")

# Load the package
library(bayestrials)
```

### Key Dependencies

bayestrials uses [brms](https://paul-buerkner.github.io/brms/) (Bayesian Regression Models using Stan) as its modeling backend. This provides:
- Robust Stan model compilation and caching
- Comprehensive family support (Gaussian, Binomial, Poisson, etc.)
- Advanced model diagnostics and convergence checking

## Quick Start

```r
library(bayestrials)

# 1. Generate synthetic clinical trial data
set.seed(123)
data <- generate_synthetic_data(
  n_observations = 200,    # Start smaller for quicker testing
  outcome_type = "continuous",
  seed = 123
)

# 2. Create model specification
model_spec <- create_model_spec(
  outcome = "outcome_continuous",
  predictors = "treatment + age_centered + baseline_severity_centered",
  family = "gaussian",
  model_name = "linear_regression"
)

# 3. Specify priors (using the improved PriorSpecification system)
priors <- PriorSpecification()
priors <- add_prior(priors, "b_treatment", "normal", 0, 2.5)
priors <- add_prior(priors, "b_age_centered", "normal", 0, 1.0)
priors <- add_prior(priors, "b_baseline_severity_centered", "normal", 0, 2.0)
priors <- add_prior(priors, "Intercept", "normal", 0, 10)
priors <- add_prior(priors, "sigma", "half_cauchy", NA, 5, 0)  # Automatically converts to student_t

# Attach priors to model
attr(model_spec, "priors") <- priors

# 4. Set up parallel processing (optional)
setup_parallel(strategy = "sequential")  # Use sequential for initial testing

# 5. Fit Bayesian model
result <- fit_model(
  model_spec = model_spec,
  data = data,
  chains = 2,      # Fewer chains for quicker testing
  iter = 1000      # Fewer iterations for quicker testing
)

# 6. Check model and diagnostics
print(result)
diagnostics <- check_diagnostics(result)
print(diagnostics)

# 7. Extract and visualize results
estimates <- extract_estimates(result)
print(estimates)

# Create visualizations
plot_posterior(result)
forest_plot <- plot_forest(estimates)
```

### Expected Output

The model should converge successfully and show:
- Treatment effect estimate with credible intervals
- Convergence diagnostics (Rhat ≈ 1.0, adequate ESS)
- Posterior distributions for all parameters

## CSV-Based Workflow

```r
# Load and prepare data
data <- load_clinical_trial("trial_data.csv")
data <- preprocess_data(data, specs = "covariate_specs.csv")

# Specify model from CSV files
model_spec <- read_model_spec("model_specs.csv", 
                             "prior_specs.csv")

# Fit model
model <- fit_model(model_spec, data)

# Generate comprehensive report
generate_report(model, template = "basic_report")
```

## Multi-Trial Workflow

```r
# Load multiple trials
trials <- load_trials(c("trial1.csv", "trial2.csv", "trial3.csv"))

# Harmonize variables across trials
harmonized_trials <- harmonize_trials(trials, mapping = "variable_mapping.csv")

# Prepare data
prep_data <- preprocess_data(harmonized_trials, specs = "covariate_specs.csv")

# Specify and fit hierarchical model
model_spec <- read_model_spec("multi_trial_model.csv", 
                            "multi_trial_priors.csv")
model <- fit_model(model_spec, prep_data)

# Analyze treatment effects across trials
trial_effects <- extract_trial_effects(model, "treatment")
plot_forest(trial_effects)

# Test for effect modifiers
modifiers <- identify_effect_modifiers(model, prep_data, 
                                    treatment_var = "treatment",
                                    candidate_modifiers = c("age", "sex", "biomarker"))
```

## Sensitivity Analysis Workflow

```r
# Specify base model
base_spec <- read_model_spec("base_model.csv", "base_priors.csv")
base_model <- fit_model(base_spec, data)

# Define sensitivity analyses
sensitivity_specs <- list(
  prior_skeptical = update_priors(base_spec, skepticism = "high"),
  prior_optimistic = update_priors(base_spec, skepticism = "low"),
  missing_mar = update_model(base_spec, missing_mechanism = "mar"),
  missing_mnar = update_model(base_spec, missing_mechanism = "mnar")
)

# Run sensitivity analyses
sensitivity_results <- run_sensitivity_analyses(base_model, 
                                              data, 
                                              sensitivity_specs)

# Visualize sensitivity
plot_sensitivity(sensitivity_results, 
                parameter = "treatment_effect")

# Generate summary report
sensitivity_report(sensitivity_results, 
                  template = "sensitivity_report")
```

## Package Structure

```
bayestrials/
├── R/
│   ├── model_spec.R             # Model specification framework
│   ├── prior_spec.R             # Prior specification system  
│   ├── model_fitting.R          # brms-based model fitting
│   ├── visualization.R          # Plotting functions
│   ├── diagnostics.R            # Convergence diagnostics
│   ├── generate_synthetic_data.R # Test data generation
│   ├── standard_priors.R        # Standard prior sets
│   ├── enhanced_diagnostics.R   # Enhanced diagnostic tools
│   └── parallel.R               # Parallel processing utilities
├── inst/
│   ├── extdata/                 # Example CSV files
│   └── rmarkdown/               # Report templates
├── tests/
│   └── testthat/                # Unit tests
└── vignettes/                   # Documentation
```

## Dependencies

### Core Statistical
- [brms](https://paul-buerkner.github.io/brms/) - Bayesian regression models using Stan
- [posterior](https://mc-stan.org/posterior/) - Tools for working with posterior distributions
- [bridgesampling](https://github.com/quentingronau/bridgesampling) - Model comparison via bridge sampling

### Data Management
- [tidyverse](https://www.tidyverse.org/) - Data manipulation and visualization
- [tibble](https://tibble.tidyverse.org/) - Modern data frames
- [dplyr](https://dplyr.tidyverse.org/) - Data manipulation
- [tidyr](https://tidyr.tidyverse.org/) - Data tidying

### Parallel Processing  
- [future](https://future.futureverse.org/) - Parallel and distributed processing
- [furrr](https://furrr.futureverse.org/) - Apply mapping functions in parallel
- [progressr](https://progressr.futureverse.org/) - Progress reporting

### Visualization & Reporting
- [ggplot2](https://ggplot2.tidyverse.org/) - Grammar of graphics
- [plotly](https://plotly.com/r/) - Interactive visualizations
- [rmarkdown](https://rmarkdown.rstudio.com/) - Dynamic documents
- [gt](https://gt.rstudio.com/) - Publication-ready tables

## Development

### Setup
```r
# Install development dependencies
devtools::install_dev_deps()

# Load package for development
devtools::load_all()

# Generate documentation
devtools::document()
```

### Testing
```r
# Run all tests
devtools::test()

# Test the complete workflow
source("test_bayestrials_workflow.R")
```

### Package Checks
```r
# Run R CMD check
devtools::check()
```

## Troubleshooting

### Common Installation Issues

**Stan Compilation Errors:**
- **macOS**: Install Xcode Command Line Tools: `xcode-select --install`
- **Windows**: Ensure Rtools is properly installed and added to PATH
- **All platforms**: Restart R after installing compiler tools

**Memory Issues:**
- Reduce sample size for initial testing: `n_observations = 100`
- Use fewer chains: `chains = 1` or `chains = 2`
- Reduce iterations: `iter = 500`

### Common Model Fitting Issues

**Convergence Problems:**
- Increase `adapt_delta`: `control = list(adapt_delta = 0.99)`
- Increase iterations: `iter = 4000`
- Check for data issues (missing values, outliers)

**Prior Specification:**
- `half_cauchy` priors are automatically converted to `student_t(3, 0, scale)`
- Use `to_brms_prior()` to preview prior conversion before fitting

## Getting Help

- **Vignettes**: Run `browseVignettes("bayestrials")` for detailed tutorials
- **Function Help**: Use `?function_name` for specific function documentation  
- **Issues**: Report bugs and feature requests on [GitHub](https://github.com/AartikSarma/bayestrials/issues)
- **Stan Resources**: [Stan User's Guide](https://mc-stan.org/users/documentation/) for underlying methodology

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use `bayestrials` in your research, please cite:

```bibtex
@software{bayestrials,
  author = {Sarma, Aartik},
  title = {BayesTrials: A Comprehensive Framework for Reproducible Bayesian Clinical Trial Analysis},
  url = {https://github.com/AartikSarma/bayestrials},
  version = {0.1.0},
  year = {2024}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.