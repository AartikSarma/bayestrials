# bayestrials: A Robust Framework for Bayesian Clinical Trial Reanalysis

`bayestrials` is a comprehensive R package designed for Bayesian reanalysis of clinical trials, combining the statistical power of `brms` with the workflow efficiency of `tidyverse`. This package enables systematic model specification, flexible parameter definition, and powerful visualization capabilities while maintaining reproducibility.

## Features

- **CSV-based Model Specification**: Define models through simple CSV files
- **Data Harmonization**: Easily work with multiple trials with different variable names
- **Flexible Prior Specification**: Customize priors or use sensible defaults
- **Systematic Sensitivity Analysis**: Assess robustness to modeling assumptions
- **Effect Modifier Identification**: Discover and validate treatment effect heterogeneity
- **Publication-Ready Visualization**: Create high-quality plots for publications
- **Reproducibility Tools**: Ensure analyses can be reproduced exactly

## Installation

You can install the development version of bayestrials from GitHub:

```r
# install.packages("devtools")
devtools::install_github("AartikSarma/bayestrials")
```

## Basic Usage

```r
library(bayestrials)

# Load and prepare data
data <- load_clinical_trial("trial_data.csv")
data <- preprocess_data(data, specs = "covariate_specs.csv")

# Specify model
model_spec <- read_model_spec("model_specs.csv", 
                             "prior_specs.csv")

# Fit model
model <- fit_model(model_spec, data)

# Check diagnostics
check_diagnostics(model)

# Analyze results
summary(model)
plot_posterior(model, variables = c("b_treatment"))

# Generate report
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

## Documentation

For more details, check out the package vignettes:

```r
browseVignettes("bayestrials")
```

## Citation

If you use `bayestrials` in your research, please cite:

```
@Manual{bayestrials,
  title = {bayestrials: A Robust Framework for Bayesian Clinical Trial Reanalysis},
  author = {Aartik Sarma},
  year = {2025},
  note = {R package version 0.1.0},
  url = {https://github.com/AartikSarma/bayestrials},
}
```

## License

MIT