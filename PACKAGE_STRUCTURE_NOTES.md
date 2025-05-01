# Package Structure Guidelines for bayestrials

## Correct Package Structure

The bayestrials package should maintain the following directory structure:

```
/bayestrials/             # Main package directory
├── R/                    # R code files
├── man/                  # Documentation files
├── inst/                 # Installed files
│   ├── extdata/          # External data for examples
│   └── ...
├── vignettes/            # Package vignettes
├── tests/                # Test files
├── DESCRIPTION           # Package description file
├── NAMESPACE             # Namespace file
├── LICENSE               # License file
├── README.md             # Readme file
└── .Rbuildignore         # Build ignore file
```

## Important Notes

1. **Avoid Nested Package Directories**: Do not create a nested directory structure like `/bayestrials/bayestrials/`. All package files should be directly in the main `/bayestrials/` directory.

2. **File Placement**: 
   - R code files go in the `R/` directory
   - Documentation files go in the `man/` directory
   - Example data files go in the `inst/extdata/` directory
   - Vignettes go in the `vignettes/` directory

3. **When Adding New Files**:
   - Always check the current directory structure
   - Make sure you're adding files to the correct level
   - Use absolute paths when necessary to avoid confusion

4. **Project File**:
   - The `.Rproj` file should be directly in the main package directory (`/bayestrials/bayestrials.Rproj`)
   - Do not create multiple project files in nested directories

## For Development

When developing the package:

1. Set your working directory to the main package directory (`/bayestrials/`)
2. Use `devtools::load_all()` to load the package for testing
3. Use `devtools::document()` to generate documentation
4. Use `devtools::build()` to build the package

## Checking Package Structure

You can verify your package structure is correct by running:

```r
setwd("/path/to/bayestrials")  # Go to main package directory
devtools::check()              # Run R CMD check
```

If the check returns errors about file paths or duplicate directories, review this guide to ensure your structure follows the recommended layout.