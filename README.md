# Set-Valued Policy Learning 

This repository provides an R package framework for Set-Valued Policy Learning, 
alongside all scripts necessary to reproduce the figures and results presented in the manuscript.

## Repository structure 
This repository follows the R package structure to handle dependencies, 
but it is primarily a reproduction archive.

* `R/`: Core functions for set-valued policy learning.
* `man/`: Auto-generated documentation for the package functions.
* `DESCRIPTION`: Metadata and R package dependencies.
* `inst/`: Reproducibility folder 
  * `toy_example/`: Scripts for synthetic data experiments and figures. 
  * `ivf_examples/`: Scripts for the IVF data application.
  * `main.R`: The entry point for real-world data applications.
  * `images/`: (Generated) Exports of paper figures.
  * `predictions/`: (Generated) Model outputs and intermediate results.


## Getting started (anonymous version)  

**1. Download the Source Code**

Download the repository as a ZIP file from
https://anonymous.4open.science/r/setValuedPolicyLearning-6B9C/
and extract it to your local machine.

**2. Install Dependencies and Load Package**

Open R or RStudio and run the following to set up the environment:
```
if (!require("devtools")) install.packages("devtools")

# Set working directory to the extracted folder
setwd("path/to/setValuedPolicyLearning-6B9C")

# Install required dependencies and load the package
devtools::install_deps(upgrade = "never")
devtools::load_all()
```

## How to reproduce synthetic setting figures 

The scripts are organized by experiment type. Running these will automatically 
generate the `inst/images/` and `inst/predictions/` folders.

```
# Execute the main toy example script
source("inst/toy_example/main_synthetic.R")
```