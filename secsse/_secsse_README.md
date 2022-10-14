# SecSSE

This folder contains the scripts for running the secsse_ml() analysis performed by Liedtke et al. 2022.

The data is divided into three taxonomic datasets: Anura, Caudata, Gymnophiona (gymno). The procedure are similar for all three. They differ primarily in that not all life history modes are represented in each set.

Each taxonomic folder contains at least the following files/directories:

* `_run_secsse_<taxon>.R` : main R script to run secsse
* `secsse_models_<taxon>.rds` : lists of starting parameters and model specifications
* `secsse_fit_summary.csv` : model fit results
* `best_model_div_rates.csv` : diversification rate estimates from best model
* `secsse_out` : folder containing all individual model results.

A fourth folder "aux_scripts" contains scripts that are called by the main "_run" scripts. They contain custom functions (`secsse_functions.R`) for preparing and processing the data, and a script that generates the transition models to be tested (`secsse_base_models.R`). Both are non-taxon specific.
  
Additionally the following global files are included:

* `<rep_modes.csv` : multi-trait coding for life history
* `amphibia.tre.tre` : the phylogeny


## Contact

All scripts were written by Christoph Liedtke (with lots of inspiration from many sources, including the package authors of the packages used). for any questions, contact me at:  `christoph.liedtke@ebd.csic.es`.
