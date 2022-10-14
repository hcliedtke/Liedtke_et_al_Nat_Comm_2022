# README

This folder contains the scripts for running the corHMM() analysis performed by Liedtke et al. 2022.

The data is divided into three taxonomic datasets: Anura, Caudata, Gymnophiona (gymno). The procedure are similar for all three. They differ primarily in that not all life history modes are represented in each set and that in the case of Anura, the ancestral state reconstructions are performed with a model that contains hidden states.

Each taxonomic folder contains at least the following files/directories:

* `<taxon>.csv` : multi-trait coding for life history
* `<taxon>.tre` : taxonomic subset of the phylogeny
* `_run_corHMM_<taxon>.R` : main R script to run corHMM
* `_corHMM_models_<taxon>.R` : auxiliary R script, called on by the main script, for constructing the different transition models/hypotheses.
* `corHMM_<taxon>_q.pdf` : output from the auxiliary models script. A visualisation of all of the transition matrices ("q") tested.
* `corHMM_<taxon>.RData` : workspace image containing all data objects
* `muhisse_fit_<taxon>.txt` : exported model fitting results (tab separated table)
* `corHMM_output` : a directory containing all models results as separate .rds files
* `corHMM_fit_<taxon>.rds` : a list object with all model results combined
* `<taxon>_fit_summary.csv` : a table of model rankings.

Each taxonomic folder also contains two sub-directories for the ancestral state reconstructions on the best performing model:
* `anc_state_joint`: contains scripts and output for reconstructing ancestral states using joint estimations.
* `anc_state_simmap`: contains scripts and output for reconstructing ancestral states using stochastic character mapping.

A fourth folder "aux_scripts" contains scripts that are called by the main "_run" scripts. They contain custom functions used by the main scripts.


## Contact

All scripts were written by Christoph Liedtke (with lots of inspiration from many sources, including the package authors of the packages used). for any questions, contact me at:  `christoph.liedtke@ebd.csic.es`.
