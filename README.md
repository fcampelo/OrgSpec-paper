# OrgSpec paper
Full reproducibility scripts and data for the paper "Organism Specific Data Sets Improve Linear B-Cell Epitope Prediction Performance". 

The scripts contained in this folder make use of the [epitopes](https://github.com/fcampelo/epitopes/tree/v0.5.1-OrgSpec-paper) package, version 0.5.1, which can be installed into R using:

`devtools::install_github("fcampelo/epitopes", ref = "v0.5.1-OrgSpec-paper")`

Our results were generated with the following setup (taking advantage of some parallel processing capabilities of the _epitopes_ package):

> R version 4.0.5 (2021-03-31)

> Platform: x86_64-apple-darwin17.0 (64-bit)

> Running under: macOS Big Sur 10.16

*****

## How to re-run the experiment:

* Download or clone this repository locally.
* For each pathogen:
    1. Set the target organism folder (under directory _Experiments_) as the working directory.    
    2. Execute routine `01_generate_datasets`
    3. Generate predictions using the benchmark predictors (ABCPred, Bepipred2, etc.) for the hold-out proteins (under subfolder _data/splits_ of your working directory) and save them to the appropriate folders under subfolder _output_.
* Run routine `Check_leakage` (under _Experiments_).
* For each pathogen:
    1. Execute routine `02_run_experiment`
* Run routine `Consolidate_results` (under _Experiments_).

*****
