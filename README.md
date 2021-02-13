# OrgSpec paper
Full reproducibility scripts and data for paper "Organism Specific Data Sets Improve Linear B-Cell Epitope Prediction Performance". 

The scripts contained in this folder make use of the [epitopes](https://github.com/fcampelo/epitopes/tree/v0.4.11-OrgSpec-paper) package, version 0.4.11, which can be installed into R using `devtools::install_github()`. 

We have tried to removed all explicit multicore parallelization from these scripts. This makes replication slower, but more robust to different platforms that you may be using. Our results were generated with the following setup (taking advantage of some parallel processing capabilities of the _epitopes_ package):

> R version 4.0.3 (2020-10-10)

> Platform: x86_64-apple-darwin17.0 (64-bit)

> Running under: macOS Big Sur 10.16

