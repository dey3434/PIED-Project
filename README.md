# PIED-Project
## About this repository
This repository contains code used to estimate Bayesian models of _Pinus edulis_ tree growth from tree cores sampled within Forest Inventory and Analysis inventory plots. Growth models developed here were used to estimate the sensitivity of _P. edulis_ to climate drivers across the entire species range. 

## Motivation
Species responses to climate changes are often modelled with species occurance data, which may underestimated responses to annual climate variability. Here we use Bayesian regression models to estimate the sensitivity of _P. edulis_ growth to spatial and time-varying climate across its core species range using 1558 tree ring time series.

## Organization of the Repository
This repository contains code that does the  model estimation, validation, and model comparison across 10 different Bayesian regression model forms for annual tree ring growth observations (for more information about the motivation for each model see the associated manuscript). Bayesian models are fit using STAN and rstan.

`run-models.R` will run all 10 growth models contained in the scripts labeled `model_[n].R` where [n] indicates models 0-9 and generate outputs. Note that the model numbering in the repository is differnt than that presented in the manuscript text, which we switched up the model numbering scheme for the sake of clarity. The models get renamed in `compare-models.R` , which summarises the model output for presentation in the manuscript. 

To recreate these model estimates, you would need the associated tree ring data and metadata, which can be found at this DOI:
Additional data used in this analysis includes PRISM 4km gridded data product; see Get_PRISM_Data.R code for PRISM extraction code.

We also fit probability of occurance models using USDA Forest Inventory and Analysis (FIA) observations (https://www.fs.usda.gov/research/products/dataandtools/datasets) of occurance and PRISM climate normals. These models are fit and outputs are generated in the `pOccurance` folder.

## Citation
If using the data or code please use the following citation

## Funding sources

