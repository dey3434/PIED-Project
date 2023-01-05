# Main Driver/Readme Script for running the PIED stan models:

# read in libraries:
# setting up libraries
library(here)
# install.packages("rstan", version = "2.19.3", repos = "http://cran.us.r-project.org")
library(rstan)
options(mc.cores = parallel::detectCores())
library(parallel) 
library(mcmcplots) ; library(lattice) ; library(MASS)
library(lme4) ; library(nlme) ; library(splines); library(MCMCpack)
library(ggplot2)
library(caret) ; library(tidyverse)
library(bayesplot)
library(here)
#library(gifski)
library(maps) 
library(tidyverse)

# use here to set wd
setwd(here())

# read in the PIED data:
source("Read_PIED_data_format.R")

# define and run stan model:
source("Run_PIED_stan_model.R")

# Make plots with output of stan model:
source("Plot_PIED_mod.R")