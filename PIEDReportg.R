---
  title: "PIED Report"
author: "Sharmila Dey, Emily Schultz"
date: "6/15/2020"
output: html_document
---
  
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
library(rstan)
options(mc.cores = parallel::detectCores())
library(parallel) 
library(mcmcplots) ; library(lattice) ; library(MASS)
library(lme4) ; library(nlme) ; library(splines); library(MCMCpack)
library(ggplot2);
library(tidyverse); library(plyr); 
data<-read.csv("data.clim copy.csv")
```

## Introduction to Stan Code

Project: To analyze local adaptation in Pinus edulis, we are using Bayesian statistics to understand changes in growth patterns in response to the changing climate.

Our Stan code uses PRISM climate data correllated with Pinus edulis data. This data is used in our "grow" function: 
  ```{r grow, echo TRUE}
grow<-na.omit(data) %>% 
  mutate_at(scale, .vars = vars(ppt_norm,ppt_yr,tmp_norm,tmp_yr,growth,solrad_an)) %>%
  arrange(PLOT.x,SUBP,name) %>%
  mutate(PlotCD=as.numeric(factor(PLOT.x, levels = unique(PLOT.x))),treeCD=as.numeric(factor(name,levels=unique(name))))
```

We then created the following objects: 
  
  - xG (predictor matrix)
- yG (log size at time t+1)
- nG (number of observations)
- plot (plot number)
- nplot (number of plots)
- K (number of predictors)
- tree (name of individual tree)
- ntree (number of unique trees)
- plotfortree (plot that each invidual tree comes from)

```{r variables, TRUE}
xG<-as.matrix(cbind(grow$ppt_norm, grow$tmp_norm, grow$ppt_yr, grow$tmp_yr, grow$DIA_prev))
##grow$ppt_norm*grow$tmp_norm,
# grow$ppt_norm*grow$tmp_yr, grow$ppt_norm*grow$ppt_yr, grow$ppt_norm*grow$DIA_prev,
#grow$ppt_yr*grow$tmp_yr, grow$ppt_yr*grow$tmp_norm, grow$ppt_yr*grow$DIA_prev, 
#grow$tmp_norm*grow$tmp_yr, grow$tmp_norm*grow$DIA_prev, grow$tmp_yr*grow$DIA_prev))
yG<-as.vector(grow$growth)
nG<-nrow(grow)
plot<-grow$PlotCD
nplot<-length(unique(grow$PlotCD))
K<-ncol(xG)
tree<-grow$treeCD
ntree<-length(unique(grow$treeCD))
plotfortree<-grow %>%
  group_by(treeCD) %>%
  summarize(Plot=mean(PlotCD))
plotfortree<-plotfortree$Plot
```

The "sink" and "cat" functions were used to initiate and manipulate the created variables. Each variable was initalized as seen below:
  ``` {r initialization, echo TRUE}
int<lower=0> K;         // N. predictors 
int<lower=0> nG;        // N. observations
matrix[nG,K] xG;        // Predictor matrix
vector[nG] yG;          // log size at time t+1 

int<lower=0> nplot;         // number of plots
int<lower=1> plot[nG];      // index for plot

int<lower=0> ntree;          // number of trees
int<lower=1> tree[nG];          // index for trees
int<lower=1> plotfortree[ntree]; //// plot index for each tree
```


## Including Plots

You can also embed plots, for example:
  
  ```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
