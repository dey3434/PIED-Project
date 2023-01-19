### Stan models for PIAL growth using technique in Ogle et al., Ecology Letters to test for climate effects
## Sharmila Dey
# 22 June 2020
# setwd("/home/work/")
#load(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_grow_coef2.rda"))
#fit_grow <- readRDS(url("https://de.cyverse.org/dl/d/888FD6F6-DAEA-46AE-BDF4-036A708990CC/log_normal_monsoonoos_pptoos.RDS"))
library(rstan)
options(mc.cores = parallel::detectCores())
library(parallel) 
library(mcmcplots) ; library(lattice) ; library(MASS)
library(lme4) ; library(nlme) ; library(splines); library(MCMCpack)
library(ggplot2)
library(caret) ;
library(tidyverse)
library(bayesplot)
library(here)
library(gifski)
library(maps)
####
library(patchwork)
library(sp)

# PIED.all <- read.csv(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_all_growth_v7.csv"))
# full.ppt.tmean.norms <- read.csv(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_all_tmean_ppt_v6.csv"))

if (file.exists(here::here("data", "PIED_data.csv"))) {
  PIED.all <- read.csv(here::here("data", "PIED_data.csv"))
} else {
  PIED.all <- read.csv(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_all_growth_v7.csv"))  
  write.csv(PIED.all, file = here::here("data", "PIED_data.csv"), row.names = FALSE)
}

if (file.exists(here::here("data", "climate_data.csv"))) {
  full.ppt.tmean.norms <- read.csv(here::here("data", "climate_data.csv"))
} else {
  full.ppt.tmean.norms <- read.csv(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_all_tmean_ppt_v6.csv"))
  write.csv(full.ppt.tmean.norms, file = here::here("data", "climate_data.csv"), row.names = FALSE)
}

grow.new <- left_join(PIED.all, full.ppt.tmean.norms, by.x = c("name", "year","LAT", "LON"),by.y = c("name", "year","lat", "lon"))


#grow.new <- read.csv("data/pied_all_tmean_ppt_v3.csv")

#grow.new <- merge(grow, newclimate, by.x = c("LON", "LAT", "name", "year"), by.y = c("lon", "lat", "name", "year"))

grow.monsoon <- na.omit(grow.new) %>% 
  mutate_at(scale, .vars = vars(Precip_JulAug, Precip_NovDecJanFebMar, tmp_norm, ppt_norm, tmp_yr)) %>%
  arrange(PLOT,SUBP,name) %>%
  mutate(PlotCD=as.numeric(factor(ST_PLT, levels = unique(ST_PLT))),treeCD=as.numeric(factor(name,levels=unique(name))),
         growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2))

set.seed(2023)
split=0.20
trainIndex <- createDataPartition(grow.monsoon$name, p=split, list=FALSE)
grow_test <- grow.monsoon[trainIndex,]
grow_train <- grow.monsoon[-trainIndex,]



xG<-as.matrix(cbind(grow_train$ppt_norm, grow_train$tmp_norm, #grow_train$Precip_JulAug, grow_train$Precip_NovDecJanFebMar,
                    grow_train$tmp_yr, grow_train$DIA_prev,
                    grow_train$ppt_norm*grow_train$tmp_norm, 
                    grow_train$ppt_norm*grow_train$tmp_yr, 
                    #grow_train$ppt_norm*grow_train$Precip_JulAug, 
                    grow_train$ppt_norm*grow_train$DIA_prev,
                    #grow_train$Precip_JulAug*grow_train$tmp_yr, grow_train$Precip_JulAug*grow_train$tmp_norm, grow_train$Precip_JulAug*grow_train$DIA_prev, 
                    grow_train$tmp_norm*grow_train$tmp_yr, 
                    grow_train$tmp_norm*grow_train$DIA_prev,
                    grow_train$tmp_yr*grow_train$DIA_prev))#, #grow_train$Precip_NovDecJanFebMar*grow_train$ppt_norm,
                    #grow_train$Precip_NovDecJanFebMar*grow_train$tmp_norm, grow_train$Precip_NovDecJanFebMar*grow_train$Precip_JulAug,
                    #grow_train$Precip_NovDecJanFebMar*grow_train$tmp_yr, grow_train$Precip_NovDecJanFebMar*grow_train$DIA_prev))
xGtest<-as.matrix(cbind(grow_test$ppt_norm, grow_test$tmp_norm, #grow_test$Precip_JulAug, grow_test$Precip_NovDecJanFebMar,
                        grow_test$tmp_yr, grow_test$DIA_prev,
                        grow_test$ppt_norm*grow_test$tmp_norm, 
                        grow_test$ppt_norm*grow_test$tmp_yr, 
                        #grow_test$ppt_norm*grow_test$Precip_JulAug, 
                        grow_test$ppt_norm*grow_test$DIA_prev,
                        #grow_test$Precip_JulAug*grow_test$tmp_yr, grow_test$Precip_JulAug*grow_test$tmp_norm, grow_test$Precip_JulAug*grow_test$DIA_prev, 
                        grow_test$tmp_norm*grow_test$tmp_yr, grow_test$tmp_norm*grow_test$DIA_prev, 
                        grow_test$tmp_yr*grow_test$DIA_prev)) #, #grow_test$Precip_NovDecJanFebMar*grow_test$ppt_norm,
                        #grow_test$Precip_NovDecJanFebMar*grow_test$tmp_norm, #grow_test$Precip_NovDecJanFebMar*grow_test$Precip_JulAug,
                        #grow_test$Precip_NovDecJanFebMar*grow_test$tmp_yr, grow_test$Precip_NovDecJanFebMar*grow_test$DIA_prev))
yG<-as.vector(grow_train$loggrowth)
yGtest<-as.vector(grow_test$loggrowth)
nG<-nrow(grow_train)
nGtest<-nrow(grow_test)
plot<-grow_train$PlotCD
nplot<-length(unique(grow_train$PlotCD))
K<-ncol(xG)
tree<-grow_train$treeCD
treetest<-grow_test$treeCD
ntree<-length(unique(grow_train$treeCD))
plotfortree<-grow_train %>%
  group_by(treeCD) %>%
  summarize(Plot=mean(PlotCD))
plotfortree<-plotfortree$Plot


sink("model_0.stan")
cat("
    data {
    
    int<lower=0> K;         // N. predictors 
    int<lower=0> nG;        // N. observations
    int<lower=0> nGtest;    // N. observations (ppc)
    matrix[nG,K] xG;        // Predictor matrix
    matrix[nGtest, K] xGtest;   // Predictor matrix (ppc)
    vector[nG] yG;          // log size at time t+1 
    
    int<lower=0> nplot;         // number of plots
    int<lower=1> plot[nG];      // index for plot
    
    int<lower=0> ntree;          // number of trees
    int<lower=1> tree[nG];          // index for trees
    int<lower=1> treetest[nGtest];  //index for trees (ppc)
    int<lower=1> plotfortree[ntree];   // plot index for each tree
    }
    
    parameters {
    
    real u_beta0;                          // intercept means
    vector[K] u_beta;                      // other coeff mean
    
    real beta0_p_tilde[nplot];                   // plot-level intercepts
    real<lower=0> s_beta0_p;               // plot-level intercept variance
    real beta0_t_tilde[ntree];                   // tree-level intercepts
    real<lower=0> s_beta0_t;               // tree-level intercept variance
    
    real<lower=0> sigma_y;                 // Residual for growth model
    
    }
    
    transformed parameters {
    real beta0_p[nplot];                   // plot-level intercepts
    real beta0_t[ntree];                   // tree-level intercepts
    
    for(p in 1:nplot){
    beta0_p[p] = u_beta0 + s_beta0_p * beta0_p_tilde[p];
    }
    
    for(t in 1:ntree){
    beta0_t[t] = beta0_p[plotfortree[t]] + s_beta0_t * beta0_t_tilde[t];
    }
    
    }
    
    model {
    vector[nG] mG;
    
    u_beta0 ~ normal(0, 100);
    beta0_p_tilde ~ normal(0,1);
    beta0_t_tilde ~ normal(0,1);
    
    u_beta ~ normal(0, 100); 
    s_beta0_p ~ cauchy(0,2.5);
    s_beta0_t ~ cauchy(0,2.5);
    
    sigma_y ~ gamma(2,0.01);
    
    //tried nesting random effects of trees within plots--had issue
    //for(p in 1:nplot){
    //beta0_p[p] ~ normal(u_beta0, s_beta0_p);
    //for(t in 1:ntree){
    //beta0_t[t] ~ normal(beta0_p[plotfortree[t]], s_beta0_t);
    //}
    //}
    
    
    // GROWTH MODEL
    
    for(n in 1:nG){
    mG[n] = beta0_t[tree[n]]+xG[n]*u_beta;
    }
    
    yG ~ normal(mG,sigma_y);
    //yG ~ gamma(mG,sigma_y);

    }
    
    generated quantities{
    vector[nGtest] yrep;
    for(n in 1:nGtest){
    yrep[n] = normal_rng(beta0_t[treetest[n]]+xGtest[n]*u_beta,sigma_y);
    }

    //for(n in 1:nGtest){
    //yrep[n] = gamma_rng(beta0_t[treetest[n]]+xGtest[n]*u_beta,sigma_y);
    //}

    }
    
    ",fill=T)

sink()


pied_dat <- list(K = K, nG = nG, nGtest = nGtest, yG = yG, xG = xG, xGtest = xGtest, plot = plot, 
                 nplot = nplot, tree = tree, treetest = treetest, ntree = ntree, 
                 plotfortree = plotfortree)
#tranG = tranG, SiteForTran = SiteForTran, Ntran_totalG=Ntran_totalG)
#indG = indG, TranForInd = TranForInd, Nind_totalG = Nind_totalG)



csvfiles <- here::here("results", paste0("log_normal_null_", 1:3, ".csv"))

if (all(file.exists(csvfiles))) {
  fit_grow <- read_stan_csv(csvfiles, col_major = TRUE) 
} else {
  fit_grow <- stan(file = 'model_0.stan', data = pied_dat, 
                   iter = 5000,
                   warmup = 1000,
                   chains = 3, cores = 8, 
                   sample_file = here::here("results", "log_normal_null"))
}
