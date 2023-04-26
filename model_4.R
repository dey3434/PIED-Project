### Stan models for PIAL growth using technique in Ogle et al., Ecology Letters to test for climate effects
## Sharmila Dey
# 22 June 2020
#setwd("/home/rstudio")
#load(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_grow_coef2.rda"))
#fit_grow <- readRDS(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/ppt_tmp_springfall_sizefix_scale_3.RDS"))
#fit_grow <- readRDS("ppt_tmp_springfall_sizefix.RDS")
#fit_grow_old <- readRDS("data/ppt_tmp_springfall.RDS")
#install.packages("rstan", version = "2.19.3", repos = "http://cran.us.r-project.org")
library(rstan)
options(mc.cores = parallel::detectCores())
library(parallel) 
library(mcmcplots) ; 
library(lattice) ; 
library(MASS)
library(lme4) ; 
library(nlme) ; 
library(splines); library(MCMCpack)
library(ggplot2)
library(caret) ; 
library(tidyverse)
library(bayesplot)
library(here)
library(gifski)
library(maps)
library(tidyr)
library(dplyr)
####
library(patchwork)
library(sp)
library(loo)
library(rstantools)


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

hist(full.ppt.tmean.norms$tmp_norm)
hist(grow.new$tmp_norm)
yearlyprecip <- as.matrix(cbind(grow.new$Precip_JulAug, grow.new$Precip_NovDecJanFebMar))
colnames(yearlyprecip) <- c("Precip JulAug", "Precip NovDecJanFebMar")
boxplot(yearlyprecip, main="Seasonal Precipitation")  

yearlytemp <- as.matrix(cbind(grow.new$Tmean_AprMayJun, grow.new$Tmean_SepOct))
colnames(yearlytemp) <- c("Tmean AprMayJun", "Tmean SepOct")
boxplot(yearlytemp, main="Seasonal Temperature")  



grow.new$Tmean_AprMayJun_l <- grow.new$Tmean_AprMayJun
grow.new$Tmean_SepOct_l <- grow.new$Tmean_SepOct
grow.new$Precip_NovDecJanFebMar_l <- grow.new$Precip_NovDecJanFebMar
grow.new$Precip_JulAug_l <- grow.new$Precip_JulAug

grow.monsoon<-na.omit(grow.new) %>% 
  mutate_at(scale, .vars = vars(tmp_norm, ppt_norm)) %>%
  
  arrange(PLOT,SUBP,name) %>%
  mutate(PlotCD=as.numeric(factor(ST_PLT, levels = unique(ST_PLT))),treeCD=as.numeric(factor(name,levels=unique(name))),
         growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2)) %>%
  group_by(PlotCD)%>%
  mutate_at(scale, .vars = vars(Precip_JulAug, Precip_NovDecJanFebMar, Tmean_AprMayJun, Tmean_SepOct)) %>%
  ungroup()


# the other way of scaling:
# grow.monsoon.old <-na.omit(grow.new) %>% 
#   mutate_at(scale, .vars = vars(Precip_JulAug, Precip_NovDecJanFebMar, Tmean_AprMayJun, Tmean_SepOct, tmp_norm, ppt_norm)) %>%
#   
#   arrange(PLOT,SUBP,name) %>%
#   mutate(PlotCD=as.numeric(factor(PLOT, levels = unique(PLOT))),treeCD=as.numeric(factor(name,levels=unique(name))),
#          growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2)) 



set.seed(2023)
split=0.20
trainIndex <- createDataPartition(grow.monsoon$name, p=split, list=FALSE)
grow_test <- grow.monsoon[trainIndex,]
grow_train <- grow.monsoon[-trainIndex,]



xG <- as.matrix(cbind(grow_train$ppt_norm, grow_train$tmp_norm, grow_train$ppt_norm*grow_train$tmp_norm,
                      grow_train$Precip_JulAug, grow_train$Precip_NovDecJanFebMar, 
                      grow_train$Tmean_AprMayJun, grow_train$Tmean_SepOct, grow_train$DIA_prev, 
                      grow_train$DIA_prev*grow_train$ppt_norm, grow_train$DIA_prev*grow_train$tmp_norm,
                      grow_train$DIA_prev*grow_train$Precip_JulAug, grow_train$DIA_prev*grow_train$Precip_NovDecJanFebMar, 
                      grow_train$DIA_prev*grow_train$Tmean_AprMayJun, grow_train$DIA_prev*grow_train$Tmean_SepOct,
                      grow_train$ppt_norm*grow_train$Precip_JulAug, grow_train$ppt_norm*grow_train$Precip_NovDecJanFebMar,
                      grow_train$ppt_norm*grow_train$Tmean_AprMayJun, grow_train$ppt_norm*grow_train$Tmean_SepOct,
                      grow_train$tmp_norm*grow_train$Precip_JulAug, grow_train$tmp_norm*grow_train$Precip_NovDecJanFebMar, 
                      grow_train$tmp_norm*grow_train$Tmean_AprMayJun, grow_train$tmp_norm*grow_train$Tmean_SepOct, 
                      grow_train$Precip_JulAug*grow_train$Precip_NovDecJanFebMar, grow_train$Precip_JulAug*grow_train$Tmean_AprMayJun, 
                      grow_train$Precip_JulAug*grow_train$Tmean_SepOct, grow_train$Precip_NovDecJanFebMar*grow_train$Tmean_AprMayJun, 
                      grow_train$Precip_NovDecJanFebMar*grow_train$Tmean_SepOct, grow_train$Tmean_AprMayJun*grow_train$Tmean_SepOct))
xGtest<-as.matrix(cbind(grow_test$ppt_norm, grow_test$tmp_norm, grow_test$ppt_norm*grow_test$tmp_norm,
                        grow_test$Precip_JulAug, grow_test$Precip_NovDecJanFebMar, 
                        grow_test$Tmean_AprMayJun, grow_test$Tmean_SepOct, grow_test$DIA_prev, 
                        grow_test$DIA_prev*grow_test$ppt_norm, grow_test$DIA_prev*grow_test$tmp_norm,
                        grow_test$DIA_prev*grow_test$Precip_JulAug, grow_test$DIA_prev*grow_test$Precip_NovDecJanFebMar, 
                        grow_test$DIA_prev*grow_test$Tmean_AprMayJun, grow_test$DIA_prev*grow_test$Tmean_SepOct,
                        grow_test$ppt_norm*grow_test$Precip_JulAug, grow_test$ppt_norm*grow_test$Precip_NovDecJanFebMar,
                        grow_test$ppt_norm*grow_test$Tmean_AprMayJun, grow_test$ppt_norm*grow_test$Tmean_SepOct,
                        grow_test$tmp_norm*grow_test$Precip_JulAug, grow_test$tmp_norm*grow_test$Precip_NovDecJanFebMar, 
                        grow_test$tmp_norm*grow_test$Tmean_AprMayJun, grow_test$tmp_norm*grow_test$Tmean_SepOct, 
                        grow_test$Precip_JulAug*grow_test$Precip_NovDecJanFebMar, grow_test$Precip_JulAug*grow_test$Tmean_AprMayJun, 
                        grow_test$Precip_JulAug*grow_test$Tmean_SepOct, grow_test$Precip_NovDecJanFebMar*grow_test$Tmean_AprMayJun, 
                        grow_test$Precip_NovDecJanFebMar*grow_test$Tmean_SepOct, grow_test$Tmean_AprMayJun*grow_test$Tmean_SepOct))
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


sink("model_4.stan")
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



csvfiles <- here::here("results", paste0("ppt_tmp_springfall_sizefix_scale_small_", 1:3, ".csv"))

if (all(file.exists(csvfiles))) {
  fit_grow <- read_stan_csv(csvfiles, col_major = TRUE) 
} else {
  fit_grow <- stan(file = 'model_4.stan', data = pied_dat, 
                   iter = 5000,
                   warmup = 1000,
                   chains = 3, cores = 8, 
                   sample_file = here::here("results", "ppt_tmp_springfall_sizefix_scale_small"))
}

# chain1 <- rstan::read_stan_csv("/home/rstudio/ppt_tmp_springfall_sizefix_scale_small_1.csv")
# chain2 <- rstan::read_stan_csv("/home/rstudio/ppt_tmp_springfall_sizefix_scale_small_2.csv")
# chain3 <- rstan::read_stan_csv("/home/rstudio/ppt_tmp_springfall_sizefix_scale_small_3.csv")

# fit_grow <- rstan::read_stan_csv(csvfiles)
# saveRDS(fit_grow, file = "ppt_tmp_springfall_sizefix_scale.RDS")


#summary<-summary(fit_grow)
#summary

# Updated code below -----
fit_grow_df <- as.data.frame(fit_grow)
plotdata<-select(fit_grow_df,"yrep[1]":"yrep[8780]")
plotdatainterval<-select(fit_grow_df, "u_beta[1]":paste0("u_beta[", ncol(xG), "]"))
# plotdatainterval<-dplyr::select(fit_grow_df, "u_beta[1]":"u_beta[28]")
colnames(plotdatainterval) <- c("MAP", "MAT","MAP*MAT", "monsoon precip", "winter precip", 
                                "spring temp", "fall temp", "tree size",
                                "tree size*MAP", "tree size*MAT", "tree size*monsoon precip",
                                "tree size*winter precip", "tree size*spring temp",
                                "tree size*fall temp", "MAP*monsoon precip", 
                                "MAP*winter precip","MAP*spring temp", "MAP*fall temp",
                                "MAT*monsoon precip", "MAT*winter precip", "MAT*spring temp",
                                "MAT*fall temp",  "monsoon precip*winter precip",
                                "monsoon precip*spring temp", "monsoon precip*fall temp",
                                "winter precip*spring temp", "winter precip*fall temp",
                                "spring temp*fall temp")

# get summaries of plotdatainterval:
if (file.exists(here::here("results", "betasummaries_model4_threechain_PIED.csv"))) {
  beta.summaries <- read.csv(here::here("results", "betasummaries_model4_threechain_PIED.csv"))
} else {
  df <- reshape2::melt(plotdatainterval)
  beta.summaries <- df %>% group_by(variable) %>%
    summarise(mean = mean(value),
              ci.lo = quantile(value, 0.025),
              ci.hi = quantile(value, 0.975))
  
  beta.summaries$allpos <- ifelse(beta.summaries$mean > 0 &  beta.summaries$ci.lo > 0 &  beta.summaries$ci.hi > 0, "yes", "no")
  beta.summaries$allneg <- ifelse(beta.summaries$mean <= 0 &  beta.summaries$ci.lo <= 0 &  beta.summaries$ci.hi <= 0, "yes", "no")
  beta.summaries$significant <- ifelse(beta.summaries$allpos == "yes" | beta.summaries$allneg == "yes", "significant", "not significant")
  write.csv(beta.summaries, here::here("results", "betasummaries_model4_threechain_PIED.csv"), row.names = FALSE)
}


####
#Validation-- This will be in a separate File
ext_fit <- rstan::extract(fit_grow)
yrep <- ext_fit$yrep
#yrep <- exp(yrep)
mean.pred <- apply(ext_fit$yrep, 2, mean)
p.o.df <- data.frame(predicted = exp(mean.pred), observed = exp(grow_test$loggrowth), error = (exp(mean.pred) - exp(grow_test$loggrowth)))
meansqrd <- (mean(p.o.df$error))^2

##### 
save(p.o.df, file = here::here("results", "model-4-pred-obs.RData"))

# ggplot(p.o.df, aes(predicted, observed)) + geom_point(alpha = 0.1) + geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
#   ylim(0, 10) + xlim(0,10)

####
p_pred_vs_observed <- ggplot(p.o.df, aes(predicted, observed)) + 
  geom_point(alpha = 0.1) + 
  geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
  ylim(0, 10) + xlim(0,10)
p_pred_vs_observed
ggsave(here::here("images", "model_4", "pred_vs_observed.png"), p_pred_vs_observed)


# Updated code below ------
sigma <- fit_grow_df[,"sigma_y"]
# Plot-level random effects (not included in the current version) ----
beta_0ps <- select(fit_grow_df, "beta0_p[1]":paste0("beta0_p[", nplot, "]"))
# Tree-level random effects ------
beta_0ts <- select(fit_grow_df, "beta0_t[1]":paste0("beta0_t[", ntree, "]"))
# Modify mu to include the tree-level and plot-level random effects ------
mu <- as.matrix(plotdatainterval) %*% t(xG) + as.matrix(beta_0ts[, tree])
# mu <- as.matrix(plotdatainterval) %*% t(xG)
ll <- matrix(0, length(sigma), length(yG))
for(i in 1:length(sigma)){
  ll[i,] <- dnorm(yG, mu[i,], sd = sigma[i], log = TRUE)
}
newll <- as.matrix(ll)
r_eff <- relative_eff(exp(ll), chain_id = rep(1:3, each = 4000), cores = 1)
leaveoneout <- loo::loo(as.matrix(ll), r_eff = r_eff, save_psis = TRUE, cores = 1)

save(ll, r_eff, leaveoneout, file = here::here("results", "model-4-loo.RData"))

yrep <- matrix(0, length(sigma), length(yG))
for(i in 1:length(sigma)){
  # Modified this to be correct ------
  yrep[i,] <- rnorm(length(yG), mu[i,], sd = sigma[i])
}
psis <- leaveoneout$psis_object
keep_obs <- sample(1:length(yG), 100)
lw <- weights(psis)
ppc1 <- ppc_loo_intervals(yG, yrep = yrep, psis_object = psis, subset = keep_obs, order = "median") 
ppc2 <- ppc_loo_pit_overlay(yG, yrep = yrep, lw = lw)
ppc3 <- ppc_loo_pit_qq(yG, yrep = yrep, lw = lw)

ggsave(here::here("images", "model_4", "ppc-plot-1.png"),
       ppc1, width = 16/3, height = 9)
ggsave(here::here("images", "model_4", "ppc-plot-2.png"),
       ppc2, width = 16/3, height = 9)
ggsave(here::here("images", "model_4", "ppc-plot-3.png"),
       ppc3, width = 16/3, height = 9)

# Commented out 03/10/2023 by JRT
waic(ll)


# yrep_mean <- apply(yrep, MARGIN = 2, FUN = mean)
# yrep_ci.low <- apply(yrep, MARGIN = 2, FUN = function(x){quantile(x, 0.025)})
# yrep_ci.high <- apply(yrep, MARGIN = 2, FUN = function(x){quantile(x, 0.975)})
# 
# p.o.df <- data.frame(ci.low = yrep_ci.low, mean = yrep_mean, ci.high = yrep_ci.high, observed = yG,
#                      ppt_norm = grow_train$ppt_norm, tmp_norm = grow_train$tmp_norm, Precip_JulAug = grow_train$Precip_JulAug,
#                      Precip_NovDecJanFebMar = grow_train$Precip_NovDecJanFebMar, Tmean_AprMayJun = grow_train$Tmean_AprMayJun,
#                      Tmean_SepOct = grow_train$Tmean_SepOct, DIA_prev = grow_train$DIA_prev, ELEV = grow_train$ELEV,
#                      SLOPE = grow_train$SLOPE, ASPECT = grow_train$ASPECT, tmp_yr = grow_train$tmp_yr, ppt_yr = grow_train$ppt_yr,
#                      year = grow_train$year, PlotCD = grow_train$PlotCD, treeCD = grow_train$treeCD,
#                      name = grow_train$name)
# p.o.df$overpredicted <- ifelse(p.o.df$observed <= p.o.df$ci.low, "over predicted", "within confidence interval")
# overpredictedpoints <- p.o.df %>% filter(observed <= ci.low)
# 
# 

# Pairs plots to evaluate--commented out here
# library(psych)
# xGpanels <- data.frame(xG[,1:7])
# colnames(xGpanels) <- c("ppt_norm", "tmp_norm", "Precip_JulAug", "Precip_NovDecJanFebMar", "Tmean_AprMayJun", "Tmean_SepOct", "DIA_prev")
# png("ls_pairs_panels.png", width = 8, height = 8, units = "in", res = 200)
# pairs.panels(xGpanels)
# dev.off()
#
# plotdataintervalpanels <- data.frame(plotdatainterval)
# colnames(plotdataintervalpanels) <- c("ppt_norm", "tmp_norm", "Precip_JA", "Precip_NDJFM", "Tmean_AMJ", "Tmean_SO", "size",
#                                       "pn_tn", "pn_PJA", "pn_size", "pn_PNDJFM", "pn_TAMJ", "pn_TSO", "tn_size", "tn_PJA",
#                                       "tn_PNDJFM", "tn_TAMJ", "tn_TSO", "size_PJA", "size_PNDJFM", "size_TAMJ",
#                                       "size_TSO", "PJA_PNDJFM", "PJA_TAMJ", "PJA_TSO", "PNDJFM_TAMJ", "PNDJFM_TSO", "TAMJ_TSO")
# png("ppc_pairs_panels.png", width = 15, height = 15, units = "in", res = 200)
# pairs.panels(plotdataintervalpanels)
# dev.off()
#
#


# posterior predictive checks
# ## Subset posterior predictive plot by size
# size_q<-quantile(grow$DIA_prev)
# sizeq1<-which(grow_test$DIA_prev<=size_q[2])
# sizeq2<-which(grow_test$DIA_prev>size_q[2] & grow_test$DIA_prev<=size_q[3])
# sizeq3<-which(grow_test$DIA_prev>size_q[3] &grow_test$DIA_prev<=size_q[4])
# sizeq4<-which(grow_test$DIA_prev>size_q[4])
#
# ppc_dens_overlay(yGtest[sizeq1], as.matrix(plotdata)[,sizeq1])
# ppc_dens_overlay(yGtest[sizeq2], as.matrix(plotdata)[,sizeq2])
# ppc_dens_overlay(yGtest[sizeq3], as.matrix(plotdata)[,sizeq3])
# ppc_dens_overlay(yGtest[sizeq4], as.matrix(plotdata)[,sizeq4])
#
# ## Subset posterior predicitve plot by ppt_norm
# ppt_norm_q<-quantile(grow$ppt_norm)
# ppt_normq1<-which(grow_test$ppt_norm<=ppt_norm_q[2])
# ppt_normq2<-which(grow_test$ppt_norm>ppt_norm_q[2] & grow_test$ppt_norm<=ppt_norm_q[3])
# ppt_normq3<-which(grow_test$ppt_norm>ppt_norm_q[3] &grow_test$ppt_norm<=ppt_norm_q[4])
# ppt_normq4<-which(grow_test$ppt_norm>ppt_norm_q[4])
#
# ppc_dens_overlay(yGtest[ppt_normq1], as.matrix(plotdata)[,ppt_normq1])
# ppc_dens_overlay(yGtest[ppt_normq2], as.matrix(plotdata)[,ppt_normq2])
# ppc_dens_overlay(yGtest[ppt_normq3], as.matrix(plotdata)[,ppt_normq3])
# ppc_dens_overlay(yGtest[ppt_normq4], as.matrix(plotdata)[,ppt_normq4])
#
# ## Subset posterior predicitve plot by ppt_norm
# tmp_norm_q<-quantile(grow$tmp_norm)
# tmp_normq1<-which(grow_test$tmp_norm<=tmp_norm_q[2])
# tmp_normq2<-which(grow_test$tmp_norm>tmp_norm_q[2] & grow_test$tmp_norm<=tmp_norm_q[3])
# tmp_normq3<-which(grow_test$tmp_norm>tmp_norm_q[3] &grow_test$tmp_norm<=tmp_norm_q[4])
# tmp_normq4<-which(grow_test$tmp_norm>tmp_norm_q[4])
#
# ppc_dens_overlay(yGtest[tmp_normq1], as.matrix(plotdata)[,tmp_normq1])
# ppc_dens_overlay(yGtest[tmp_normq2], as.matrix(plotdata)[,tmp_normq2])
# ppc_dens_overlay(yGtest[tmp_normq3], as.matrix(plotdata)[,tmp_normq3])
# ppc_dens_overlay(yGtest[tmp_normq4], as.matrix(plotdata)[,tmp_normq4])
#
#
# ##generating a plotting function
# ##ppcdensoverlay
# make_plot <- function() {
#   for (i in min(grow_test$year):max(grow_test$year)) {
#     year<-which(grow_test$year == i)
#     p = ppc_dens_overlay(yGtest[year], as.matrix(plotdata)[,year]) +
#       theme(
#         plot.title = element_text(size = rel(2.5), legend.text = element_text(size = 16),
#                                   axis.text.x = element_text(size = 12),
#                                   legend.key.size = unit(1.2, "lines")
#         ) + xlim(-6.91, 3.96) +
#           ggtitle(
#             paste(i)
#           ))
#     print(p)
#   }
# }
#
#
# if (!file.exists(here::here("images", "ppc_year-animation.gif"))) {
#
#   gifski::save_gif(
#     make_plot(),
#     gif_file = here::here("images", "ppc_year-animation.gif"),
#     progress = FALSE,
#     delay = 0.5,
#     height = 360, width = 640, units = "px"
#   )
# }
#
# ##density function for climate
# make_ppt_plot <- function() {
#   for (i in min(grow_test$year):max(grow_test$year)) {
#     year <- i
#     p = ggplot(grow_train[grow_train$year == year,], aes(x = ppt_yr )) + geom_density() +
#       theme(
#         plot.title = element_text(size = rel(2.5)),  legend.text = element_text(size = 16),
#         axis.text = element_text(size = 12),
#         legend.key.size = unit(1.2, "lines")
#       ) + xlim(min(grow_train$ppt_yr), max(grow_train$ppt_yr)) +
#       ggtitle(
#         paste(i)
#       )
#     print(p)
#   }
# }
#
# if (!file.exists(here::here("images", "ppt_year-animation.gif"))) {
#
#   gifski::save_gif(
#     make_ppt_plot(),
#     gif_file = here::here("images", "ppt_year-animation.gif"),
#     progress = FALSE,
#     delay = 0.5,
#     height = 360, width = 640, units = "px"
#   )
# }
#
#
# make_ppt_plot <- function() {
#   for (i in min(grow_test$year):max(grow_test$year)) {
#     year <- i
#     p = ggplot(grow_train[grow_train$year == year,], aes(x = tmp_yr )) + geom_density() +
#       theme(
#         plot.title = element_text(size = rel(2.5)),  legend.text = element_text(size = 16),
#         axis.text = element_text(size = 12),
#         legend.key.size = unit(1.2, "lines")
#       ) + xlim(min(grow_train$tmp_yr), max(grow_train$tmp_yr)) +
#       ggtitle(
#         paste(i)
#       )
#     print(p)
#   }
# }
#
# if (!file.exists(here::here("images", "tmp_year-animation.gif"))) {
#
#   gifski::save_gif(
#     make_ppt_plot(),
#     gif_file = here::here("images", "tmp_year-animation.gif"),
#     progress = FALSE,
#     delay = 0.5,
#     height = 360, width = 640, units = "px"
#   )
# }


#MCMC Intervals Plots--Keep these

# figure 2 in the Manuscript
pdf(here::here("images", "model_4", "plotdatainterval_mcmc_intervals.pdf"), height = 6, width = 5) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
mcmc_intervals(plotdatainterval, prob = 0.5, prob_outer = 0.9, point_size = 2) # put the code that makes one of the plots in here
dev.off() # "device off" tells R to stop printing stuff to the pdf

pdf(here::here("images", "model_4", "plotdatainterval_mcmc_intervals_1.3.pdf"), height = 6, width = 4) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
mcmc_intervals(plotdatainterval, prob = 0.5, prob_outer = 0.9, point_size = 2) # put the code that makes one of the plots in here
dev.off() # "device off" tells R to stop printing stuff to the pdf


#MCMC Area Plot
pdf(here::here("images", "model_4", "plotdatainterval_mcmc_areas.pdf"), height = 6, width = 10) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
color_scheme_set("purple")
mcmc_areas(plotdatainterval, prob = 0.8)
dev.off() # "device off" tells R to stop printing stuff to the pdf

#MCMC Traces

pdf(here::here("images", "model_4", "plotdatainterval_mcmc_traces.pdf"), height = 12, width = 6) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
color_scheme_set("viridis")
mcmc_trace(plotdatainterval,
           facet_args = list(nrow = 10, ncol = 3))
dev.off() # "device off" tells R to stop printing stuff to the pdf

#mcmcplot(As.mcmc.list(fit_grow))


#### Paper figures below

#Effects Plots--set a theme so we dont have to keep setting it
mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               legend.text=element_text(size=11),legend.title=element_text(size=12),
               legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
               axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
               axis.line.x = element_line(color="black", linewidth = 0.3),
               axis.line.y = element_line(color="black", linewidth = 0.3))






#### Individual Response Plots

#plotting individual tree responses to time-varying climate variables, holding tree size constant
# we make these for:
# Spring temperature, colored by MAP and with MAT bins (figure in main text)
# Fall temperature, colored by MAP and with MAT bins (S4)
# Winter Precipitation, colored by MAP and with MAT bins (S5)
# Monsoon Precipitation, colored by MAP and iwth MAT bins (S6)

#### Spring Temperature Tmean_AprMayJun
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)

grow.monsoon$tmp_norm_q <- ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.25), "0-25% quantile MAT",
                                  ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.50), "25-50% quantile MAT",
                                         ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.75), "50-75% quantile MAT",
                                                "75-100% quantile MAT")))

grow.monsoon$ppt_norm_q <- ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.25), "0-25% quantile MAP",
                                  ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.50), "25-50% quantile MAP",
                                         ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.75), "50-75% quantile MAP",
                                                "75-100% quantile MAP")))

ind.samples <- unique(grow.monsoon[,c("LATLONbin", "tmp_norm_q", "ppt_norm_q", "treeCD")])


#map of LATLONbin
all_states <- map_data("state")
states <- subset(all_states, region %in% c("arizona", "colorado", "utah", "new mexico"))
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata<-states
mapdata<-data.frame(mapdata)
ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
  geom_point(data = grow.monsoon, aes(x = LON, y = LAT, color = LATLONbin)) + theme_bw()


#Tmean_AprMayJun individual response function by ppt_norm and tmp_norm bins
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "tmp_norm_q", "ppt_norm_q", "treeCD")])

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)

  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  x <- 20#mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  #ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  # growthpredictionTmeanAprMayJun_pnorm <- matrix(NA, length(plotdatainterval$MAP), length(Tmean_AprMayJun))

  # for(i in 1:length(plotdatainterval$MAP)){
  #   growthpredictionTmeanAprMayJun_pnorm[i,] <- plotdatainterval[i,"MAP"]*ppt_norm +
  #     plotdatainterval[i,"MAT"]*tmp_norm +
  #     plotdatainterval[i, "MAP*MAT"]*tmp_norm*ppt_norm +
  #     plotdatainterval[i,"monsoon precip"]*Precip_JulAug +
  #     plotdatainterval[i,"winter precip"]*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"spring temp"]*Tmean_AprMayJun +
  #     plotdatainterval[i,"fall temp"]*Tmean_SepOct +
  #     plotdatainterval[i,"tree size"]*x +
  #     plotdatainterval[i,"MAP*monsoon precip"]*ppt_norm*Precip_JulAug +
  #     plotdatainterval[i,"tree size*MAP"]*ppt_norm*x +
  #     plotdatainterval[i,"MAP*winter precip"]*ppt_norm*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"MAP*spring temp"]*ppt_norm*Tmean_AprMayJun +
  #     plotdatainterval[i,"MAP*fall temp"]*ppt_norm*Tmean_SepOct +
  #     plotdatainterval[i,"tree size*MAT"]*tmp_norm*x +
  #     plotdatainterval[i,"MAT*monsoon precip"]*tmp_norm*Precip_JulAug +
  #     plotdatainterval[i,"MAT*winter precip"]*tmp_norm*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"MAT*spring temp"]*tmp_norm*Tmean_AprMayJun +
  #     plotdatainterval[i,"MAT*fall temp"]*tmp_norm*Tmean_SepOct +
  #     plotdatainterval[i,"tree size*Precip_JulAug"]*x*Precip_JulAug +
  #     plotdatainterval[i,"tree size*winter precip"]*x*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"tree size*spring temp"]*x*Tmean_AprMayJun +
  #     plotdatainterval[i,"tree size*fall temp"]*x*Tmean_SepOct +
  #     plotdatainterval[i,"monsoon precip*winter precip"]*Precip_JulAug*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"monsoon precip*spring temp"]*Precip_JulAug*Tmean_AprMayJun +
  #     plotdatainterval[i,"monsoon precip*fall temp"]*Precip_JulAug*Tmean_SepOct +
  #     plotdatainterval[i,"winter precip*spring temp"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
  #     plotdatainterval[i,"winter precip*fall temp"]*Precip_NovDecJanFebMar*Tmean_SepOct +
  #     plotdatainterval[i,"spring temp*fall temp"]*Tmean_AprMayJun*Tmean_SepOct
  # }
  n_cols <- length(Tmean_AprMayJun)
  n_rows <- length(plotdatainterval$MAP)
  growthpredictionTmeanAprMayJun_pnorm <- matrix(plotdatainterval[["MAP"]], n_rows, n_cols) * ppt_norm +
    
    matrix(plotdatainterval[["MAT"]], n_rows, n_cols) * tmp_norm +
    matrix(plotdatainterval[["MAP*MAT"]], n_rows, n_cols) * tmp_norm*ppt_norm +
    matrix(plotdatainterval[["monsoon precip"]], n_rows, n_cols) * Precip_JulAug +
    matrix(plotdatainterval[["winter precip"]], n_rows, n_cols) * Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["spring temp"]], n_rows, n_cols) * matrix(Tmean_AprMayJun, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["fall temp"]], n_rows, n_cols) * Tmean_SepOct +
    matrix(plotdatainterval[["tree size"]], n_rows, n_cols) * x +
    matrix(plotdatainterval[["MAP*monsoon precip"]], n_rows, n_cols) * ppt_norm*Precip_JulAug +
    matrix(plotdatainterval[["tree size*MAP"]], n_rows, n_cols) * ppt_norm*x +
    matrix(plotdatainterval[["MAP*winter precip"]], n_rows, n_cols) * ppt_norm*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["MAP*spring temp"]], n_rows, n_cols) * ppt_norm*matrix(Tmean_AprMayJun, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["MAP*fall temp"]], n_rows, n_cols) * ppt_norm*Tmean_SepOct +
    matrix(plotdatainterval[["tree size*MAT"]], n_rows, n_cols) * tmp_norm*x +
    matrix(plotdatainterval[["MAT*monsoon precip"]], n_rows, n_cols) * tmp_norm*Precip_JulAug +
    matrix(plotdatainterval[["MAT*winter precip"]], n_rows, n_cols) * tmp_norm*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["MAT*spring temp"]], n_rows, n_cols) * tmp_norm*matrix(Tmean_AprMayJun, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["MAT*fall temp"]], n_rows, n_cols) * tmp_norm*Tmean_SepOct +
    matrix(plotdatainterval[["tree size*monsoon precip"]], n_rows, n_cols) * x*Precip_JulAug +
    matrix(plotdatainterval[["tree size*winter precip"]], n_rows, n_cols) * x*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["tree size*spring temp"]], n_rows, n_cols) * x*matrix(Tmean_AprMayJun, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["tree size*fall temp"]], n_rows, n_cols) * x*Tmean_SepOct +
    matrix(plotdatainterval[["monsoon precip*winter precip"]], n_rows, n_cols) * Precip_JulAug*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["monsoon precip*spring temp"]], n_rows, n_cols) * Precip_JulAug*matrix(Tmean_AprMayJun, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["monsoon precip*fall temp"]], n_rows, n_cols) * Precip_JulAug*Tmean_SepOct +
    matrix(plotdatainterval[["winter precip*spring temp"]], n_rows, n_cols) * Precip_NovDecJanFebMar*matrix(Tmean_AprMayJun, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["winter precip*fall temp"]], n_rows, n_cols) * Precip_NovDecJanFebMar*Tmean_SepOct +
    matrix(plotdatainterval[["spring temp*fall temp"]], n_rows, n_cols) * matrix(Tmean_AprMayJun, n_rows, n_cols, byrow = TRUE)*Tmean_SepOct
  Tmean_AprMayJun_prediction_trpnorm <- exp(growthpredictionTmeanAprMayJun_pnorm)
  ci.Tmean_AprMayJunpnorm <- apply(Tmean_AprMayJun_prediction_trpnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJunpnorm.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, ppt_norm = ppt_norm, median = ci.Tmean_AprMayJunpnorm[2,], ci.low = ci.Tmean_AprMayJunpnorm[1,], ci.high = ci.Tmean_AprMayJunpnorm[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_pnormint <- rbind(ci.Tmean_AprMayJunpnorm.df)
  print(ind.samples[j,])
  return(Tmean_AprMayJun_pnormint)
}

#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
# #color by group
# ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# #color by LATLONbin
# ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)

colnames(merged.response.samples)[3] <- "MAP"

#color by ppt_norm & use tmp_norm bins
png(here::here("images", "model_4", "individual_response_Tmean_APRMAYJUN_MAP.png"), height = 5, width = 6, units = "in", res = 300) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = MAP, group = ci.group)) + 
  geom_line(alpha = 0.5) +
  mytheme + 
  ylab("Predicted Growth") + 
  scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4") + 
  facet_wrap(~tmp_norm_q) +
  ylab("Predicted Growth (mm)") + 
  xlab("Spring Temperature Anomaly" ) + 
  ylim(0, 3)
dev.off()




#### monsoon precip
#Precip_JulAug and tmp_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)

grow.monsoon$tmp_norm_q <- ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.25), "0-25% quantile MAT",
                                  ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.50), "25-50% quantile MAT",
                                         ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.75), "50-75% quantile MAT",
                                                "75-100% quantile MAT")))

grow.monsoon$ppt_norm_q <- ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.25), "0-25% quantile MAP",
                                  ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.50), "25-50% quantile MAP",
                                         ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.75), "50-75% quantile MAP",
                                                "75-100% quantile MAP")))

ind.samples <- unique(grow.monsoon[,c("LATLONbin", "tmp_norm_q", "ppt_norm_q", "treeCD", "ppt_norm")])


get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)

  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  x <- 20
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  # growthpredictionPrecipJulAug_tmpnorm <- matrix(NA, length(plotdatainterval$`monsoon precip`), length(Precip_JulAug))
  n_cols <- length(Precip_JulAug)
  n_rows <- length(plotdatainterval$MAP)
  growthpredictionPrecipJulAug_tmpnorm <- matrix(plotdatainterval[["MAT"]], n_rows, n_cols) * tmp_norm +
    matrix(plotdatainterval[["MAP*MAT"]], n_rows, n_cols) * tmp_norm*ppt_norm +
    matrix(plotdatainterval[["monsoon precip"]], n_rows, n_cols) * matrix(Precip_JulAug, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["winter precip"]], n_rows, n_cols) * Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["spring temp"]], n_rows, n_cols) * Tmean_AprMayJun +
    matrix(plotdatainterval[["fall temp"]], n_rows, n_cols) * Tmean_SepOct +
    matrix(plotdatainterval[["tree size"]], n_rows, n_cols) * x +
    matrix(plotdatainterval[["MAP*monsoon precip"]], n_rows, n_cols) * ppt_norm*matrix(Precip_JulAug, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["tree size*MAP"]], n_rows, n_cols) * ppt_norm*x +
    matrix(plotdatainterval[["MAP*winter precip"]], n_rows, n_cols) * ppt_norm*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["MAP*spring temp"]], n_rows, n_cols) * ppt_norm*Tmean_AprMayJun +
    matrix(plotdatainterval[["MAP*fall temp"]], n_rows, n_cols) * ppt_norm*Tmean_SepOct +
    matrix(plotdatainterval[["tree size*MAT"]], n_rows, n_cols) * tmp_norm*x +
    matrix(plotdatainterval[["MAT*monsoon precip"]], n_rows, n_cols) * tmp_norm*matrix(Precip_JulAug, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["MAT*winter precip"]], n_rows, n_cols) * tmp_norm*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["MAT*spring temp"]], n_rows, n_cols) * tmp_norm*Tmean_AprMayJun +
    matrix(plotdatainterval[["MAT*fall temp"]], n_rows, n_cols) * tmp_norm*Tmean_SepOct +
    matrix(plotdatainterval[["tree size*monsoon precip"]], n_rows, n_cols) * x*matrix(Precip_JulAug, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["tree size*winter precip"]], n_rows, n_cols) * x*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["tree size*spring temp"]], n_rows, n_cols) * x*Tmean_AprMayJun +
    matrix(plotdatainterval[["tree size*fall temp"]], n_rows, n_cols) * x*Tmean_SepOct +
    matrix(plotdatainterval[["monsoon precip*winter precip"]], n_rows, n_cols) * matrix(Precip_JulAug, n_rows, n_cols, byrow = TRUE)*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["monsoon precip*spring temp"]], n_rows, n_cols) * matrix(Precip_JulAug, n_rows, n_cols, byrow = TRUE)*Tmean_AprMayJun +
    matrix(plotdatainterval[["monsoon precip*fall temp"]], n_rows, n_cols) * matrix(Precip_JulAug, n_rows, n_cols, byrow = TRUE)*Tmean_SepOct +
    matrix(plotdatainterval[["winter precip*spring temp"]], n_rows, n_cols) * Precip_NovDecJanFebMar*Tmean_AprMayJun +
    matrix(plotdatainterval[["winter precip*fall temp"]], n_rows, n_cols) * Precip_NovDecJanFebMar*Tmean_SepOct +
    matrix(plotdatainterval[["spring temp*fall temp"]], n_rows, n_cols) * Tmean_AprMayJun*Tmean_SepOct
  # for(i in 1:length(plotdatainterval$MAP)){
  #   growthpredictionPrecipJulAug_tmpnorm[i,] <- plotdatainterval[i,"MAP"]*ppt_norm +
  #     plotdatainterval[i,"MAT"]*tmp_norm +
  #     plotdatainterval[i, "MAP*MAT"]*tmp_norm*ppt_norm +
  #     plotdatainterval[i,"monsoon precip"]*Precip_JulAug +
  #     plotdatainterval[i,"winter precip"]*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"spring temp"]*Tmean_AprMayJun +
  #     plotdatainterval[i,"fall temp"]*Tmean_SepOct +
  #     plotdatainterval[i,"tree size"]*x +
  #     plotdatainterval[i,"MAP*monsoon precip"]*ppt_norm*Precip_JulAug +
  #     plotdatainterval[i,"tree size*MAP"]*ppt_norm*x +
  #     plotdatainterval[i,"MAP*winter precip"]*ppt_norm*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"MAP*spring temp"]*ppt_norm*Tmean_AprMayJun +
  #     plotdatainterval[i,"MAP*fall temp"]*ppt_norm*Tmean_SepOct +
  #     plotdatainterval[i,"tree size*MAT"]*tmp_norm*x +
  #     plotdatainterval[i,"MAT*monsoon precip"]*tmp_norm*Precip_JulAug +
  #     plotdatainterval[i,"MAT*winter precip"]*tmp_norm*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"MAT*spring temp"]*tmp_norm*Tmean_AprMayJun +
  #     plotdatainterval[i,"MAT*fall temp"]*tmp_norm*Tmean_SepOct +
  #     plotdatainterval[i,"tree size*Precip_JulAug"]*x*Precip_JulAug +
  #     plotdatainterval[i,"tree size*winter precip"]*x*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"tree size*spring temp"]*x*Tmean_AprMayJun +
  #     plotdatainterval[i,"tree size*fall temp"]*x*Tmean_SepOct +
  #     plotdatainterval[i,"monsoon precip*winter precip"]*Precip_JulAug*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"monsoon precip*spring temp"]*Precip_JulAug*Tmean_AprMayJun +
  #     plotdatainterval[i,"monsoon precip*fall temp"]*Precip_JulAug*Tmean_SepOct +
  #     plotdatainterval[i,"winter precip*spring temp"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
  #     plotdatainterval[i,"winter precip*fall temp"]*Precip_NovDecJanFebMar*Tmean_SepOct +
  #     plotdatainterval[i,"spring temp*fall temp"]*Tmean_AprMayJun*Tmean_SepOct
  # }
  Precip_JulAug_prediction_trtmpnorm <- exp(growthpredictionPrecipJulAug_tmpnorm)
  ci.Precip_JulAugtmpnorm <- apply(Precip_JulAug_prediction_trtmpnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_JulAugtmpnorm.df <- data.frame(Precip_JulAug = Precip_JulAug, tmp_norm = tmp_norm, median = ci.Precip_JulAugtmpnorm[2,], ci.low = ci.Precip_JulAugtmpnorm[1,], ci.high = ci.Precip_JulAugtmpnorm[3,], ci.group = tree.subset$treeCD)
  Precip_JulAug_tmpnormint <- rbind(ci.Precip_JulAugtmpnorm.df)
  print(ind.samples[j,])
  Precip_JulAug_tmpnormint
}
#get.ind.tmp.response(i = 6)
Precip_JulAug_tree_response <- list()
Precip_JulAug_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_JulAug_tree_response.df <- do.call(rbind, Precip_JulAug_tree_response)
merged.response.samples <- merge(Precip_JulAug_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
colnames(merged.response.samples)[3] <- "MAT"
colnames(merged.response.samples)[10] <- "MAP"

#color by ppt_norm & use tmp_norm bins
png(here::here("images", "model_4", "individual_response_Precip_JulAug_MAP.png"), height = 5, width = 6, units = "in", res = 300) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = MAP, group = ci.group)) + geom_line(alpha = 0.5) +
  mytheme + ylab("Predicted Growth") + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")+ facet_wrap(~tmp_norm_q)+
  ylab("Predicted Growth (mm)")+xlab("Scaled Monsoon Precipitation")+ylim(0, 3)
dev.off()


#### winter precip:
#Precip_NovDecJanFebMar and tmp_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)

grow.monsoon$tmp_norm_q <- ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.25), "0-25% quantile MAT",
                                  ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.50), "25-50% quantile MAT",
                                         ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.75), "50-75% quantile MAT",
                                                "75-100% quantile MAT")))

grow.monsoon$ppt_norm_q <- ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.25), "0-25% quantile MAP",
                                  ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.50), "25-50% quantile MAP",
                                         ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.75), "50-75% quantile MAP",
                                                "75-100% quantile MAP")))

ind.samples <- unique(grow.monsoon[,c("LATLONbin", "tmp_norm_q", "ppt_norm_q", "treeCD", "ppt_norm")])


get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)

  Precip_NovDecJanFebMarrng <- range(tree.grow$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
  Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.1)
  x <- 20
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  n_cols <- length(Precip_NovDecJanFebMar)
  n_rows <- length(plotdatainterval$MAP)
  growthpredictionPrecipNovDecJanFebMar_tmpnorm <- matrix(plotdatainterval[["MAT"]], n_rows, n_cols) * tmp_norm +
    matrix(plotdatainterval[["MAP*MAT"]], n_rows, n_cols) * tmp_norm*ppt_norm +
    matrix(plotdatainterval[["monsoon precip"]], n_rows, n_cols) * Precip_JulAug +
    matrix(plotdatainterval[["winter precip"]], n_rows, n_cols) * matrix(Precip_NovDecJanFebMar, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["spring temp"]], n_rows, n_cols) * Tmean_AprMayJun +
    matrix(plotdatainterval[["fall temp"]], n_rows, n_cols) * Tmean_SepOct +
    matrix(plotdatainterval[["tree size"]], n_rows, n_cols) * x +
    matrix(plotdatainterval[["MAP*monsoon precip"]], n_rows, n_cols) * ppt_norm*Precip_JulAug +
    matrix(plotdatainterval[["tree size*MAP"]], n_rows, n_cols) * ppt_norm*x +
    matrix(plotdatainterval[["MAP*winter precip"]], n_rows, n_cols) * ppt_norm*matrix(Precip_NovDecJanFebMar, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["MAP*spring temp"]], n_rows, n_cols) * ppt_norm*Tmean_AprMayJun +
    matrix(plotdatainterval[["MAP*fall temp"]], n_rows, n_cols) * ppt_norm*Tmean_SepOct +
    matrix(plotdatainterval[["tree size*MAT"]], n_rows, n_cols) * tmp_norm*x +
    matrix(plotdatainterval[["MAT*monsoon precip"]], n_rows, n_cols) * tmp_norm*Precip_JulAug +
    matrix(plotdatainterval[["MAT*winter precip"]], n_rows, n_cols) * tmp_norm*matrix(Precip_NovDecJanFebMar, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["MAT*spring temp"]], n_rows, n_cols) * tmp_norm*Tmean_AprMayJun +
    matrix(plotdatainterval[["MAT*fall temp"]], n_rows, n_cols) * tmp_norm*Tmean_SepOct +
    matrix(plotdatainterval[["tree size*monsoon precip"]], n_rows, n_cols) * x*Precip_JulAug +
    matrix(plotdatainterval[["tree size*winter precip"]], n_rows, n_cols) * x*matrix(Precip_NovDecJanFebMar, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["tree size*spring temp"]], n_rows, n_cols) * x*Tmean_AprMayJun +
    matrix(plotdatainterval[["tree size*fall temp"]], n_rows, n_cols) * x*Tmean_SepOct +
    matrix(plotdatainterval[["monsoon precip*winter precip"]], n_rows, n_cols) * Precip_JulAug*matrix(Precip_NovDecJanFebMar, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["monsoon precip*spring temp"]], n_rows, n_cols) * Precip_JulAug*Tmean_AprMayJun +
    matrix(plotdatainterval[["monsoon precip*fall temp"]], n_rows, n_cols) * Precip_JulAug*Tmean_SepOct +
    matrix(plotdatainterval[["winter precip*spring temp"]], n_rows, n_cols) * matrix(Precip_NovDecJanFebMar, n_rows, n_cols, byrow = TRUE)*Tmean_AprMayJun +
    matrix(plotdatainterval[["winter precip*fall temp"]], n_rows, n_cols) * matrix(Precip_NovDecJanFebMar, n_rows, n_cols, byrow = TRUE)*Tmean_SepOct +
    matrix(plotdatainterval[["spring temp*fall temp"]], n_rows, n_cols) * Tmean_AprMayJun*Tmean_SepOct
  
  # growthpredictionPrecipNovDecJanFebMar_tmpnorm <- matrix(NA, length(plotdatainterval$`winter precip`), length(Precip_NovDecJanFebMar))
  # 
  # for(i in 1:length(plotdatainterval$MAP)){
  #   growthpredictionPrecipNovDecJanFebMar_tmpnorm [i,] <- plotdatainterval[i,"MAP"]*ppt_norm +
  #     plotdatainterval[i,"MAT"]*tmp_norm +
  #     plotdatainterval[i, "MAP*MAT"]*tmp_norm*ppt_norm +
  #     plotdatainterval[i,"monsoon precip"]*Precip_JulAug +
  #     plotdatainterval[i,"winter precip"]*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"spring temp"]*Tmean_AprMayJun +
  #     plotdatainterval[i,"fall temp"]*Tmean_SepOct +
  #     plotdatainterval[i,"tree size"]*x +
  #     plotdatainterval[i,"MAP*monsoon precip"]*ppt_norm*Precip_JulAug +
  #     plotdatainterval[i,"tree size*MAP"]*ppt_norm*x +
  #     plotdatainterval[i,"MAP*winter precip"]*ppt_norm*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"MAP*spring temp"]*ppt_norm*Tmean_AprMayJun +
  #     plotdatainterval[i,"MAP*fall temp"]*ppt_norm*Tmean_SepOct +
  #     plotdatainterval[i,"tree size*MAT"]*tmp_norm*x +
  #     plotdatainterval[i,"MAT*monsoon precip"]*tmp_norm*Precip_JulAug +
  #     plotdatainterval[i,"MAT*winter precip"]*tmp_norm*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"MAT*spring temp"]*tmp_norm*Tmean_AprMayJun +
  #     plotdatainterval[i,"MAT*fall temp"]*tmp_norm*Tmean_SepOct +
  #     plotdatainterval[i,"tree size*Precip_JulAug"]*x*Precip_JulAug +
  #     plotdatainterval[i,"tree size*winter precip"]*x*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"tree size*spring temp"]*x*Tmean_AprMayJun +
  #     plotdatainterval[i,"tree size*fall temp"]*x*Tmean_SepOct +
  #     plotdatainterval[i,"monsoon precip*winter precip"]*Precip_JulAug*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"monsoon precip*spring temp"]*Precip_JulAug*Tmean_AprMayJun +
  #     plotdatainterval[i,"monsoon precip*fall temp"]*Precip_JulAug*Tmean_SepOct +
  #     plotdatainterval[i,"winter precip*spring temp"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
  #     plotdatainterval[i,"winter precip*fall temp"]*Precip_NovDecJanFebMar*Tmean_SepOct +
  #     plotdatainterval[i,"spring temp*fall temp"]*Tmean_AprMayJun*Tmean_SepOct
  # }
  Precip_NovDecJanFebMar_prediction_trtmpnorm <- exp(growthpredictionPrecipNovDecJanFebMar_tmpnorm)
  ci.Precip_NovDecJanFebMartmpnorm <- apply(Precip_NovDecJanFebMar_prediction_trtmpnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_NovDecJanFebMartmpnorm.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, tmp_norm = tmp_norm, median = ci.Precip_NovDecJanFebMartmpnorm[2,], ci.low = ci.Precip_NovDecJanFebMartmpnorm[1,], ci.high = ci.Precip_NovDecJanFebMartmpnorm[3,], ci.group = tree.subset$treeCD)
  Precip_NovDecJanFebMar_tmpnormint <- rbind(ci.Precip_NovDecJanFebMartmpnorm.df)
  print(ind.samples[j,])
  Precip_NovDecJanFebMar_tmpnormint
}
#get.ind.tmp.response(i = 6)
Precip_NovDecJanFebMar_tree_response <- list()
Precip_NovDecJanFebMar_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_NovDecJanFebMar_tree_response.df <- do.call(rbind, Precip_NovDecJanFebMar_tree_response)
merged.response.samples <- merge(Precip_NovDecJanFebMar_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)

colnames(merged.response.samples)[10] <- "MAP"
#color by ppt_norm & use tmp_norm bins
png(here::here("images", "model_4", "individual_response_Precip_NovDecJanFebMar_MAP.png"), height = 5, width = 6, units = "in", res = 300) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = MAP, group = ci.group)) + geom_line(alpha = 0.5) +
  mytheme + ylab("Predicted Growth") + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")+ facet_wrap(~tmp_norm_q)+
  ylab("Predicted Growth (mm)")+xlab("Scaled Winter Precipitation")+ylim(0, 3)
dev.off()


#### Fall temperature
#Tmean_SepOct and tmp_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

grow.monsoon$tmp_norm_q <- ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.25), "0-25% quantile MAT",
                                  ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.50), "25-50% quantile MAT",
                                         ifelse(grow.monsoon$tmp_norm <= quantile(grow.monsoon$tmp_norm, 0.75), "50-75% quantile MAT",
                                                "75-100% quantile MAT")))

grow.monsoon$ppt_norm_q <- ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.25), "0-25% quantile MAP",
                                  ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.50), "25-50% quantile MAP",
                                         ifelse(grow.monsoon$ppt_norm <= quantile(grow.monsoon$ppt_norm, 0.75), "50-75% quantile MAP",
                                                "75-100% quantile MAP")))

ind.samples <- unique(grow.monsoon[,c("LATLONbin", "tmp_norm_q", "ppt_norm_q", "treeCD", "ppt_norm")])


get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)

  Tmean_SepOctrng <- range(tree.grow$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.1)
  x <- 20
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  n_cols <- length(Precip_NovDecJanFebMar)
  n_rows <- length(plotdatainterval$MAP)
  growthpredictionTmeanSepOct_tmpnorm <- matrix(plotdatainterval[["MAT"]], n_rows, n_cols) * tmp_norm +
    matrix(plotdatainterval[["MAP*MAT"]], n_rows, n_cols) * tmp_norm*ppt_norm +
    matrix(plotdatainterval[["monsoon precip"]], n_rows, n_cols) * Precip_JulAug +
    matrix(plotdatainterval[["winter precip"]], n_rows, n_cols) * Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["spring temp"]], n_rows, n_cols) * Tmean_AprMayJun +
    matrix(plotdatainterval[["fall temp"]], n_rows, n_cols) * matrix(Tmean_SepOct, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["tree size"]], n_rows, n_cols) * x +
    matrix(plotdatainterval[["MAP*monsoon precip"]], n_rows, n_cols) * ppt_norm*Precip_JulAug +
    matrix(plotdatainterval[["tree size*MAP"]], n_rows, n_cols) * ppt_norm*x +
    matrix(plotdatainterval[["MAP*winter precip"]], n_rows, n_cols) * ppt_norm*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["MAP*spring temp"]], n_rows, n_cols) * ppt_norm*Tmean_AprMayJun +
    matrix(plotdatainterval[["MAP*fall temp"]], n_rows, n_cols) * ppt_norm*matrix(Tmean_SepOct, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["tree size*MAT"]], n_rows, n_cols) * tmp_norm*x +
    matrix(plotdatainterval[["MAT*monsoon precip"]], n_rows, n_cols) * tmp_norm*Precip_JulAug +
    matrix(plotdatainterval[["MAT*winter precip"]], n_rows, n_cols) * tmp_norm*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["MAT*spring temp"]], n_rows, n_cols) * tmp_norm*Tmean_AprMayJun +
    matrix(plotdatainterval[["MAT*fall temp"]], n_rows, n_cols) * tmp_norm*matrix(Tmean_SepOct, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["tree size*monsoon precip"]], n_rows, n_cols) * x*Precip_JulAug +
    matrix(plotdatainterval[["tree size*winter precip"]], n_rows, n_cols) * x*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["tree size*spring temp"]], n_rows, n_cols) * x*Tmean_AprMayJun +
    matrix(plotdatainterval[["tree size*fall temp"]], n_rows, n_cols) * x*matrix(Tmean_SepOct, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["monsoon precip*winter precip"]], n_rows, n_cols) * Precip_JulAug*Precip_NovDecJanFebMar +
    matrix(plotdatainterval[["monsoon precip*spring temp"]], n_rows, n_cols) * Precip_JulAug*Tmean_AprMayJun +
    matrix(plotdatainterval[["monsoon precip*fall temp"]], n_rows, n_cols) * Precip_JulAug*matrix(Tmean_SepOct, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["winter precip*spring temp"]], n_rows, n_cols) * Precip_NovDecJanFebMar*Tmean_AprMayJun +
    matrix(plotdatainterval[["winter precip*fall temp"]], n_rows, n_cols) * Precip_NovDecJanFebMar*matrix(Tmean_SepOct, n_rows, n_cols, byrow = TRUE) +
    matrix(plotdatainterval[["spring temp*fall temp"]], n_rows, n_cols) * Tmean_AprMayJun*matrix(Tmean_SepOct, n_rows, n_cols, byrow = TRUE)
  # growthpredictionTmeanSepOct_tmpnorm <- matrix(NA, length(plotdatainterval$`fall temp`), length(Tmean_SepOct))
  # 
  # for(i in 1:length(plotdatainterval$MAP)){
  #   growthpredictionTmeanSepOct_tmpnorm [i,] <- plotdatainterval[i,"MAP"]*ppt_norm +
  #     plotdatainterval[i,"MAT"]*tmp_norm +
  #     plotdatainterval[i, "MAP*MAT"]*tmp_norm*ppt_norm +
  #     plotdatainterval[i,"monsoon precip"]*Precip_JulAug +
  #     plotdatainterval[i,"winter precip"]*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"spring temp"]*Tmean_AprMayJun +
  #     plotdatainterval[i,"fall temp"]*Tmean_SepOct +
  #     plotdatainterval[i,"tree size"]*x +
  #     plotdatainterval[i,"MAP*monsoon precip"]*ppt_norm*Precip_JulAug +
  #     plotdatainterval[i,"tree size*MAP"]*ppt_norm*x +
  #     plotdatainterval[i,"MAP*winter precip"]*ppt_norm*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"MAP*spring temp"]*ppt_norm*Tmean_AprMayJun +
  #     plotdatainterval[i,"MAP*fall temp"]*ppt_norm*Tmean_SepOct +
  #     plotdatainterval[i,"tree size*MAT"]*tmp_norm*x +
  #     plotdatainterval[i,"MAT*monsoon precip"]*tmp_norm*Precip_JulAug +
  #     plotdatainterval[i,"MAT*winter precip"]*tmp_norm*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"MAT*spring temp"]*tmp_norm*Tmean_AprMayJun +
  #     plotdatainterval[i,"MAT*fall temp"]*tmp_norm*Tmean_SepOct +
  #     plotdatainterval[i,"tree size*Precip_JulAug"]*x*Precip_JulAug +
  #     plotdatainterval[i,"tree size*winter precip"]*x*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"tree size*spring temp"]*x*Tmean_AprMayJun +
  #     plotdatainterval[i,"tree size*fall temp"]*x*Tmean_SepOct +
  #     plotdatainterval[i,"monsoon precip*winter precip"]*Precip_JulAug*Precip_NovDecJanFebMar +
  #     plotdatainterval[i,"monsoon precip*spring temp"]*Precip_JulAug*Tmean_AprMayJun +
  #     plotdatainterval[i,"monsoon precip*fall temp"]*Precip_JulAug*Tmean_SepOct +
  #     plotdatainterval[i,"winter precip*spring temp"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
  #     plotdatainterval[i,"winter precip*fall temp"]*Precip_NovDecJanFebMar*Tmean_SepOct +
  #     plotdatainterval[i,"spring temp*fall temp"]*Tmean_AprMayJun*Tmean_SepOct
  # }
  Tmean_SepOct_prediction_trtmpnorm <- exp(growthpredictionTmeanSepOct_tmpnorm)
  ci.Tmean_SepOcttmpnorm <- apply(Tmean_SepOct_prediction_trtmpnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_SepOcttmpnorm.df <- data.frame(Tmean_SepOct = Tmean_SepOct, tmp_norm = tmp_norm, median = ci.Tmean_SepOcttmpnorm[2,], ci.low = ci.Tmean_SepOcttmpnorm[1,], ci.high = ci.Tmean_SepOcttmpnorm[3,], ci.group = tree.subset$treeCD)
  Tmean_SepOct_tmpnormint <- rbind(ci.Tmean_SepOcttmpnorm.df)
  print(ind.samples[j,])
  Tmean_SepOct_tmpnormint
}
#get.ind.tmp.response(i = 6)
Tmean_SepOct_tree_response <- list()
Tmean_SepOct_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_SepOct_tree_response.df <- do.call(rbind, Tmean_SepOct_tree_response)
merged.response.samples <- merge(Tmean_SepOct_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)

colnames(merged.response.samples)[10]<- "MAP"
#color by ppt_norm & use tmp_norm bins
png(here::here("images", "model_4", "individual_response_Tmean_SepOct_MAP.png"), height = 5, width = 6, units = "in", res = 300) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = MAP, group = ci.group)) + geom_line(alpha = 0.5) +
  mytheme + ylab("Predicted Growth") + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")+ facet_wrap(~tmp_norm_q)+
  ylab("Predicted Growth (mm)")+xlab("Scaled Fall Temperature")+ylim(0, 3)
dev.off()


#####  Extra plots ####
# mapping out MAP and MAT in space
#ggplot map of tmp_norm
all_states <- map_data("state")
states <- subset(all_states, region %in% c("arizona", "colorado", "utah", "new mexico"))
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata<-states
mapdata<-data.frame(mapdata)
ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
  geom_point(data = grow.monsoon, aes(x = LON, y = LAT, color = tmp_norm))+ scale_color_gradient2(low = "#4575b4", mid = "#fddbc7", high = "#b2182b")+
  theme_bw()


#ggplot map of ppt_norm

ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
  geom_point(data = grow.monsoon, aes(x = LON, y = LAT, color = ppt_norm)) + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4") +
  theme_bw()


#ggplots of monsoon climate
monsoonxtmpyr <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = tmp_yr)) + geom_point()
monsoonxtmpyr

monsoonxtmpnorm <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = tmp_norm)) + geom_point()
monsoonxtmpnorm

monsoonxpptnorm <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = ppt_norm)) + geom_point()
monsoonxpptnorm

monsoonxsize <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = DIA_prev)) + geom_point()
monsoonxsize

coolxmonsoon <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = Precip_DecJanFeb)) + geom_point()
coolxmonsoon

pfallxmonsoon <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = Precip_SepOctNov)) + geom_point()
pfallxmonsoon

fsummerxmonsoon <-ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = Precip_MarAprMay)) + geom_point()
fsummerxmonsoon







#### the following could be used to make a whole bunch of basic interaction/population response plots:
## The main thing that needs to be changed on all of these is the names of the parameters for model 4
# #effect of ppt_norm
# ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for ppt_normrng
# ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.25)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# growthpredictionpptnorm <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm))
#
#
# pfun_plotdatainterval <- function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
#   for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
#     growthpredictionpptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#       plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#   }
#   growthpredictionpptnorm
# }
# ppt_norm_prediction <- pfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug,
#                                              Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
# ppt_norm_prediction_tr <- exp(ppt_norm_prediction)
# ci.ppt_norm <- apply(ppt_norm_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.ppt_norm.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_norm[2,], ci.low = ci.ppt_norm[1,], ci.high = ci.ppt_norm[3,])
# ggplot() +
#   geom_ribbon(data = ci.ppt_norm.df, aes(x = ppt_norm, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") +
#   geom_line(data = ci.ppt_norm.df, aes(x = ppt_norm, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) +
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "tmp_norm")]), aes(x = ppt_norm))
#
#
# #effect of tmp_norm
# tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
# tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.35)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# growthpredictiontmpnorm <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm))
# tfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
#   for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
#     growthpredictiontmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#       plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#   }
#   growthpredictiontmpnorm
# }
# tmp_norm_prediction <- tfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug,
#                                              Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
# tmp_norm_prediction_tr <- exp(tmp_norm_prediction)
# ci.tmp_norm <- apply(tmp_norm_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_norm.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_norm[2,], ci.low = ci.tmp_norm[1,], ci.high = ci.tmp_norm[3,])
# ggplot() +
#   geom_ribbon(data = ci.tmp_norm.df, aes(x = tmp_norm, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") +
#   geom_line(data = ci.tmp_norm.df, aes(x = tmp_norm, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) +
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "tmp_norm")]), aes(x = tmp_norm))
#
#
# #effect of Precip_JulAug
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# growthpredictionPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
# pjafun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
#   for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#     growthpredictionPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#       plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#   }
#   growthpredictionPrecipJulAug
# }
# PrecipJulAug_prediction <- pjafun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug,
#                                                    Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
# PrecipJulAug_prediction_tr <- exp(PrecipJulAug_prediction)
# ci.Precip_JulAug <- apply(PrecipJulAug_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAug.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAug[2,], ci.low = ci.Precip_JulAug[1,], ci.high = ci.Precip_JulAug[3,])
# ggplot() +
#   geom_ribbon(data = ci.Precip_JulAug.df, aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") +
#   geom_line(data = ci.Precip_JulAug.df, aes(x = Precip_JulAug, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) +
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "Precip_JulAug")]), aes(x = Precip_JulAug))
#
#
# #effect of Precip_NovDecJanFebMar
# Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
# Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# growthpredictionPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar))
# pndjfmfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
#   for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
#     growthpredictionPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#       plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#   }
#   growthpredictionPrecipNovDecJanFebMar
# }
# PrecipNovDecJanFebMar_prediction <- pndjfmfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug,
#                                                                Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
# PrecipNovDecJanFebMar_prediction_tr <- exp(PrecipNovDecJanFebMar_prediction)
# ci.Precip_NovDecJanFebMar <- apply(PrecipNovDecJanFebMar_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_NovDecJanFebMar.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_NovDecJanFebMar[3,])
# ggplot() +
#   geom_ribbon(data = ci.Precip_NovDecJanFebMar.df, aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") +
#   geom_line(data = ci.Precip_NovDecJanFebMar.df, aes(x = Precip_NovDecJanFebMar, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) +
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "Precip_NovDecJanFebMar")]), aes(x = Precip_NovDecJanFebMar))
#
# #effect of Tmean_AprMayJun
# Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# growthpredictionTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
# tamjfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
#   for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
#     growthpredictionTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#       plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#   }
#   growthpredictionTmeanAprMayJun
# }
# TmeanAprMayJun_prediction <- tamjfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug,
#                                                       Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
# TmeanAprMayJun_prediction_tr <- exp(TmeanAprMayJun_prediction)
# ci.Tmean_AprMayJun <- apply(TmeanAprMayJun_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_AprMayJun.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJun[2,], ci.low = ci.Tmean_AprMayJun[1,], ci.high = ci.Tmean_AprMayJun[3,])
# ggplot() +
#   geom_ribbon(data = ci.Tmean_AprMayJun.df, aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") +
#   geom_line(data = ci.Tmean_AprMayJun.df, aes(x = Tmean_AprMayJun, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) +
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "Tmean_AprMayJun")]), aes(x = Tmean_AprMayJun))
#
#
# #effect of Tmean_SepOct
# Tmean_SepOctrng <- range(grow_train$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# growthpredictionTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct))
# tsofun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
#   for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
#     growthpredictionTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#       plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#   }
#   growthpredictionTmeanSepOct
# }
# TmeanSepOct_prediction <- tsofun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug,
#                                                   Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
# TmeanSepOct_prediction_tr <- exp(TmeanSepOct_prediction)
# ci.Tmean_SepOct <- apply(TmeanSepOct_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_SepOct.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_SepOct[2,], ci.low = ci.Tmean_SepOct[1,], ci.high = ci.Tmean_SepOct[3,])
# ggplot() +
#   geom_ribbon(data = ci.Tmean_SepOct.df, aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") +
#   geom_line(data = ci.Tmean_SepOct.df, aes(x = Tmean_SepOct, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) +
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "Tmean_SepOct")]), aes(x = Tmean_SepOct))
#
#
# #effect of DIA_prev
# sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
# size <- seq(sizerng[1], sizerng[2], by = 0.5)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# Tmean_AprMayJun <- mean(grow_train$Precip_AprMayJun)
# growthpredictionsize <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size))
# sizefun_plotdatainterval<-function(size,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_SepOct,Tmean_AprMayJun){
#   for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
#     growthpredictionsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#       plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#       plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#       plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#   }
#   growthpredictionsize
# }
# size_prediction <- sizefun_plotdatainterval(size = size, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_SepOct = Tmean_SepOct, Tmean_AprMayJun = Tmean_AprMayJun)
# size_prediction_tr <- exp(size_prediction)
# ci.size <- apply(size_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.size.df <- data.frame(size = size, median = ci.size[2,], ci.low = ci.size[1,], ci.high = ci.size[3,])
# ggplot() +
#   geom_ribbon(data = ci.size.df, aes(x = size, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") +
#   geom_line(data = ci.size.df, aes(x = size, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) +
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "DIA_prev")]), aes(x = size))
#
#
# #Precip_JulAug and Tmean_AprMayJun
# Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# Precip_JulAug_range <- quantile(grow_train$Precip_JulAug, c(0.2, 0.8))
# growthpredictionTmeanAprMayJun_highPrecipJulAug <- growthpredictionTmeanAprMayJun_lowPrecipJulAug <- growthpredictionTmeanAprMayJun_midPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
#   growthpredictionTmeanAprMayJun_highPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug_range[2] +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug_range[2]*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionTmeanAprMayJun_midPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionTmeanAprMayJun_lowPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug_range[1] +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug_range[1]*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Tmean_AprMayJun_prediction_trlow <- exp(growthpredictionTmeanAprMayJun_lowPrecipJulAug)
# Tmean_AprMayJun_prediction_trmid <- exp(growthpredictionTmeanAprMayJun_midPrecipJulAug)
# Tmean_AprMayJun_prediction_trhigh <- exp(growthpredictionTmeanAprMayJun_highPrecipJulAug)
# ci.Tmean_AprMayJunhigh <- apply(Tmean_AprMayJun_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_AprMayJunmid <- apply(Tmean_AprMayJun_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunlow <- apply(Tmean_AprMayJun_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highPrecip_JulAug")
# ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midPrecip_JulAug")
# ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowPrecip_JulAug")
# Tmean_AprMayJun_Precip_JulAugint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
# ggplot(data = Tmean_AprMayJun_Precip_JulAugint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Tmean_AprMayJun_Precip_JulAugint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
# #Precip_JulAug and Tmean_SepOct
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# Tmean_SepOct_range <- quantile(grow_train$Tmean_SepOct, c(0.2, 0.8))
# growthpredictionPrecipJulAug_highTmeanSepOct <- growthpredictionPrecipJulAug_lowTmeanSepOct <- growthpredictionPrecipJulAug_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[2]
#
#   growthpredictionPrecipJulAug_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[1]
# }
# Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowTmeanSepOct)
# Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midTmeanSepOct)
# Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highTmeanSepOct)
# ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highTmean_SepOct")
# ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midTmean_SepOct")
# ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowTmean_SepOct")
# Precip_JulAug_Tmean_SepOctint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
# ggplot(data = Precip_JulAug_Tmean_SepOctint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Precip_JulAug_Tmean_SepOctint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
#   geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #Precip_JulAug and Precip_NovDecJanFebMar
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# Precip_NovDecJanFebMar_range <- quantile(grow_train$Precip_NovDecJanFebMar, c(0.2, 0.8))
# growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar <- growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar <- growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar_range[2]*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar_range[1]*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Precip_JulAug_predictionpndjfm_trlow <- exp(growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar)
# Precip_JulAug_predictionpndjfm_trmid <- exp(growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar)
# Precip_JulAug_predictionpndjfm_trhigh <- exp(growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar)
# ci.Precip_JulAughigh_Precip_NovDecJanFebMar <- apply(Precip_JulAug_predictionpndjfm_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid_Precip_NovDecJanFebMar <- apply(Precip_JulAug_predictionpndjfm_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow_Precip_NovDecJanFebMar <- apply(Precip_JulAug_predictionpndjfm_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh_Precip_NovDecJanFebMar.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh_Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_JulAughigh_Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_JulAughigh_Precip_NovDecJanFebMar[3,], ci.group = "highPrecip_NovDecJanFebMar")
# ci.Precip_JulAugmid_Precip_NovDecJanFebMar.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid_Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_JulAugmid_Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_JulAugmid_Precip_NovDecJanFebMar[3,], ci.group = "midPrecip_NovDecJanFebMar")
# ci.Precip_JulAuglow_Precip_NovDecJanFebMar.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow_Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_JulAuglow_Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_JulAuglow_Precip_NovDecJanFebMar[3,], ci.group = "lowPrecip_NovDecJanFebMar")
# Precip_JulAug_Precip_NovDecJanFebMarint <- rbind(ci.Precip_JulAughigh_Precip_NovDecJanFebMar.df, ci.Precip_JulAugmid_Precip_NovDecJanFebMar.df, ci.Precip_JulAuglow_Precip_NovDecJanFebMar.df)
# ggplot(data = Precip_JulAug_Precip_NovDecJanFebMarint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Precip_JulAug_Precip_NovDecJanFebMarint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
#   geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #Precip_JulAug and tmp_norm
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictionPrecipJulAug_hightnorm <- growthpredictionPrecipJulAug_lowtnorm <- growthpredictionPrecipJulAug_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowtnorm)
# Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midtnorm)
# Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_hightnorm)
# ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "hightnorm")
# ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midtnorm")
# ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowtnorm")
# Precip_JulAug_tnormint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
# ggplot(data = Precip_JulAug_tnormint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Precip_JulAug_tnormint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
#   geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
#
# #Precip_JulAug and ppt_norm
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
# growthpredictionPrecipJulAug_highpnorm <- growthpredictionPrecipJulAug_lowpnorm <- growthpredictionPrecipJulAug_midpnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_highpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[2]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_midpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_lowpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[1]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowpnorm)
# Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midpnorm)
# Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highpnorm)
# ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highpnorm")
# ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midpnorm")
# ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowpnorm")
# Precip_JulAug_pnormint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
# ggplot(data = Precip_JulAug_pnormint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Precip_JulAug_pnormint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
#   geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #Precip_JulAug and size
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
# growthpredictionPrecipJulAug_highsize <- growthpredictionPrecipJulAug_lowsize <- growthpredictionPrecipJulAug_midsize <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipJulAug_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowsize)
# Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midsize)
# Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highsize)
# ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highsize")
# ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midsize")
# ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowsize")
# Precip_JulAug_sizeint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
# ggplot(data = Precip_JulAug_sizeint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Precip_JulAug_sizeint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
#   geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #ppt_norm and tmp_norm
# ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
# ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictionpnorm_hightnorm <- growthpredictionpnorm_lowtnorm <- growthpredictionpnorm_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm))
#
# for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
#   growthpredictionpnorm_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpnorm_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpnorm_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowtnorm)
# ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midtnorm)
# ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_hightnorm)
# ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normhigh[2,], ci.low = ci.ppt_normhigh[1,], ci.high = ci.ppt_normhigh[3,], ci.group = "hightnorm")
# ci.ppt_normmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normmid[2,], ci.low = ci.ppt_normmid[1,], ci.high = ci.ppt_normmid[3,], ci.group = "midtnorm")
# ci.ppt_normlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normlow[2,], ci.low = ci.ppt_normlow[1,], ci.high = ci.ppt_normlow[3,], ci.group = "lowtnorm")
# ppt_norm_pnormint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
# ggplot(data = ppt_norm_pnormint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = ppt_norm_pnormint, aes(x = ppt_norm, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
#   geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #ppt_norm and Precip_DecJanFeb interaction
# Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
# Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
# growthpredictionpptnorm_highPrecipNovDecJanFebMar <- growthpredictionpptnorm_lowPrecipNovDecJanFebMar <- growthpredictionpptnorm_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
#   growthpredictionpptnorm_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[2]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpptnorm_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpptnorm_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[1]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Precip_NovDecJanFebMar_prediction_trlow <- exp(growthpredictionpptnorm_lowPrecipNovDecJanFebMar)
# Precip_NovDecJanFebMar_prediction_trmid <- exp(growthpredictionpptnorm_midPrecipNovDecJanFebMar)
# Precip_NovDecJanFebMar_prediction_trhigh <- exp(growthpredictionpptnorm_highPrecipNovDecJanFebMar)
# ci.Precip_NovDecJanFebMarhigh <- apply(Precip_NovDecJanFebMar_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_NovDecJanFebMarmid <- apply(Precip_NovDecJanFebMar_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_NovDecJanFebMarlow <- apply(Precip_NovDecJanFebMar_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_NovDecJanFebMarhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "highppt_norm")
# ci.Precip_NovDecJanFebMarmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midppt_norm")
# ci.Precip_NovDecJanFebMarlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowppt_norm")
# Precip_NovDecJanFebMar_ppt_normint <- rbind(ci.Precip_NovDecJanFebMarhigh.df, ci.Precip_NovDecJanFebMarmid.df, ci.Precip_NovDecJanFebMarlow.df)
# ggplot(data = Precip_NovDecJanFebMar_ppt_normint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Precip_NovDecJanFebMar_ppt_normint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #ppt_norm and Tmean_AprMayJun interaction
# Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# ppt_norm <- mean(grow_train$ppt_norm)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
# growthpredictionpnorm_highTmeanAprMayJun <- growthpredictionpnorm_lowTmeanAprMayJun <- growthpredictionpnorm_midTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
#   growthpredictionpnorm_highTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[2]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpnorm_midTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpnorm_lowTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[1]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowTmeanAprMayJun)
# ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midTmeanAprMayJun)
# ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highTmeanAprMayJun)
# ci.Tmean_AprMayJunhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_AprMayJunmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highppt_norm")
# ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midppt_norm")
# ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowppt_norm")
# ppt_norm_Tmean_AprMayJunint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
# ggplot(data = ppt_norm_Tmean_AprMayJunint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = ppt_norm_Tmean_AprMayJunint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #ppt_norm and Tmean_SepOct interaction
# Tmean_SepOctrng <- range(grow_train$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# ppt_norm <- mean(grow_train$ppt_norm)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
# growthpredictionpnorm_highTmeanSepOct <- growthpredictionpnorm_lowTmeanSepOct <- growthpredictionpnorm_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
#   growthpredictionpnorm_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[2]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpnorm_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpnorm_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[1]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowTmeanSepOct)
# ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midTmeanSepOct)
# ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highTmeanSepOct)
# ci.Tmean_SepOcthigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_SepOctmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_SepOctlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_SepOcthigh.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_SepOcthigh[2,], ci.low = ci.Tmean_SepOcthigh[1,], ci.high = ci.Tmean_SepOcthigh[3,], ci.group = "highppt_norm")
# ci.Tmean_SepOctmid.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_SepOctmid[2,], ci.low = ci.Tmean_SepOctmid[1,], ci.high = ci.Tmean_SepOctmid[3,], ci.group = "midppt_norm")
# ci.Tmean_SepOctlow.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_SepOctlow[2,], ci.low = ci.Tmean_SepOctlow[1,], ci.high = ci.Tmean_SepOctlow[3,], ci.group = "lowppt_norm")
# ppt_norm_Tmean_SepOctint <- rbind(ci.Tmean_SepOcthigh.df, ci.Tmean_SepOctmid.df, ci.Tmean_SepOctlow.df)
# ggplot(data = ppt_norm_Tmean_AprMayJunint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = ppt_norm_Tmean_SepOctint, aes(x = Tmean_SepOct, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #ppt_norm and size
# ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
# ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
# growthpredictionpnorm_highsize <- growthpredictionpnorm_lowsize <- growthpredictionpnorm_midsize <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm))
#
# for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
#   growthpredictionpnorm_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpnorm_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionpnorm_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowsize)
# ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midsize)
# ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highsize)
# ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normhigh[2,], ci.low = ci.ppt_normhigh[1,], ci.high = ci.ppt_normhigh[3,], ci.group = "highsize")
# ci.ppt_normmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normmid[2,], ci.low = ci.ppt_normmid[1,], ci.high = ci.ppt_normmid[3,], ci.group = "midTmeansize")
# ci.ppt_normlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normlow[2,], ci.low = ci.ppt_normlow[1,], ci.high = ci.ppt_normlow[3,], ci.group = "lowTmeansize")
# ppt_norm_sizeint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
# ggplot(data = ppt_norm_sizeint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = ppt_norm_sizeint, aes(x = ppt_norm, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #tmp_norm and Precip_NovDecJanFebMar interaction
# Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
# Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictiontnorm_highPrecipNovDecJanFebMar <- growthpredictiontnorm_lowPrecipNovDecJanFebMar <- growthpredictiontnorm_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
#   growthpredictiontnorm_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowPrecipNovDecJanFebMar)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midPrecipNovDecJanFebMar)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highPrecipNovDecJanFebMar)
# ci.Precip_NovDecJanFebMarhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_NovDecJanFebMarmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_NovDecJanFebMarlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_NovDecJanFebMarhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "hightmp_norm")
# ci.Precip_NovDecJanFebMarmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midtmp_norm")
# ci.Precip_NovDecJanFebMarlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowtmp_norm")
# tmp_norm_Precip_NovDecJanFebMarint <- rbind(ci.Precip_NovDecJanFebMarhigh.df, ci.Precip_NovDecJanFebMarmid.df, ci.Precip_NovDecJanFebMarlow.df)
# ggplot(data = tmp_norm_Precip_NovDecJanFebMarint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = tmp_norm_Precip_NovDecJanFebMarint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #tmp_norm and Precip_JulAug interaction
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictiontnorm_highPrecipJulAug <- growthpredictiontnorm_lowPrecipJulAug <- growthpredictiontnorm_midPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictiontnorm_highPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_midPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_lowPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowPrecipJulAug)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midPrecipJulAug)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highPrecipJulAug)
# ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normhigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "hightnorm")
# ci.tmp_normmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midtnorm")
# ci.tmp_normlow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowtnorm")
# tmp_norm_Precip_JulAugint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = tmp_norm_Precip_JulAugint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = tmp_norm_Precip_JulAugint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #tmp_norm and Tmean_AprMayJun interaction
# Tmean_AprMayJun <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictiontnorm_highTmeanAprMayJun <- growthpredictiontnorm_lowTmeanAprMayJun <- growthpredictiontnorm_midTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
#   growthpredictiontnorm_highTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_midTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_lowTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowTmeanAprMayJun)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midTmeanAprMayJun)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highTmeanAprMayJun)
# ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "hightnorm")
# ci.tmp_normmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midtnorm")
# ci.tmp_normlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowtnorm")
# tmp_norm_Precip_JulAugint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = tmp_norm_Precip_JulAugint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = tmp_norm_Precip_JulAugint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #tmp_norm and Tmean_SepOct interaction
# Tmean_SepOct <- range(grow_train$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictiontnorm_highTmeanSepOct <- growthpredictiontnorm_lowTmeanSepOct <- growthpredictiontnorm_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
#   growthpredictiontnorm_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowTmeanSepOct)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midTmeanSepOct)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highTmeanSepOct)
# ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normhigh.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "hightnorm")
# ci.tmp_normmid.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midtnorm")
# ci.tmp_normlow.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowtnorm")
# tmp_norm_Tmean_SepOctint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = tmp_norm_Tmean_SepOctint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = tmp_norm_Tmean_SepOctint, aes(x = Tmean_SepOct, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #tmp_norm and size interaction
# tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
# tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
# growthpredictiontnorm_highsize <- growthpredictiontnorm_lowsize <- growthpredictiontnorm_midsize <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm))
#
# for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
#   growthpredictiontnorm_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictiontnorm_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowsize)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midsize)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highsize)
# ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highsize")
# ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midsize")
# ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowsize")
# tmp_norm_sizeint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = tmp_norm_sizeint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = tmp_norm_sizeint, aes(x = tmp_norm, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #Precip_NovDecJanFebMar and Tmean_DecJanFeb interaction
# Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
# Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# Tmean_SepOct_range <- quantile(grow_train$Tmean_SepOct, c(0.2, 0.8))
# growthpredictionPrecipNovDecJanFebMar_highTmeanSepOct <- growthpredictionPrecipNovDecJanFebMar_lowTmeanSepOct <- growthpredictionPrecipNovDecJanFebMar_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
#   growthpredictionPrecipNovDecJanFebMar_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[2]
#
#   growthpredictionPrecipNovDecJanFebMar_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipNovDecJanFebMar_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[1]
# }
# Precip_NovDecJanFebMar_prediction_trlow <- exp(growthpredictionPrecipNovDecJanFebMar_lowTmeanSepOct)
# Precip_NovDecJanFebMar_prediction_trmid <- exp(growthpredictionPrecipNovDecJanFebMar_midTmeanSepOct)
# Precip_NovDecJanFebMar_prediction_trhigh <- exp(growthpredictionPrecipNovDecJanFebMar_highTmeanSepOct)
# ci.Precip_NovDecJanFebMarhigh <- apply(Precip_NovDecJanFebMar_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_NovDecJanFebMarmid <- apply(Precip_NovDecJanFebMar_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_NovDecJanFebMarlow <- apply(Precip_NovDecJanFebMar_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_NovDecJanFebMarhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "highTmeanSepOct")
# ci.Precip_NovDecJanFebMarmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midTmeanSepOct")
# ci.Precip_NovDecJanFebMarlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowTmeanSepOct")
# Precip_NovDecJanFebMar_Tmean_SepOctint <- rbind(ci.Precip_NovDecJanFebMarhigh.df, ci.Precip_NovDecJanFebMarmid.df, ci.Precip_NovDecJanFebMarlow.df)
# ggplot(data = Precip_NovDecJanFebMar_Tmean_SepOctint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Precip_NovDecJanFebMar_Tmean_SepOctint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
#
# #Precip_NovDecJanFebMar and Tmean_AprMayJun interaction
# Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# Precip_NovDecJanFebMar_range <- quantile(grow_train$Precip_NovDecJanFebMar, c(0.2, 0.8))
# growthpredictionTmeanAprMayJun_highPrecipNovDecJanFebMar <- growthpredictionTmeanAprMayJun_lowPrecipNovDecJanFebMar <- growthpredictionTmeanAprMayJun_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
#   growthpredictionTmeanAprMayJun_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar_range[2]*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionTmeanAprMayJun_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionTmeanAprMayJun_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar_range[1]*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Tmean_AprMayJun_prediction_trlow <- exp(growthpredictionTmeanAprMayJun_lowPrecipNovDecJanFebMar)
# Tmean_AprMayJun_prediction_trmid <- exp(growthpredictionTmeanAprMayJun_midPrecipNovDecJanFebMar)
# Tmean_AprMayJun_prediction_trhigh <- exp(growthpredictionTmeanAprMayJun_highPrecipNovDecJanFebMar)
# ci.Tmean_AprMayJunhigh <- apply(Tmean_AprMayJun_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_AprMayJunmid <- apply(Tmean_AprMayJun_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunlow <- apply(Tmean_AprMayJun_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highPrecip_NovDecJanFebMar")
# ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midPrecip_NovDecJanFebMar")
# ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowPrecip_NovDecJanFebMar")
# Tmean_AprMayJun_Precip_NovDecJanFebMarint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
# ggplot(data = Tmean_AprMayJun_Precip_NovDecJanFebMarint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Tmean_AprMayJun_Precip_NovDecJanFebMarint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #Precip_NovDecJanFebMar and size interaction
# Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
# Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
# growthpredictionPrecipNovDecJanFebMar_highsize <- growthpredictionPrecipNovDecJanFebMar_lowsize <- growthpredictionPrecipNovDecJanFebMar_midsize <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar))
#
# for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
#   growthpredictionPrecipNovDecJanFebMar_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipNovDecJanFebMar_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionPrecipNovDecJanFebMar_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Precip_NovDecJanFebMar_prediction_trlow <- exp(growthpredictionPrecipNovDecJanFebMar_lowsize)
# Precip_NovDecJanFebMar_prediction_trmid <- exp(growthpredictionPrecipNovDecJanFebMar_midsize)
# Precip_NovDecJanFebMar_prediction_trhigh <- exp(growthpredictionPrecipNovDecJanFebMar_highsize)
# ci.Precip_NovDecJanFebMarhigh <- apply(Precip_NovDecJanFebMar_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_NovDecJanFebMarmid <- apply(Precip_NovDecJanFebMar_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_NovDecJanFebMarlow <- apply(Precip_NovDecJanFebMar_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_NovDecJanFebMarhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "highsize")
# ci.Precip_NovDecJanFebMarmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midsize")
# ci.Precip_NovDecJanFebMarlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowsize")
# Precip_NovDecJanFebMar_sizeint <- rbind(ci.Precip_NovDecJanFebMarhigh.df, ci.Precip_NovDecJanFebMarmid.df, ci.Precip_NovDecJanFebMarlow.df)
# ggplot(data = Precip_NovDecJanFebMar_sizeint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Precip_NovDecJanFebMar_sizeint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
#
# #Tmean_AprMayJun and Tmean_SepOct interaction
# Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# Tmean_SepOct_range <- quantile(grow_train$Tmean_SepOct, c(0.2, 0.8))
# growthpredictionTmeanAprMayJun_highTmeanSepOct <- growthpredictionTmeanAprMayJun_lowTmeanSepOct <- growthpredictionTmeanAprMayJun_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
#   growthpredictionTmeanAprMayJun_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[2]
#
#   growthpredictionTmeanAprMayJun_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionTmeanAprMayJun_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[1]
# }
# Tmean_AprMayJun_prediction_trlow <- exp(growthpredictionTmeanAprMayJun_lowTmeanSepOct)
# Tmean_AprMayJun_prediction_trmid <- exp(growthpredictionTmeanAprMayJun_midTmeanSepOct)
# Tmean_AprMayJun_prediction_trhigh <- exp(growthpredictionTmeanAprMayJun_highTmeanSepOct)
# ci.Tmean_AprMayJunhigh <- apply(Tmean_AprMayJun_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_AprMayJunmid <- apply(Tmean_AprMayJun_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunlow <- apply(Tmean_AprMayJun_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highTmeanSepOct")
# ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midTmeanSepOct")
# ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowTmeanSepOct")
# Tmean_AprMayJun_Tmean_SepOctint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
# ggplot(data = Tmean_AprMayJun_Tmean_SepOctint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Tmean_AprMayJun_Tmean_SepOctint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #Tmean_AprMayJun and size interaction
# Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
# growthpredictionTmeanAprMayJun_highsize <- growthpredictionTmeanAprMayJun_lowsize <- growthpredictionTmeanAprMayJun_midsize <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
#   growthpredictionTmeanAprMayJun_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionTmeanAprMayJun_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionTmeanAprMayJun_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# Tmean_AprMayJun_prediction_trlow <- exp(growthpredictionTmeanAprMayJun_lowsize)
# Tmean_AprMayJun_prediction_trmid <- exp(growthpredictionTmeanAprMayJun_midsize)
# Tmean_AprMayJun_prediction_trhigh <- exp(growthpredictionTmeanAprMayJun_highsize)
# ci.Tmean_AprMayJunhigh <- apply(Tmean_AprMayJun_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_AprMayJunmid <- apply(Tmean_AprMayJun_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunlow <- apply(Tmean_AprMayJun_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highsize")
# ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midsize")
# ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowsize")
# Tmean_AprMayJun_sizeint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
# ggplot(data = Tmean_AprMayJun_sizeint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# ggplot(data = Tmean_AprMayJun_sizeint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
#   scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
#   geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
#   scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
#
# #Tmean_AprMayJun and size
# Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 7.5)
# Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# size <- mean(grow_train$DIA_prev)
# size_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
# growthpredictionsize_highTmeanAprMayJun <- growthpredictionsize_lowTmeanAprMayJun <- growthpredictionsize_midTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
#   growthpredictionsize_highTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionsize_midTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_lowTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# size_prediction_trlow <- exp(growthpredictionsize_lowTmeanDecJanFeb)
# size_prediction_trmid <- exp(growthpredictionsize_midTmeanDecJanFeb)
# size_prediction_trhigh <- exp(growthpredictionsize_highTmeanDecJanFeb)
# ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highTmeanDecJanFeb")
# ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midTmeanDecJanFeb")
# ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowTmeanDecJanFeb")
# size_Tmean_DecJanFebint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
# ggplot(data = size_Tmean_DecJanFebint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #size and Tmean_SepOct
# Tmean_SepOctrng <- range(grow_train$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_SepOct <- seq(sizerng[1], sizerng[2], by = 7.5)
# size <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
# Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
# size_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
# growthpredictionsize_highTmeanSepOct <- growthpredictionsize_lowTmeanSepOct <- growthpredictionsize_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct))
#
# for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
#   growthpredictionsize_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size_range[2]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size_range[2]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
#
#   growthpredictionsize_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
#     plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size_range[1]*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size_range[1]*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +
#     plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct +
#     plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
# }
# size_prediction_trlow <- exp(growthpredictionsize_lowTmeanJulAug)
# size_prediction_trmid <- exp(growthpredictionsize_midTmeanJulAug)
# size_prediction_trhigh <- exp(growthpredictionsize_highTmeanJulAug)
# ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highTmeanJulAug")
# ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midTmeanJulAug")
# ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowTmeanJulAug")
# size_Tmean_JulAugint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
# ggplot(data = size_Tmean_JulAugint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #size and Precip_JulAug
# sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
# size <- seq(sizerng[1], sizerng[2], by = 7.5)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Precip_JulAug_range <- quantile(grow_train$Precip_JulAug, c(0.2, 0.8))
# growthpredictionsize_highPrecipJulAug <- growthpredictionsize_lowPrecipJulAug <- growthpredictionsize_midPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size))
#
# for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
#   growthpredictionsize_highPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug_range[2]*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug_range[2]*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug_range[2] +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_midPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_lowPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug_range[1]*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug_range[1]*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug_range[1] +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
# }
# size_prediction_trlow <- exp(growthpredictionsize_lowPrecipJulAug)
# size_prediction_trmid <- exp(growthpredictionsize_midPrecipJulAug)
# size_prediction_trhigh <- exp(growthpredictionsize_highPrecipJulAug)
# ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highPrecipJulAug")
# ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midPrecipJulAug")
# ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowPrecipJulAug")
# size_Precip_JulAugint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
# ggplot(data = size_Precip_JulAugint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
#
# #size and ppt_norm
# sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
# size <- seq(sizerng[1], sizerng[2], by = 7.5)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
# growthpredictionsize_highpnorm <- growthpredictionsize_lowpnorm <- growthpredictionsize_midpnorm <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size))
#
# for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
#   growthpredictionsize_highpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_midpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_lowpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
# }
# size_prediction_trlow <- exp(growthpredictionsize_lowpnorm)
# size_prediction_trmid <- exp(growthpredictionsize_midpnorm)
# size_prediction_trhigh <- exp(growthpredictionsize_highpnorm)
# ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highpnorm")
# ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midpnorm")
# ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowpnorm")
# size_pnormint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
# ggplot(data = size_pnormint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #size and tmp_norm
# sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
# size <- seq(sizerng[1], sizerng[2], by = 7.5)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictionsize_hightnorm <- growthpredictionsize_lowtnorm <- growthpredictionsize_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size))
#
# for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
#   growthpredictionsize_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[2] + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm_range[2] + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[1] + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm_range[1] + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
# }
# size_prediction_trlow <- exp(growthpredictionsize_lowtnorm)
# size_prediction_trmid <- exp(growthpredictionsize_midtnorm)
# size_prediction_trhigh <- exp(growthpredictionsize_hightnorm)
# ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "hightnorm")
# ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midtnorm")
# ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowtnorm")
# size_tnormint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
# ggplot(data = size_tnormint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
# #size and Precip_DecJanFeb
# sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
# size <- seq(sizerng[1], sizerng[2], by = 7.5)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Precip_DecJanFeb_range <- quantile(grow_train$Precip_DecJanFeb, c(0.2, 0.8))
# growthpredictionsize_highPrecipDecJanFeb <- growthpredictionsize_lowPrecipDecJanFeb <- growthpredictionsize_midPrecipDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size))
#
# for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
#   growthpredictionsize_highPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[2]*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[2]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[2]*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[2]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[2]*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_midPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
#
#   growthpredictionsize_lowPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[1]*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[1]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[1]*Tmean_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[1]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[1]*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
# }
# size_prediction_trlow <- exp(growthpredictionsize_lowPrecipDecJanFeb)
# size_prediction_trmid <- exp(growthpredictionsize_midPrecipDecJanFeb)
# size_prediction_trhigh <- exp(growthpredictionsize_highPrecipDecJanFeb)
# ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highPrecipDecJanFeb")
# ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midPrecipDecJanFeb")
# ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowPrecipDecJanFeb")
# size_PrecipDecJanFebint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
# ggplot(data = size_PrecipDecJanFebint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) +
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#
#
