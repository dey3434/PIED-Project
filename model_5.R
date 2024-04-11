### Stan models for PIAL growth using technique in Ogle et al., Ecology Letters to test for climate effects
## Sharmila Dey
# 22 June 2020
# setwd("/home/rstudio")
#load(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_grow_coef2.rda"))
#readRDS("log_normal_fg_7_24_20.RDS")
library(rstan)
options(mc.cores = parallel::detectCores())
library(parallel) 
library(mcmcplots) ; library(lattice) ; library(MASS)
library(lme4) ; library(nlme) ; library(splines); library(MCMCpack)
library(ggplot2)
library(caret) ; library(tidyverse)
library(bayesplot)
library(here)
library(gifski)
library(maps)
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

grow.new$Tmean_AprMayJun_l <- grow.new$Tmean_AprMayJun
grow.new$Tmean_SepOct_l <- grow.new$Tmean_SepOct
grow.new$Precip_NovDecJanFebMar_l <- grow.new$Precip_NovDecJanFebMar
grow.new$Precip_JulAug_l <- grow.new$Precip_JulAug

grow.monsoon<-na.omit(grow.new) %>% 
  mutate_at(scale, .vars = vars(Precip_JulAug, Precip_NovDecJanFebMar, Tmean_AprMayJun, Tmean_SepOct, ppt_norm, tmp_norm)) %>%
  arrange(PLOT,SUBP,name) %>%
  mutate(PlotCD=as.numeric(factor(ST_PLT, levels = unique(ST_PLT))),treeCD=as.numeric(factor(name,levels=unique(name))),
         growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2)) %>%
  group_by(PlotCD)%>%
    mutate_at(scale, .vars = vars(Precip_JulAug_l, Precip_NovDecJanFebMar_l, Tmean_AprMayJun_l, Tmean_SepOct_l)) %>%
    ungroup()
  

set.seed(2023)
split=0.20
trainIndex <- createDataPartition(grow.monsoon$name, p=split, list=FALSE)
grow_test <- grow.monsoon[trainIndex,]
grow_train <- grow.monsoon[-trainIndex,]



xG<-as.matrix(cbind(grow_train$ppt_norm, grow_train$tmp_norm, grow_train$Precip_JulAug, grow_train$Precip_NovDecJanFebMar, 
                    grow_train$Tmean_AprMayJun, grow_train$Tmean_SepOct,
                    grow_train$Precip_JulAug_l, grow_train$Precip_NovDecJanFebMar_l,
                    grow_train$Tmean_AprMayJun_l, grow_train$Tmean_SepOct_l, grow_train$DIA_prev, 
                    grow_train$ppt_norm*grow_train$tmp_norm, grow_train$ppt_norm*grow_train$Precip_JulAug, 
                    grow_train$ppt_norm*grow_train$Precip_NovDecJanFebMar, grow_train$ppt_norm*grow_train$Tmean_AprMayJun,
                    grow_train$ppt_norm*grow_train$Tmean_SepOct,grow_train$ppt_norm*grow_train$DIA_prev,
                    grow_train$tmp_norm*grow_train$Precip_JulAug, grow_train$tmp_norm*grow_train$Precip_NovDecJanFebMar,
                    grow_train$tmp_norm*grow_train$Tmean_AprMayJun, grow_train$tmp_norm*grow_train$Tmean_SepOct, 
                    grow_train$tmp_norm*grow_train$DIA_prev, grow_train$Precip_JulAug*grow_train$Precip_NovDecJanFebMar,
                    grow_train$Precip_JulAug*grow_train$Tmean_AprMayJun, grow_train$Precip_JulAug*grow_train$Tmean_SepOct,
                    grow_train$Precip_JulAug*grow_train$DIA_prev, grow_train$Precip_NovDecJanFebMar*grow_train$Tmean_AprMayJun, 
                    grow_train$Precip_NovDecJanFebMar*grow_train$Tmean_SepOct, grow_train$Precip_NovDecJanFebMar*grow_train$DIA_prev, 
                    grow_train$Tmean_AprMayJun*grow_train$Tmean_SepOct, grow_train$Tmean_AprMayJun*grow_train$DIA_prev,
                    grow_train$Tmean_SepOct*grow_train$DIA_prev,
                    grow_train$Precip_JulAug_l*grow_train$ppt_norm, grow_train$Precip_JulAug_l*grow_train$tmp_norm,
                    grow_train$Precip_JulAug_l*grow_train$Precip_NovDecJanFebMar_l, grow_train$Precip_JulAug_l*grow_train$Tmean_AprMayJun_l,
                    grow_train$Precip_JulAug_l*grow_train$Tmean_SepOct_l, grow_train$Precip_JulAug_l*grow_train$DIA_prev,
                    grow_train$Precip_NovDecJanFebMar_l*grow_train$ppt_norm, grow_train$Precip_NovDecJanFebMar_l*grow_train$tmp_norm,
                    grow_train$Precip_NovDecJanFebMar_l*grow_train$Tmean_AprMayJun_l, grow_train$Precip_NovDecJanFebMar_l*grow_train$Tmean_SepOct_l,
                    grow_train$Precip_NovDecJanFebMar_l*grow_train$DIA_prev,
                    grow_train$Tmean_AprMayJun_l*grow_train$ppt_norm, grow_train$Tmean_AprMayJun_l*grow_train$tmp_norm,
                    grow_train$Tmean_AprMayJun_l*grow_train$Tmean_SepOct_l, grow_train$Tmean_AprMayJun_l*grow_train$DIA_prev,
                    grow_train$Tmean_SepOct_l*grow_train$ppt_norm, grow_train$Tmean_SepOct_l*grow_train$tmp_norm,
                    grow_train$Tmean_SepOct_l*grow_train$DIA_prev))
xGtest<-as.matrix(cbind(grow_test$ppt_norm, grow_test$tmp_norm, grow_test$Precip_JulAug, grow_test$Precip_NovDecJanFebMar, 
                        grow_test$Tmean_AprMayJun, grow_test$Tmean_SepOct,
                        grow_test$Precip_JulAug_l, grow_test$Precip_NovDecJanFebMar_l,
                        grow_test$Tmean_AprMayJun_l, grow_test$Tmean_SepOct_l, grow_test$DIA_prev, 
                        grow_test$ppt_norm*grow_test$tmp_norm, grow_test$ppt_norm*grow_test$Precip_JulAug, 
                        grow_test$ppt_norm*grow_test$Precip_NovDecJanFebMar, grow_test$ppt_norm*grow_test$Tmean_AprMayJun,
                        grow_test$ppt_norm*grow_test$Tmean_SepOct,grow_test$ppt_norm*grow_test$DIA_prev,
                        grow_test$tmp_norm*grow_test$Precip_JulAug, grow_test$tmp_norm*grow_test$Precip_NovDecJanFebMar,
                        grow_test$tmp_norm*grow_test$Tmean_AprMayJun, grow_test$tmp_norm*grow_test$Tmean_SepOct, 
                        grow_test$tmp_norm*grow_test$DIA_prev, grow_test$Precip_JulAug*grow_test$Precip_NovDecJanFebMar,
                        grow_test$Precip_JulAug*grow_test$Tmean_AprMayJun, grow_test$Precip_JulAug*grow_test$Tmean_SepOct,
                        grow_test$Precip_JulAug*grow_test$DIA_prev, grow_test$Precip_NovDecJanFebMar*grow_test$Tmean_AprMayJun, 
                        grow_test$Precip_NovDecJanFebMar*grow_test$Tmean_SepOct, grow_test$Precip_NovDecJanFebMar*grow_test$DIA_prev, 
                        grow_test$Tmean_AprMayJun*grow_test$Tmean_SepOct, grow_test$Tmean_AprMayJun*grow_test$DIA_prev,
                        grow_test$Tmean_SepOct*grow_test$DIA_prev,
                        grow_test$Precip_JulAug_l*grow_test$ppt_norm, grow_test$Precip_JulAug_l*grow_test$tmp_norm,
                        grow_test$Precip_JulAug_l*grow_test$Precip_NovDecJanFebMar_l, grow_test$Precip_JulAug_l*grow_test$Tmean_AprMayJun_l,
                        grow_test$Precip_JulAug_l*grow_test$Tmean_SepOct_l, grow_test$Precip_JulAug_l*grow_test$DIA_prev,
                        grow_test$Precip_NovDecJanFebMar_l*grow_test$ppt_norm, grow_test$Precip_NovDecJanFebMar_l*grow_test$tmp_norm,
                        grow_test$Precip_NovDecJanFebMar_l*grow_test$Tmean_AprMayJun_l, grow_test$Precip_NovDecJanFebMar_l*grow_test$Tmean_SepOct_l,
                        grow_test$Precip_NovDecJanFebMar_l*grow_test$DIA_prev,
                        grow_test$Tmean_AprMayJun_l*grow_test$ppt_norm, grow_test$Tmean_AprMayJun_l*grow_test$tmp_norm,
                        grow_test$Tmean_AprMayJun_l*grow_test$Tmean_SepOct_l, grow_test$Tmean_AprMayJun_l*grow_test$DIA_prev,
                        grow_test$Tmean_SepOct_l*grow_test$ppt_norm, grow_test$Tmean_SepOct_l*grow_test$tmp_norm,
                        grow_test$Tmean_SepOct_l*grow_test$DIA_prev))
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


sink("stancode/model_5.stan")
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



csvfiles <- here::here("results", paste0("megamodel_", 1:3, ".csv"))

if (all(file.exists(csvfiles))) {
  fit_grow <- read_stan_csv(csvfiles, col_major = TRUE) 
} else {
  fit_grow <- stan(file = 'stancode/model_5.stan', data = pied_dat, 
                   iter = 5000,
                   warmup = 1000,
                   chains = 3, cores = 8, 
                   sample_file = here::here("results", "megamodel"))
}

# chain1 <- rstan::read_stan_csv("/home/rstudio/megamodel_1.csv")
# chain2 <- rstan::read_stan_csv("/home/rstudio/megamodel_2.csv")
# chain3 <- rstan::read_stan_csv("/home/rstudio/megamodel_3.csv")

# fit_grow <- rstan::read_stan_csv(csvfiles)
# saveRDS(fit_grow, file = "megamodel.RDS")
summary<-summary(fit_grow)
summary

# Updated code below -----
fit_grow_df <- as.data.frame(fit_grow)
plotdata<-select(fit_grow_df,"yrep[1]":"yrep[8780]")
plotdatainterval<-select(fit_grow_df, "u_beta[1]":paste0("u_beta[", ncol(xG), "]"))
# plotdatainterval<-select(fit_grow_df, "u_beta[1]":"u_beta[28]")
colnames(plotdatainterval) <- c("u_beta_ppt_norm", "u_beta_tmp_norm", "u_beta_Precip_JulAug", "u_beta_Precip_NovDecJanFebMar", 
                                "u_beta_Tmean_AprMayJun", "u_beta_Tmean_JulAug", "u_beta_Tmean_SepOct", "u_beta_Tmean_NovDecJanFebMar",
                                "u_beta_DIA_prev", "u_beta_ppt_norm_tmp_norm", "u_beta_ppt_norm_Precip_JulAug", 
                                "u_beta_ppt_norm_Precip_NovDecJanFebMar", "u_beta_ppt_norm_Tmean_AprMayJun",
                                "u_beta_ppt_norm_Tmean_JulAug", "u_beta_ppt_norm_Tmean_SepOct",
                                "u_beta_ppt_norm_Tmean_NovDecJanFebMar", "u_beta_ppt_norm_DIA_prev",
                                "u_beta_tmp_norm_Precip_JulAug", "u_beta_tmp_norm_Precip_NovDecJanFebMar",
                                "u_beta_tmp_norm_Tmean_AprMayJun", "u_beta_tmp_norm_Tmean_JulAug",
                                "u_beta_tmp_norm_Tmean_SepOct", "u_beta_tmp_norm_Tmean_NovDecJanFebMar",
                                "u_beta_tmp_norm_DIA_prev", "u_beta_Precip_JulAug_Precip_NovDecJanFebMar",
                                "u_beta_Precip_JulAug_Tmean_AprMayJun", "u_beta_Precip_JulAug_Tmean_JulAug",
                                "u_beta_Precip_JulAug_Tmean_SepOct", "u_beta_Precip_JulAug_Tmean_NovDecJanFebMar",
                                "u_beta_Precip_JulAug_DIA_prev", "u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun", 
                                "u_beta_Precip_NovDecJanFebMar_Tmean_JulAug", "u_beta_Precip_NovDecJanFebMar_Tmean_SepOct", 
                                "u_beta_Precip_NovDecJanFebMar_Tmean_NovDecJanFebMar", "u_beta_Precip_NovDecJanFebMar_DIA_prev", 
                                "u_beta_Tmean_AprMayJun_Tmean_JulAug", "u_beta_Tmean_AprMayJun_Tmean_SepOct", 
                                "u_beta_Tmean_AprMayJun_Tmean_NovDecJanFebMar", "u_beta_Tmean_AprMayJun_DIA_prev", 
                                "u_beta_Tmean_JulAug_Tmean_SepOct", "u_beta_Tmean_JulAug_Tmean_NovDecJanFebMar", 
                                "u_beta_Tmean_JulAug_DIA_prev", "u_beta_Tmean_SepOct_Tmean_NovDecJanFebMar", 
                                "u_beta_Tmean_SepOct_DIA_prev", "u_beta_Tmean_NovDecJanFebMar_DIA_prev")

# get summaries of plotdatainterval:

df <- reshape2::melt(plotdatainterval)
beta.summaries <- df %>% group_by(variable) %>%
  summarise(mean = mean(value),
            ci.lo = quantile(value, 0.025),
            ci.hi = quantile(value, 0.975))

beta.summaries$allpos <- ifelse(beta.summaries$mean > 0 &  beta.summaries$ci.lo > 0 &  beta.summaries$ci.hi > 0, "yes", "no")
beta.summaries$allneg <- ifelse(beta.summaries$mean <= 0 &  beta.summaries$ci.lo <= 0 &  beta.summaries$ci.hi <= 0, "yes", "no")
beta.summaries$significant <- ifelse(beta.summaries$allpos == "yes" | beta.summaries$allneg == "yes", "significant", "not significant")
write.csv(beta.summaries, "betasummaries_model5_threechain_PIED.csv", row.names = FALSE)


####
#Validation-- This will be in a separate File
ext_fit <- rstan::extract(fit_grow)
yrep <- ext_fit$yrep
#yrep <- exp(yrep)
mean.pred <- apply(ext_fit$yrep, 2, mean)
p.o.df <- data.frame(predicted = exp(mean.pred), observed = exp(grow_test$loggrowth), error = (exp(mean.pred) - exp(grow_test$loggrowth)))
meansqrd <- (mean(p.o.df$error))^2

##### 
save(p.o.df, file = here::here("results", "model-5-pred-obs.RData"))

# ggplot(p.o.df, aes(predicted, observed)) + geom_point(alpha = 0.1) + geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
#   ylim(0, 10) + xlim(0,10)

####
p_pred_vs_observed <- ggplot(p.o.df, aes(predicted, observed)) + 
  geom_point(alpha = 0.1) + 
  geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
  ylim(0, 10) + xlim(0,10)
p_pred_vs_observed
ggsave(here::here("images", "model_5", "pred_vs_observed.png"), p_pred_vs_observed)


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

save(ll, r_eff, leaveoneout, file = here::here("results", "model-5-loo.RData"))

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

ggsave(here::here("images", "model_5", "ppc-plot-1.png"),
       ppc1, width = 16/3, height = 9)
ggsave(here::here("images", "model_5", "ppc-plot-2.png"),
       ppc2, width = 16/3, height = 9)
ggsave(here::here("images", "model_5", "ppc-plot-3.png"),
       ppc3, width = 16/3, height = 9)


# # Commented out 03/10/2023 by JRT

# ppc_dens_overlay(yGtest, as.matrix(plotdata))
# 
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
#         print(p)
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
# 
# #MCMC Intervals Plots
# mcmc_intervals(plotdatainterval, prob = 0.5)
# pdf("plotdatainterval_mcmc_intervals.pdf", height = 6, width = 10) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
# mcmc_intervals(plotdatainterval, prob = 0.5) # put the code that makes one of the plots in here
# dev.off() # "device off" tells R to stop printing stuff to the pdf
# 
# 
# #MCMC Area Plot
# pdf("plotdatainterval_mcmc_areas.pdf", height = 8, width = 20) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
# color_scheme_set("purple")
# mcmc_areas(plotdatainterval, prob = 0.8)
# dev.off() # "device off" tells R to stop printing stuff to the pdf
# 
# #MCMC Traces
# pdf("plotdatainterval_mcmc_traces.pdf", height = 6, width = 30) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
# color_scheme_set("mix-blue-red")
# mcmc_trace(plotdatainterval,
#            facet_args = list(nrow = 2))
# #pars = c("alpha", "sigma")
# dev.off() # "device off" tells R to stop printing stuff to the pdf
# 
# mcmcplot(As.mcmc.list(fit_grow))
# 
# 
# png("plotdatainterval.png", height = 20, width = 24, units = "in", res = 200)
# pairs.panels(as.matrix(plotdatainterval))
# dev.off()
# 
# 
# 
# #Effects Plots
# mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                panel.background = element_blank(), axis.line = element_line(colour = "black"),
#                legend.text=element_text(size=11),legend.title=element_text(size=12),
#                legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
#                axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
#                axis.line.x = element_line(color="black", size = 0.3),
#                axis.line.y = element_line(color="black", size = 0.3))
# 
# #effect of ppt_norm
# ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for ppt_normrng
# ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.25)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# growthpredictionpptnorm <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm))
# pfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_DecJanFeb,Tmean_JulAug,Tmean_DecJanFeb){
#   for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
#     growthpredictionpptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   }
#   growthpredictionpptnorm
# }
# ppt_norm_prediction <- pfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_DecJanFeb = Precip_DecJanFeb, Tmean_JulAug = Tmean_JulAug, Tmean_DecJanFeb = Tmean_DecJanFeb)
# ppt_norm_prediction_tr <- exp(ppt_norm_prediction)
# ci.ppt_norm <- apply(ppt_norm_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.ppt_norm.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_norm[2,], ci.low = ci.ppt_norm[1,], ci.high = ci.ppt_norm[3,])
# ggplot() + 
#   geom_ribbon(data = ci.ppt_norm.df, aes(x = ppt_norm, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
#   geom_line(data = ci.ppt_norm.df, aes(x = ppt_norm, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "tmp_norm")]), aes(x = ppt_norm, color = tmp_norm))
# 
# 
# #effect of tmp_norm
# tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
# tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.35)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# growthpredictiontmpnorm <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm))
# tfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_DecJanFeb,Tmean_JulAug,Tmean_DecJanFeb){
#   for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
#     growthpredictiontmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   }
#   growthpredictiontmpnorm
# }
# tmp_norm_prediction <- tfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_DecJanFeb = Precip_DecJanFeb, Tmean_JulAug = Tmean_JulAug, Tmean_DecJanFeb = Tmean_DecJanFeb)
# tmp_norm_prediction_tr <- exp(tmp_norm_prediction)
# ci.tmp_norm <- apply(tmp_norm_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_norm.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_norm[2,], ci.low = ci.tmp_norm[1,], ci.high = ci.tmp_norm[3,])
# ggplot() + 
#   geom_ribbon(data = ci.tmp_norm.df, aes(x = tmp_norm, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
#   geom_line(data = ci.tmp_norm.df, aes(x = tmp_norm, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "tmp_norm")]), aes(x = ppt_norm, color = tmp_norm))
# 
# #effect of Precip_JulAug
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# growthpredictionPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
# pjafun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_DecJanFeb,Tmean_JulAug,Tmean_DecJanFeb){
#   for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#     growthpredictionPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   }
#   growthpredictionPrecipJulAug
# }
# PrecipJulAug_prediction <- pjafun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_DecJanFeb = Precip_DecJanFeb, Tmean_JulAug = Tmean_JulAug, Tmean_DecJanFeb = Tmean_DecJanFeb)
# PrecipJulAug_prediction_tr <- exp(PrecipJulAug_prediction)
# ci.Precip_JulAug <- apply(PrecipJulAug_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAug.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAug[2,], ci.low = ci.Precip_JulAug[1,], ci.high = ci.Precip_JulAug[3,])
# ggplot() + 
#   geom_ribbon(data = ci.Precip_JulAug.df, aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
#   geom_line(data = ci.Precip_JulAug.df, aes(x = Precip_JulAug, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "Precip_JulAug")]), aes(x = Precip_JulAug, color = ppt_norm))
# 
# 
# #effect of Precip_DecJanFeb
# Precip_DecJanFebrng <- range(grow_train$Precip_DecJanFeb,na.rm = TRUE) #setting range for tmp_normrng
# Precip_DecJanFeb <- seq(Precip_DecJanFebrng[1], Precip_DecJanFebrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# growthpredictionPrecipDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_Precip_DecJanFeb), length(Precip_DecJanFeb))
# pdjffun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_DecJanFeb,Tmean_JulAug,Tmean_DecJanFeb){
#   for(i in 1:length(plotdatainterval$u_beta_Precip_DecJanFeb)){
#     growthpredictionPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   }
#   growthpredictionPrecipDecJanFeb
# }
# PrecipDecJanFeb_prediction <- pdjffun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_DecJanFeb = Precip_DecJanFeb, Tmean_JulAug = Tmean_JulAug, Tmean_DecJanFeb = Tmean_DecJanFeb)
# PrecipDecJanFeb_prediction_tr <- exp(PrecipDecJanFeb_prediction)
# ci.Precip_DecJanFeb <- apply(PrecipDecJanFeb_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_DecJanFeb.df <- data.frame(Precip_DecJanFeb = Precip_DecJanFeb, median = ci.Precip_DecJanFeb[2,], ci.low = ci.Precip_DecJanFeb[1,], ci.high = ci.Precip_DecJanFeb[3,])
# ggplot() + 
#   geom_ribbon(data = ci.Precip_DecJanFeb.df, aes(x = Precip_DecJanFeb, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
#   geom_line(data = ci.Precip_DecJanFeb.df, aes(x = Precip_DecJanFeb, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "Precip_DecJanFeb")]), aes(x = Precip_DecJanFeb, color = ppt_norm))
# 
# #effect of Tmean_JulAug
# Tmean_JulAugrng <- range(grow_train$Tmean_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_JulAug <- seq(Tmean_JulAugrng[1], Tmean_JulAugrng[2], by = 0.3)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# growthpredictionTmeanJulAug <- matrix(NA, length(plotdatainterval$u_beta_Tmean_JulAug), length(Tmean_JulAug))
# tjafun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_DecJanFeb,Tmean_JulAug,Tmean_DecJanFeb){
#   for(i in 1:length(plotdatainterval$u_beta_Tmean_JulAug)){
#     growthpredictionTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   }
#   growthpredictionTmeanJulAug
# }
# Tmean_JulAug_prediction <- tjafun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_DecJanFeb = Precip_DecJanFeb, Tmean_JulAug = Tmean_JulAug, Tmean_DecJanFeb = Tmean_DecJanFeb)
# Tmean_JulAug_prediction_tr <- exp(Tmean_JulAug_prediction)
# ci.Tmean_JulAug <- apply(Tmean_JulAug_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_JulAug.df <- data.frame(Tmean_JulAug = Tmean_JulAug, median = ci.Tmean_JulAug[2,], ci.low = ci.Tmean_JulAug[1,], ci.high = ci.Tmean_JulAug[3,])
# ggplot() + 
#   geom_ribbon(data = ci.Tmean_JulAug.df, aes(x = Tmean_JulAug, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
#   geom_line(data = ci.Tmean_JulAug.df, aes(x = Tmean_JulAug, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "tmp_norm", "Tmean_JulAug")]), aes(x = Tmean_JulAug, color = tmp_norm))
# 
# 
# #effect of Tmean_DecJanFeb
# Tmean_DecJanFebrng <- range(grow_train$Tmean_DecJanFeb,na.rm = TRUE) #setting range for tmp_normrng
# Tmean_DecJanFeb <- seq(Tmean_DecJanFebrng[1], Tmean_DecJanFebrng[2], by = 0.3)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# growthpredictionTmeanDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_Tmean_DecJanFeb), length(Tmean_DecJanFeb))
# tdjffun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_DecJanFeb,Tmean_JulAug,Tmean_DecJanFeb){
#   for(i in 1:length(plotdatainterval$u_beta_Tmean_DecJanFeb)){
#     growthpredictionTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   }
#   growthpredictionTmeanDecJanFeb
# }
# Tmean_DecJanFeb_prediction <- tdjffun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_DecJanFeb = Precip_DecJanFeb, Tmean_JulAug = Tmean_JulAug, Tmean_DecJanFeb = Tmean_DecJanFeb)
# Tmean_DecJanFeb_prediction_tr <- exp(Tmean_DecJanFeb_prediction)
# ci.Tmean_DecJanFeb <- apply(Tmean_DecJanFeb_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Tmean_DecJanFeb.df <- data.frame(Tmean_DecJanFeb = Tmean_DecJanFeb, median = ci.Tmean_DecJanFeb[2,], ci.low = ci.Tmean_DecJanFeb[1,], ci.high = ci.Tmean_DecJanFeb[3,])
# ggplot() + 
#   geom_ribbon(data = ci.Tmean_DecJanFeb.df, aes(x = Tmean_DecJanFeb, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
#   geom_line(data = ci.Tmean_DecJanFeb.df, aes(x = Tmean_DecJanFeb, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "Tmean_DecJanFeb")]), aes(x = Tmean_DecJanFeb, color = ppt_norm))
# 
# 
# #effect of DIA_prev
# sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for sizerng
# size <- seq(sizerng[1], sizerng[2], by = 9)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# growthpredictionsize <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size))
# sizefun_plotdatainterval<-function(size,tmp_norm,ppt_norm,Precip_JulAug,Precip_DecJanFeb,Tmean_JulAug,Tmean_DecJanFeb){
#   for(i in 1:length(plotdatainterval$u_beta_Tmean_JulAug)){
#     growthpredictionsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   }
#   growthpredictionsize
# }
# size_prediction <- sizefun_plotdatainterval(size = size, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_DecJanFeb = Precip_DecJanFeb, Tmean_JulAug = Tmean_JulAug, Tmean_DecJanFeb = Tmean_DecJanFeb)
# size_prediction_tr <- exp(size_prediction)
# ci.size <- apply(size_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.size.df <- data.frame(size = size, median = ci.size[2,], ci.low = ci.size[1,], ci.high = ci.size[3,])
# ggplot() + 
#   geom_ribbon(data = ci.size.df, aes(x = size, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
#   geom_line(data = ci.size.df, aes(x = size, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
#   geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "DIA_prev")]), aes(x = size, color = ppt_norm))
# 
# ggplot(data = ci.size.df, aes(x = size, y = median)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
#   geom_line(color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#  
# 
# #Precip_JulAug and Tmean_JulAug
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Tmean_JulAug_range <- quantile(grow_train$Tmean_JulAug, c(0.2, 0.8))
# growthpredictionPrecipJulAug_highTmeanJulAug <- growthpredictionPrecipJulAug_lowTmeanJulAug <- growthpredictionPrecipJulAug_midTmeanJulAug <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
# 
#   for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#     growthpredictionPrecipJulAug_highTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug_range[2] + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug_range[2] + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug_range[2]*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug_range[2]*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug_range[2]*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug_range[2]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug_range[2]*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#     
#     growthpredictionPrecipJulAug_midTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#     
#     growthpredictionPrecipJulAug_lowTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#       plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug_range[1] + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#       plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#       plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#       plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug_range[1] + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#       plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug_range[1]*ppt_norm +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug_range[1]*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug_range[1]*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug_range[1]*Tmean_DecJanFeb +
#       plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug_range[1]*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#       plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#     }
# Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowTmeanJulAug)
# Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midTmeanJulAug)
# Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highTmeanJulAug)
# ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highTmean_JulAug")
# ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midTmean_JulAug")
# ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowTmean_JulAug")
# Precip_JulAug_Tmean_JulAugint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
# ggplot(data = Precip_JulAug_Tmean_JulAugint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# 
# #Precip_JulAug and Tmean_DecJanFeb
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Tmean_DecJanFeb_range <- quantile(grow_train$Tmean_DecJanFeb, c(0.2, 0.8))
# growthpredictionPrecipJulAug_highTmeanDecJanFeb <- growthpredictionPrecipJulAug_lowTmeanDecJanFeb <- growthpredictionPrecipJulAug_midTmeanDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_highTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[2] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[2]*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[2]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[2]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[2]*x
#   
#   growthpredictionPrecipJulAug_midTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionPrecipJulAug_lowTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[1] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[1]*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[1]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[1]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[1]*x
# }
# Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowTmeanDecJanFeb)
# Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midTmeanDecJanFeb)
# Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highTmeanDecJanFeb)
# ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highTmean_DecJanFeb")
# ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midTmean_DecJanFeb")
# ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowTmean_DecJanFeb")
# Precip_JulAug_Tmean_DecJanFebint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
# ggplot(data = Precip_JulAug_Tmean_DecJanFebint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# #Precip_JulAug and Precip_DecJanFeb
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Precip_DecJanFeb_range <- quantile(grow_train$Precip_DecJanFeb, c(0.2, 0.8))
# growthpredictionPrecipJulAug_highPrecipDecJanFeb <- growthpredictionPrecipJulAug_lowPrecipDecJanFeb <- growthpredictionPrecipJulAug_midPrecipDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_highPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[2]*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[2]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[2]*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[2]*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[2]*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionPrecipJulAug_midPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionPrecipJulAug_lowPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[1]*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[1]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[1]*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[1]*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[1]*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
# }
# Precip_JulAug_predictionpdjf_trlow <- exp(growthpredictionPrecipJulAug_lowPrecipDecJanFeb)
# Precip_JulAug_predictiondpjf_trmid <- exp(growthpredictionPrecipJulAug_midPrecipDecJanFeb)
# Precip_JulAug_predictionpdjf_trhigh <- exp(growthpredictionPrecipJulAug_highPrecipDecJanFeb)
# ci.Precip_JulAughigh_Precip_DecJanFeb <- apply(Precip_JulAug_predictionpdjf_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid_Precip_DecJanFeb <- apply(Precip_JulAug_predictiondpjf_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow_Precip_DecJanFeb <- apply(Precip_JulAug_predictionpdjf_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh_Precip_DecJanFeb.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh_Precip_DecJanFeb[2,], ci.low = ci.Precip_JulAughigh_Precip_DecJanFeb[1,], ci.high = ci.Precip_JulAughigh_Precip_DecJanFeb[3,], ci.group = "highPrecip_DecJanFeb")
# ci.Precip_JulAugmid_Precip_DecJanFeb.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid_Precip_DecJanFeb[2,], ci.low = ci.Precip_JulAugmid_Precip_DecJanFeb[1,], ci.high = ci.Precip_JulAugmid_Precip_DecJanFeb[3,], ci.group = "midPrecip_DecJanFeb")
# ci.Precip_JulAuglow_Precip_DecJanFeb.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow_Precip_DecJanFeb[2,], ci.low = ci.Precip_JulAuglow_Precip_DecJanFeb[1,], ci.high = ci.Precip_JulAuglow_Precip_DecJanFeb[3,], ci.group = "lowPrecip_DecJanFeb")
# Precip_JulAug_Precip_DecJanFebint <- rbind(ci.Precip_JulAughigh_Precip_DecJanFeb.df, ci.Precip_JulAugmid_Precip_DecJanFeb.df, ci.Precip_JulAuglow_Precip_DecJanFeb.df)
# ggplot(data = Precip_JulAug_Precip_DecJanFebint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# #Precip_JulAug and tmp_norm
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# ppt_norm <- mean(grow_train$ppt_norm)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictionPrecipJulAug_hightnorm <- growthpredictionPrecipJulAug_lowtnorm <- growthpredictionPrecipJulAug_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[2] + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm_range[2] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm_range[2] + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm_range[2] + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionPrecipJulAug_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionPrecipJulAug_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[1] + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm_range[1] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm_range[1] + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm_range[1] + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
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
# 
# 
# #ppt_norm and tmp_norm
# ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
# ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
# growthpredictionpnorm_hightnorm <- growthpredictionpnorm_lowtnorm <- growthpredictionpnorm_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
#   growthpredictionpnorm_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[2] + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm_range[2] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm_range[2] + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm_range[2] + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionpnorm_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionpnorm_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[1] + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm_range[1] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm_range[1] + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm_range[1] + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
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
# 
# 
# #ppt_norm and Precip_DecJanFeb interaction
# ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
# ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Precip_DecJanFeb_range <- quantile(grow_train$Precip_DecJanFeb, c(0.2, 0.8))
# growthpredictionpnorm_highPrecipDecJanFeb <- growthpredictionpnorm_lowPrecipDecJanFeb <- growthpredictionpnorm_midPrecipDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
#   growthpredictionpnorm_highPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[2]*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[2]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[2]*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[2]*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[2]*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionpnorm_midPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionpnorm_lowPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[1]*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[1]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[1]*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[1]*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[1]*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
# }
# ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowPrecipDecJanFeb)
# ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midPrecipDecJanFeb)
# ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highPrecipDecJanFeb)
# ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normhigh[2,], ci.low = ci.ppt_normhigh[1,], ci.high = ci.ppt_normhigh[3,], ci.group = "highPrecipDecJanFeb")
# ci.ppt_normmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normmid[2,], ci.low = ci.ppt_normmid[1,], ci.high = ci.ppt_normmid[3,], ci.group = "midPrecipDecJanFeb")
# ci.ppt_normlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normlow[2,], ci.low = ci.ppt_normlow[1,], ci.high = ci.ppt_normlow[3,], ci.group = "lowPrecipDecJanFeb")
# ppt_norm_Precip_DecJanFebint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
# ggplot(data = ppt_norm_Precip_DecJanFebint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# 
# #ppt_norm and Tmean_JulAug interaction
# ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
# ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Tmean_JulAug_range <- quantile(grow_train$Tmean_JulAug, c(0.2, 0.8))
# growthpredictionpnorm_highTmeanJulAug <- growthpredictionpnorm_lowTmeanJulAug <- growthpredictionpnorm_midTmeanJulAug <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
#   growthpredictionpnorm_highTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug_range[2] + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug_range[2] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug_range[2]*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug_range[2]*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug_range[2]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug_range[2]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug_range[2]*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionpnorm_midTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionpnorm_lowTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug_range[1] + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug_range[1] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug_range[1]*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug_range[1]*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug_range[1]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug_range[1]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug_range[1]*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
# }
# ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowTmeanJulAug)
# ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midTmeanJulAug)
# ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highTmeanJulAug)
# ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normhigh[2,], ci.low = ci.ppt_normhigh[1,], ci.high = ci.ppt_normhigh[3,], ci.group = "highTmeanJulAug")
# ci.ppt_normmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normmid[2,], ci.low = ci.ppt_normmid[1,], ci.high = ci.ppt_normmid[3,], ci.group = "midTmeanJulAug")
# ci.ppt_normlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normlow[2,], ci.low = ci.ppt_normlow[1,], ci.high = ci.ppt_normlow[3,], ci.group = "lowTmeanJulAug")
# ppt_norm_Tmean_JulAugint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
# ggplot(data = ppt_norm_Tmean_JulAugint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# #ppt_norm and Tmean_DecJanFeb interaction
# ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
# ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# tmp_norm <- mean(grow_train$tmp_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Tmean_DecJanFeb_range <- quantile(grow_train$Tmean_DecJanFeb, c(0.2, 0.8))
# growthpredictionpnorm_highTmeanDecJanFeb <- growthpredictionpnorm_lowTmeanDecJanFeb <- growthpredictionpnorm_midTmeanDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
#   growthpredictionpnorm_highTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[2] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[2]*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[2]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[2]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[2]*x
#   
#   growthpredictionpnorm_midTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionpnorm_lowTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[1] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[1]*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[1]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[1]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[1]*x
# }
# ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowTmeanDecJanFeb)
# ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midTmeanDecJanFeb)
# ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highTmeanDecJanFeb)
# ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.ppt_normhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normhigh[2,], ci.low = ci.ppt_normhigh[1,], ci.high = ci.ppt_normhigh[3,], ci.group = "highTmeanDecJanFeb")
# ci.ppt_normmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normmid[2,], ci.low = ci.ppt_normmid[1,], ci.high = ci.ppt_normmid[3,], ci.group = "midTmeanDecJanFeb")
# ci.ppt_normlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normlow[2,], ci.low = ci.ppt_normlow[1,], ci.high = ci.ppt_normlow[3,], ci.group = "lowTmeanDecJanFeb")
# ppt_norm_Tmean_DecJanFebint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
# ggplot(data = ppt_norm_Tmean_DecJanFebint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# 
# #tmp_norm and Precip_DecJanFeb interaction
# tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
# tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Precip_DecJanFeb_range <- quantile(grow_train$Precip_DecJanFeb, c(0.2, 0.8))
# growthpredictiontnorm_highPrecipDecJanFeb <- growthpredictiontnorm_lowPrecipDecJanFeb <- growthpredictiontnorm_midPrecipDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
#   growthpredictiontnorm_highPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[2]*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[2]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[2]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[2]*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[2]*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[2]*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictiontnorm_midPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictiontnorm_lowPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[1]*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[1]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[1]*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[1]*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[1]*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[1]*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowPrecipDecJanFeb)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midPrecipDecJanFeb)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highPrecipDecJanFeb)
# ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highPrecipDecJanFeb")
# ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midPrecipDecJanFeb")
# ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowPrecipDecJanFeb")
# tmp_norm_Precip_DecJanFebint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = tmp_norm_Precip_DecJanFebint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# 
# #tmp_norm and Precip_JulAug interaction
# tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
# tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Precip_JulAug_range <- quantile(grow_train$Precip_JulAug, c(0.2, 0.8))
# growthpredictiontnorm_highPrecipJulAug <- growthpredictiontnorm_lowPrecipJulAug <- growthpredictiontnorm_midPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
#   growthpredictiontnorm_highPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug_range[2]*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug_range[2]*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug_range[2] +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug_range[2] + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug_range[2] + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictiontnorm_midPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictiontnorm_lowPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug_range[1]*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug_range[1]*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug_range[1] +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug_range[1] + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug_range[1] + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowPrecipJulAug)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midPrecipJulAug)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highPrecipJulAug)
# ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highPrecipJulAug")
# ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midPrecipJulAug")
# ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowPrecipJulAug")
# tmp_norm_Precip_JulAugint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = tmp_norm_Precip_JulAugint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# 
# #tmp_norm and Tmean_JulAug interaction
# tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
# tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Tmean_JulAug_range <- quantile(grow_train$Tmean_JulAug, c(0.2, 0.8))
# growthpredictiontnorm_highTmeanJulAug <- growthpredictiontnorm_lowTmeanJulAug <- growthpredictiontnorm_midTmeanJulAug <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
#   growthpredictiontnorm_highTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug_range[2] + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug_range[2] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug_range[2]*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug_range[2]*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug_range[2]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug_range[2]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug_range[2]*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictiontnorm_midTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictiontnorm_lowTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug_range[1] + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug_range[1] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug_range[1]*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug_range[1]*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug_range[1]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug_range[1]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug_range[1]*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowTmeanJulAug)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midTmeanJulAug)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highTmeanJulAug)
# ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highTmeanJulAug")
# ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midTmeanJulAug")
# ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowTmeanJulAug")
# tmp_norm_Tmean_JulAugint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = tmp_norm_Tmean_JulAugint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# 
# #tmp_norm and Tmean_DecJanFeb interaction
# tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
# tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# Precip_JulAug <- mean(grow_train$Precip_JulAug)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Tmean_DecJanFeb_range <- quantile(grow_train$Tmean_JulAug, c(0.2, 0.8))
# growthpredictiontnorm_highTmeanDecJanFeb <- growthpredictiontnorm_lowTmeanDecJanFeb <- growthpredictiontnorm_midTmeanDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
#   growthpredictiontnorm_highTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[2] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[2]*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[2]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[2]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[2]*x
#   
#   growthpredictiontnorm_midTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictiontnorm_lowTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[1] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[1]*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[1]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[1]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[1]*x
# }
# tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowTmeanDecJanFeb)
# tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midTmeanDecJanFeb)
# tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highTmeanDecJanFeb)
# ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highTmeanDecJanFeb")
# ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midTmeanDecJanFeb")
# ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowTmeanDecJanFeb")
# tmp_norm_Tmean_DecJanFebint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = tmp_norm_Tmean_DecJanFebint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# 
# #Precip_JulAug and Tmean_DecJanFeb interaction
# Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
# Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
# x <- mean(grow_train$DIA_prev)
# tmp_norm <- mean(grow_train$tmp_norm)
# ppt_norm <- mean(grow_train$ppt_norm)
# Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
# Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
# Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
# Tmean_DecJanFeb_range <- quantile(grow_train$Tmean_JulAug, c(0.2, 0.8))
# growthpredictionPrecipJulAug_highTmeanDecJanFeb <- growthpredictionPrecipJulAug_lowTmeanDecJanFeb <- growthpredictionPrecipJulAug_midTmeanDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
# 
# for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
#   growthpredictionPrecipJulAug_highTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[2] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[2] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[2]*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[2]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[2]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[2]*x
#   
#   growthpredictionPrecipJulAug_midTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*x
#   
#   growthpredictionPrecipJulAug_lowTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
#     plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
#     plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
#     plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
#     plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x+
#     plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[1] + 
#     plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*x + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[1] +
#     plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*x + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[1]*ppt_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[1]*tmp_norm + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[1]*Precip_JulAug + 
#     plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[1]*x
# }
# PrecipJulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowTmeanDecJanFeb)
# PrecipJulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midTmeanDecJanFeb)
# PrecipJulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highTmeanDecJanFeb)
# ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
# ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
# ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highTmeanDecJanFeb")
# ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midTmeanDecJanFeb")
# ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowTmeanDecJanFeb")
# Precip_JulAug_Tmean_DecJanFebint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
# ggplot(data = Precip_JulAug_Tmean_DecJanFebint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
#   geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
# 
# 
# #ggplot map of tmp_norm
# all_states <- map_data("state")
# states <- subset(all_states, region %in% c("arizona"))
# coordinates(states)<-~long+lat
# class(states)
# proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
# mapdata<-states
# mapdata<-data.frame(mapdata)
# ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
#   geom_point(data = grow, aes(x = LON, y = LAT, color = tmp_norm))+ scale_color_gradient2(low = "blue", mid = "white", high = "red")+
#   theme_bw()
# 
# 
# #ggplot map of ppt_norm
# all_states <- map_data("state")
# states <- subset(all_states, region %in% c("arizona"))
# coordinates(states)<-~long+lat
# class(states)
# proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
# mapdata<-states
# mapdata<-data.frame(mapdata)
# ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
#   geom_point(data = grow, aes(x = LON, y = LAT, color = ppt_norm))+ scale_color_gradient2(low = "red", mid = "white", high = "blue")+
#   theme_bw()
# 
# 
# #ggplots of monsoon climate
# monsoonxtmpyr <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = tmp_yr)) + geom_point()
# monsoonxtmpyr
# 
# monsoonxtmpnorm <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = tmp_norm)) + geom_point()
# monsoonxtmpnorm
# 
# monsoonxpptnorm <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = ppt_norm)) + geom_point()
# monsoonxpptnorm
# 
# monsoonxsize <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = DIA_prev)) + geom_point()
# monsoonxsize
# 
# coolxmonsoon <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = Precip_DecJanFeb)) + geom_point()
# coolxmonsoon
# 
# pfallxmonsoon <- ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = Precip_SepOct)) + geom_point()
# pfallxmonsoon
# 
# fsummerxmonsoon <-ggplot(data = grow.monsoon, aes(x = Precip_JulAug, y = Precip_MarAprMay)) + geom_point()
# fsummerxmonsoon
# 
# 
# 
# 
# 
# line_type<-c(NA,"solid","dashed","dotted","dotdash")
# quartiles<-c(NA,"1st quartile","2nd quartile","3rd quartile","4th quartile")
# 
# q=quantile(climdata.scaled$ppt_norm)
# norms=c(0)
# for (i in 2:5){
#   
#   n=mean(c(q[i-1],q[i]))
#   norms[i]<-n
# }
# 
# climdata.scaled$ppt_norm_cat[climdata.scaled$ppt_norm<=q[2]]<-quartiles[2]
# climdata.scaled$ppt_norm_cat[climdata.scaled$ppt_norm>q[2] & climdata.scaled$ppt_norm<=q[3]]<-quartiles[3]
# climdata.scaled$ppt_norm_cat[climdata.scaled$ppt_norm>q[3] & climdata.scaled$ppt_norm<=q[4]]<-quartiles[4]
# climdata.scaled$ppt_norm_cat[climdata.scaled$ppt_norm>q[4] & climdata.scaled$ppt_norm<=q[5]]<-quartiles[5]
# 
# pptyr_plot <- ggplot(data = climdata.scaled, aes(x = ppt_yr, y = growth, col = ppt_norm_cat)) + 
#   geom_point()+
#   stat_function(fun=pfun_interall,args=list(temp=mean(climdata.scaled$tmp_yr),norm=norms[2]),
#                 aes(linetype=quartiles[2]),size=1,colour="black")+
#   stat_function(fun=pfun_interall,args=list(temp=mean(climdata.scaled$tmp_yr),norm=norms[3]),
#                 aes(linetype=quartiles[3]),size=1,colour="black")+
#   stat_function(fun=pfun_interall,args=list(temp=mean(climdata.scaled$tmp_yr),norm=norms[4]),
#                 aes(linetype=quartiles[4]),size=1,colour="black")+
#   stat_function(fun=pfun_interall,args=list(temp=mean(climdata.scaled$tmp_yr),norm=norms[5]),
#                 aes(linetype=quartiles[5]),size=1,colour="black")+
#   scale_colour_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
#                       values=c("1st quartile"="#41ae76","2nd quartile"="#238b45",
#                                "3rd quartile"="#006d2c","4th quartile"="#00441b"),
#                       labels=quartiles[2:5])+
#   scale_linetype_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
#                         values=c("1st quartile"="solid","2nd quartile"="dashed",
#                                  "3rd quartile"="dotdash","4th quartile"="dotted"),labels=quartiles[2:5])+
#   mytheme+labs(x="Annual precipitation",y="Growth")
# pptyr_plot  
# 
# tmpyr_plot <- ggplot(data = climdata.scaled, aes(x = tmp_yr, y = growth, col = ppt_norm_cat)) + 
#   geom_point()+
#   stat_function(fun=tfun_interall,args=list(ppt=mean(climdata.scaled$ppt_yr),norm=norms[2]),
#                 aes(linetype=quartiles[2]),size=0.75,colour="black")+
#   stat_function(fun=tfun_interall,args=list(ppt=mean(climdata.scaled$ppt_yr),norm=norms[3]),
#                 aes(linetype=quartiles[3]),size=0.75,colour="black")+
#   stat_function(fun=tfun_interall,args=list(ppt=mean(climdata.scaled$ppt_yr),norm=norms[4]),
#                 aes(linetype=quartiles[4]),size=0.75,colour="black")+
#   stat_function(fun=tfun_interall,args=list(ppt=mean(climdata.scaled$ppt_yr),norm=norms[5]),
#                 aes(linetype=quartiles[5]),size=0.75,colour="black")+
#   scale_colour_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
#                       values=c("1st quartile"="#41ae76","2nd quartile"="#238b45",
#                                "3rd quartile"="#006d2c","4th quartile"="#00441b"),
#                       labels=quartiles[2:5])+
#   scale_linetype_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
#                         values=c("1st quartile"="solid","2nd quartile"="dashed",
#                                  "3rd quartile"="dotdash","4th quartile"="dotted"),labels=quartiles[2:5])+
#   mytheme+labs(x="Annual temperature",y="Growth")
# tmpyr_plot  
# 
# ppt_violin <- ggplot(data = climdata.scaled, aes(x = ppt_norm_cat, y = ppt_yr)) + 
#   geom_violin()+
#   mytheme+labs(x="PPT norm",y="Annual precipitation")
# ppt_violin
# 
# ppt_hist <- ggplot(data = climdata.scaled, aes(x = ppt_yr,fill = ppt_norm_cat)) + 
#   geom_histogram()+
#   scale_fill_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
#                     values=c("1st quartile"="#41ae76","2nd quartile"="#238b45",
#                              "3rd quartile"="#006d2c","4th quartile"="#00441b"),
#                     labels=quartiles[2:5])+
#   mytheme+labs(x="Annual precipitation",y="Count")
# ppt_hist
# 
# ggsave("pptyr_plot.png",pptyr_plot)
# ggsave("tmpyr_plot.png",tmpyr_plot)
# ggsave("ppt_violin.png",ppt_violin)
# ggsave("ppt_hist.png",ppt_hist)
# 
# 
# 
# 
# write.csv(summary$summary,"./piedsummarylong.csv")
# 
# pied_grow_coef2<-summary$summary
# pied_grow_coef2
# save(pied_grow_coef2, grow, file="./pied_grow_coef2.rda")
# 
# save(pied_grow_coef, grow, file="./pied_grow_coef1.rda")
# 
# pial_post_grow<-fit_grow_df
# write.csv(pial_post_grow,"./Ogle/pial_post_grow.csv",row.names=F)
# 
# pial_grow_coef<-read.csv("pial_grow_ogle_coef1.csv",row.names=1)
# 
# #Check residuals
# pred<-grow_coef[1,2]+grow_coef[2,2]*grow$Size_t+grow_coef[28,2]*(grow$scovar%*%grow_coef[3:27,2])+grow_coef[29,2]*(grow$scovar%*%grow_coef[3:27,2])^2
# resid<-log(grow$Growth+0.001)*10-pial_grow_coef[83:1533,1]   #pred
# hist(resid)
# qqnorm(resid)
# qqline(resid)
# 
# #Plot of temperature coefficients
# mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                panel.background = element_blank(), axis.line = element_line(colour = "black"),
#                legend.text=element_text(size=8),axis.text=element_text(size=8),
#                axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),
#                axis.line.x = element_line(color="black", size = 0.3),
#                axis.line.y = element_line(color="black", size = 0.3))
# 
# ##Growth
# grow$Growth_1000<-grow$Size_t1-grow$Size_t
# q<-quantile(grow$Size_t,seq(0,1,length=4))
# Size=(q[1:3]+q[2:4])/2
# Dens=median(grow$Density)
# grow1<-subset(grow,grow$Size_t<q[2])
# grow2<-subset(grow,grow$Size_t>=q[2] & grow$Size_t<q[3])
# grow3<-subset(grow,grow$Size_t>q[3])
# 
# spei_min<-apply(grow$scovar,2,min)
# spei_max<-apply(grow$scovar,2,max)
# spei<-matrix(NA,50,25)
# for(i in 1:25){
#   spei[,i]<-seq(spei_min[i],spei_max[i],length=50)
# }
# 
# ncuts=25
# 
# Temperature1<-as.matrix(grow1$scovar*pial_grow_coef[9:33,1])
# Temperature1<-as.data.frame(.rowSums(Temperature1,nrow(Temperature1),25))
# Temperature1_2<-as.matrix((grow1$scovar^2)*pial_grow_coef[9:33,1])
# Temperature1_2<-as.data.frame(.rowSums(Temperature1_2,nrow(Temperature1),25))
# 
# names(Temperature1)<-"Temps"
# names(Temperature1_2)<-"Temps"
# chopsize1<-cut(Temperature1$Temps,ncuts)
# growbinned1<-as.vector(sapply(split(grow1$Growth_1000,chopsize1),mean))
# tempbinned1<-as.vector(sapply(split(Temperature1$Temps,chopsize1),mean))
# tempbinned1_2<-as.vector(sapply(split(Temperature1_2$Temps,chopsize1),mean))
# 
# g_binned1<-as.data.frame(cbind(tempbinned1,tempbinned1_2,growbinned1))
# names(g_binned1)<-c("temperature","temperature2","grow")
# g_binned1$size<-Size[1]
# 
# grow1<-as.data.frame(cbind(Temperature1,Temperature1_2,grow1$Growth_1000))
# names(grow1)<-c("temperature","temperature2","grow")
# grow1$size<-Size[1]
# 
# Temperature2<-as.matrix(grow2$scovar*pial_grow_coef[9:33,1])
# Temperature2<-as.data.frame(.rowSums(Temperature2,nrow(Temperature2),25))
# Temperature2_2<-as.matrix((grow2$scovar^2)*pial_grow_coef[9:33,1])
# Temperature2_2<-as.data.frame(.rowSums(Temperature2_2,nrow(Temperature2),25))
# 
# names(Temperature2)<-"Temps"
# names(Temperature2_2)<-"Temps"
# chopsize2<-cut(Temperature2$Temps,ncuts)
# growbinned2<-as.vector(sapply(split(grow2$Growth_1000,chopsize2),mean))
# tempbinned2<-as.vector(sapply(split(Temperature2$Temps,chopsize2),mean))
# tempbinned2_2<-as.vector(sapply(split(Temperature2_2$Temps,chopsize2),mean))
# 
# g_binned2<-as.data.frame(cbind(tempbinned2,tempbinned2_2,growbinned2))
# names(g_binned2)<-c("temperature","temperature2","grow")
# g_binned2$size<-Size[2]
# 
# grow2<-as.data.frame(cbind(Temperature2,Temperature2_2,grow2$Growth_1000))
# names(grow2)<-c("temperature","temperature2","grow")
# grow2$size<-Size[2]
# 
# Temperature3<-as.matrix(grow3$scovar*pial_grow_coef[9:33,1])
# Temperature3<-as.data.frame(.rowSums(Temperature3,nrow(Temperature3),25))
# Temperature3_2<-as.matrix((grow3$scovar^2)*pial_grow_coef[9:33,1])
# Temperature3_2<-as.data.frame(.rowSums(Temperature3_2,nrow(Temperature3),25))
# 
# names(Temperature3)<-"Temps"
# names(Temperature3_2)<-"Temps"
# chopsize3<-cut(Temperature3$Temps,ncuts)
# growbinned3<-as.vector(sapply(split(grow3$Growth_1000,chopsize3),mean))
# tempbinned3<-as.vector(sapply(split(Temperature3$Temps,chopsize3),mean))
# tempbinned3_2<-as.vector(sapply(split(Temperature3_2$Temps,chopsize3),mean))
# 
# g_binned3<-as.data.frame(cbind(tempbinned3,tempbinned3_2,growbinned3))
# names(g_binned3)<-c("temperature","temperature2","grow")
# g_binned3$size<-Size[3]
# 
# grow3<-as.data.frame(cbind(Temperature3,Temperature3_2,grow3$Growth_1000))
# names(grow3)<-c("temperature","temperature2","grow")
# grow3$size<-Size[3]
# 
# grow_plot_data<-rbind(grow1,grow2,grow3)
# grow_binned_plot_data<-rbind(g_binned1,g_binned2,g_binned3)
# 
# g_fun<-function(spei,size){
#   spei1<-(spei/100)%*%pial_grow_coef[9:33,1]
#   spei2<-((spei/100)^2)%*%pial_grow_coef[9:33,1]
#   size=size
#   grow<-(pial_grow_coef[1,1]+pial_grow_coef[8,1]*size+
#            pial_grow_coef[34,1]*spei1+pial_grow_coef[35,1]*spei2+
#            pial_grow_coef[36,1]*spei1*size+pial_grow_coef[37,1]*spei2*size)-size
#   return(data.frame(spei=spei1*100,grow=grow,size=size,row.names=NULL))
# }
# 
# Growth1<-g_fun(spei=spei,size=Size[1])
# Growth2<-g_fun(spei=spei,size=Size[2])
# Growth3<-g_fun(spei=spei,size=Size[3])
# 
# Growth<-rbind(Growth1,Growth2,Growth3)
# 
# grow_plot<-ggplot(grow_binned_plot_data,aes(x=temperature,y=grow,col=size))+
#   labs(x="SPEI",y="Growth (cm)")+
#   geom_point()+
#   geom_line(data=Growth,aes(x=spei,y=grow,group=size))+
#   scale_colour_viridis("Size")+
#   coord_cartesian(xlim=c(min(grow_plot_data$temperature,na.rm=T),max(grow_plot_data$temperature,na.rm=T)))+
#   #ylim=c((min(grow_plot_data$size,na.rm=T))+5,(max(grow_plot_data$size,na.rm=T)+1)))+
#   mytheme
# 
# png(file="./Output/grow_spei.png",400,360,type="cairo")
# grow_plot
# dev.off()
