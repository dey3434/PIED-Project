### Stan models for PIAL growth using technique in Ogle et al., Ecology Letters to test for climate effects
## Sharmila Dey
# 22 June 2020
# setwd("/home/rstudio")
#load(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_grow_coef2.rda"))
#fit_grow <- readRDS(url("https://de.cyverse.org/dl/d/3C958A2C-6439-47C9-AFB9-2EC4A4C3824A/ppt_tmp_springfall_sizefix.RDS"))
#fit_grow <- readRDS("ppt_tmp_springfall_sizefix.RDS")
#fit_grow_old <- readRDS("data/ppt_tmp_springfall.RDS")
#install.packages("rstan", version = "2.19.3", repos = "http://cran.us.r-project.org")
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

grow.monsoon<-na.omit(grow.new) %>% 
  mutate_at(scale, .vars = vars(Precip_JulAug, Precip_NovDecJanFebMar, Tmean_AprMayJun, Tmean_SepOct)) %>%
  arrange(PLOT,SUBP,name) %>%
  mutate(PlotCD=as.numeric(factor(ST_PLT, levels = unique(ST_PLT))),treeCD=as.numeric(factor(name,levels=unique(name))),
         growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2))

set.seed(2023)
split=0.20
trainIndex <- createDataPartition(grow.monsoon$name, p=split, list=FALSE)
grow_test <- grow.monsoon[trainIndex,]
grow_train <- grow.monsoon[-trainIndex,]



xG<-as.matrix(cbind(grow_train$Precip_JulAug, grow_train$Precip_NovDecJanFebMar, 
                    grow_train$Tmean_AprMayJun, grow_train$Tmean_SepOct, grow_train$DIA_prev,
                    grow_train$DIA_prev*grow_train$Precip_JulAug,
                    grow_train$DIA_prev*grow_train$Precip_NovDecJanFebMar, grow_train$DIA_prev*grow_train$Tmean_AprMayJun,
                    grow_train$DIA_prev*grow_train$Tmean_SepOct, grow_train$Precip_JulAug*grow_train$Precip_NovDecJanFebMar,
                    grow_train$Precip_JulAug*grow_train$Tmean_AprMayJun, grow_train$Precip_JulAug*grow_train$Tmean_SepOct,
                    grow_train$Precip_NovDecJanFebMar*grow_train$Tmean_AprMayJun, grow_train$Precip_NovDecJanFebMar*grow_train$Tmean_SepOct,
                    grow_train$Tmean_AprMayJun*grow_train$Tmean_SepOct))
xGtest<-as.matrix(cbind(grow_test$Precip_JulAug, grow_test$Precip_NovDecJanFebMar, 
                        grow_test$Tmean_AprMayJun, grow_test$Tmean_SepOct, grow_test$DIA_prev,
                        grow_test$DIA_prev*grow_test$Precip_JulAug,
                        grow_test$DIA_prev*grow_test$Precip_NovDecJanFebMar, grow_test$DIA_prev*grow_test$Tmean_AprMayJun,
                        grow_test$DIA_prev*grow_test$Tmean_SepOct, grow_test$Precip_JulAug*grow_test$Precip_NovDecJanFebMar,
                        grow_test$Precip_JulAug*grow_test$Tmean_AprMayJun, grow_test$Precip_JulAug*grow_test$Tmean_SepOct,
                        grow_test$Precip_NovDecJanFebMar*grow_test$Tmean_AprMayJun, grow_test$Precip_NovDecJanFebMar*grow_test$Tmean_SepOct,
                        grow_test$Tmean_AprMayJun*grow_test$Tmean_SepOct))
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


sink("model_8.stan")
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



csvfiles <- here::here("results", paste0("ppt_tmp_springfall_sizefix_", 1:3, ".csv"))

if (all(file.exists(csvfiles))) {
  fit_grow <- read_stan_csv(csvfiles, col_major = TRUE) 
} else {
  fit_grow <- stan(file = 'model_8.stan', data = pied_dat, 
                   iter = 5000,
                   warmup = 1000,
                   chains = 3, cores = 8, 
                   sample_file = here::here("results", "ppt_tmp_springfall_sizefix"))
}


# chain1 <- rstan::read_stan_csv("/home/rstudio/ppt_tmp_springfall_sizefix_1.csv")
# chain2 <- rstan::read_stan_csv("/home/rstudio/ppt_tmp_springfall_sizefix_2.csv")
# chain3 <- rstan::read_stan_csv("/home/rstudio/ppt_tmp_springfall_sizefix_3.csv")

# fit_grow <- rstan::read_stan_csv(csvfiles)

# not sure the name is correct here
# saveRDS(fit_grow, file = "ppt_tmp_springfall_sizefix.RDS")
summary<-summary(fit_grow)
summary

plotdata<-select(as.data.frame(fit_grow),"yrep[1]":"yrep[8780]")
plotdatainterval<-select(as.data.frame(fit_grow), "u_beta[1]":"u_beta[28]")
colnames(plotdatainterval) <- c("u_beta_Precip_JulAug", "u_beta_Precip_NovDecJanFebMar", 
                                "u_beta_Tmean_AprMayJun", "u_beta_Tmean_SepOct", "u_beta_DIA_prev",
                                "u_beta_DIA_prev_Precip_JulAug",
                                "u_beta_DIA_prev_Precip_NovDecJanFebMar", "u_beta_DIA_prev_Tmean_AprMayJun",
                                "u_beta_DIA_prev_Tmean_SepOct", "u_beta_Precip_JulAug_Precip_NovDecJanFebMar",
                                "u_beta_Precip_JulAug_Tmean_AprMayJun", "u_beta_Precip_JulAug_Tmean_SepOct",
                                "u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun", "u_beta_Precip_NovDecJanFebMar_Tmean_SepOct",
                                "u_beta_Tmean_AprMayJun_Tmean_SepOct")

# get summaries of plotdatainterval:

df <- reshape2::melt(plotdatainterval)
beta.summaries <- df %>% group_by(variable) %>%
  summarise(mean = mean(value),
            ci.lo = quantile(value, 0.025),
            ci.hi = quantile(value, 0.975))

beta.summaries$allpos <- ifelse(beta.summaries$mean > 0 &  beta.summaries$ci.lo > 0 &  beta.summaries$ci.hi > 0, "yes", "no")
beta.summaries$allneg <- ifelse(beta.summaries$mean <= 0 &  beta.summaries$ci.lo <= 0 &  beta.summaries$ci.hi <= 0, "yes", "no")
beta.summaries$significant <- ifelse(beta.summaries$allpos == "yes" | beta.summaries$allneg == "yes", "significant", "not significant")
write.csv(beta.summaries, "betasummaries_model8_threechain_PIED.csv", row.names = FALSE)

####
#Validation-- This will be in a separate File
ext_fit <- rstan::extract(fit_grow)
yrep <- ext_fit$yrep
#yrep <- exp(yrep)
mean.pred <- apply(ext_fit$yrep, 2, mean)
p.o.df <- data.frame(predicted = exp(mean.pred), observed = exp(grow_test$loggrowth), error = (exp(mean.pred) - exp(grow_test$loggrowth)))
meansqrd <- (mean(p.o.df$error))^2

##### 
save(p.o.df, file = here::here("results", "model-8-pred-obs.RData"))

# ggplot(p.o.df, aes(predicted, observed)) + geom_point(alpha = 0.1) + geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
#   ylim(0, 10) + xlim(0,10)

####
p_pred_vs_observed <- ggplot(p.o.df, aes(predicted, observed)) + 
  geom_point(alpha = 0.1) + 
  geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
  ylim(0, 10) + xlim(0,10)
p_pred_vs_observed
ggsave(here::here("images", "model_8", "pred_vs_observed.png"), p_pred_vs_observed)



sigma <- as.data.frame(fit_grow)[,"sigma_y"]
mu <- as.matrix(plotdatainterval) %*% t(xG)
ll <- matrix(0, length(sigma), length(yG))
for(i in 1:length(sigma)){
  ll[i,] <- dnorm(yG, mu[i,], sd = sigma[i], log = TRUE)
}
newll <- as.matrix(ll)
r_eff <- relative_eff(exp(ll), chain_id = rep(1:3, each = 4000), cores = 1)
leaveoneout <- loo::loo(as.matrix(ll), r_eff = r_eff, save_psis = TRUE, cores = 1)

save(ll, r_eff, leaveoneout, file = here::here("results", "model-8-loo.RData"))

yrep <- matrix(0, length(sigma), length(yG))
for(i in 1:length(sigma)){
  yrep[i,] <- rnorm(yG, mu[i,], sd = sigma[i])
}
psis <- leaveoneout$psis_object
keep_obs <- sample(1:length(yG), 100)
lw <- weights(psis)
ppc1 <- ppc_loo_intervals(yG, yrep = yrep, psis_object = psis, subset = keep_obs, order = "median") 
ppc2 <- ppc_loo_pit_overlay(yG, yrep = yrep, lw = lw)
ppc3 <- ppc_loo_pit_qq(yG, yrep = yrep, lw = lw)

ggsave(here::here("images", "model_8", "ppc-plot-1.png"),
       ppc1, width = 16/3, height = 9)
ggsave(here::here("images", "model_8", "ppc-plot-2.png"),
       ppc2, width = 16/3, height = 9)
ggsave(here::here("images", "model_8", "ppc-plot-3.png"),
       ppc3, width = 16/3, height = 9)

ppc_dens_overlay(yGtest, as.matrix(plotdata))

ext_fit <- rstan::extract(fit_grow)
yrep <- ext_fit$yrep
#yrep <- exp(yrep)
mean.pred <- apply(ext_fit$yrep, 2, median)
p.o.df <- data.frame(predicted = exp(mean.pred), observed = exp(grow_test$loggrowth), error = (exp(mean.pred) - exp(grow_test$loggrowth)))
meansqrd <- (mean(p.o.df$error))^2
ggplot(p.o.df, aes(predicted, observed)) + geom_point(alpha = 0.1) + geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
  ylim(0, 10) + xlim(0,10)


#Validation
sigma <- as.data.frame(fit_grow)[,"sigma_y"]
mu <- as.matrix(plotdatainterval) %*% t(xG)
ll <- matrix(0, length(sigma), length(yG))
for(i in 1:length(sigma)){
  ll[i,] <- dnorm(yG, mu[i,], sd = sigma[i], log = TRUE)
}
newll <- as.matrix(ll)
r_eff <- relative_eff(exp(ll), chain_id = rep(1:3, each = 4000), cores = 8)
leaveoneout <- loo(as.matrix(ll), r_eff = r_eff, save_psis = TRUE, cores = 8)

yrep <- matrix(0, length(sigma), length(yG))
for(i in 1:length(sigma)){
  yrep[i,] <- rnorm(yG, mu[i,], sd = sigma[i])
}
psis <- leaveoneout$psis_object
keep_obs <- sample(1:length(yG), 100)
lw <- weights(psis)
ppc_loo_intervals(yG, yrep = yrep, psis_object = psis, subset = keep_obs, order = "median") 
ppc_loo_pit_overlay(yG, yrep = yrep, lw = lw)
ppc_loo_pit_qq(yG, yrep = yrep, lw = lw)


yrep_mean <- apply(yrep, MARGIN = 2, FUN = mean)
yrep_ci.low <- apply(yrep, MARGIN = 2, FUN = function(x){quantile(x, 0.025)})
yrep_ci.high <- apply(yrep, MARGIN = 2, FUN = function(x){quantile(x, 0.975)})

p.o.df <- data.frame(ci.low = yrep_ci.low, mean = yrep_mean, ci.high = yrep_ci.high, observed = yG,
                     ppt_norm = grow_train$ppt_norm, tmp_norm = grow_train$tmp_norm, Precip_JulAug = grow_train$Precip_JulAug, 
                     Precip_NovDecJanFebMar = grow_train$Precip_NovDecJanFebMar, Tmean_AprMayJun = grow_train$Tmean_AprMayJun, 
                     Tmean_SepOct = grow_train$Tmean_SepOct, DIA_prev = grow_train$DIA_prev, ELEV = grow_train$ELEV, 
                     SLOPE = grow_train$SLOPE, ASPECT = grow_train$ASPECT, tmp_yr = grow_train$tmp_yr, ppt_yr = grow_train$ppt_yr, 
                     year = grow_train$year, PlotCD = grow_train$PlotCD, treeCD = grow_train$treeCD,
                     name = grow_train$name)
p.o.df$overpredicted <- ifelse(p.o.df$observed <= p.o.df$ci.low, "over predicted", "within confidence interval")
overpredictedpoints <- p.o.df %>% filter(observed <= ci.low)

ggplot(data = p.o.df, aes(x = ELEV, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = tmp_norm, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = ppt_norm, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = Precip_JulAug, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = Precip_NovDecJanFebMar, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = Tmean_AprMayJun, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = Tmean_SepOct, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = DIA_prev, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = SLOPE, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = ASPECT, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = tmp_yr, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = ppt_yr, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = year, y = observed, color = overpredicted)) + geom_point(size = 0.5)
ggplot(data = p.o.df, aes(x = name, y = observed, color = overpredicted)) + geom_point(size = 0.5)



## Subset posterior predictive plot by size
size_q<-quantile(grow$DIA_prev)
sizeq1<-which(grow_test$DIA_prev<=size_q[2])
sizeq2<-which(grow_test$DIA_prev>size_q[2] & grow_test$DIA_prev<=size_q[3])
sizeq3<-which(grow_test$DIA_prev>size_q[3] &grow_test$DIA_prev<=size_q[4])
sizeq4<-which(grow_test$DIA_prev>size_q[4])

ppc_dens_overlay(yGtest[sizeq1], as.matrix(plotdata)[,sizeq1])
ppc_dens_overlay(yGtest[sizeq2], as.matrix(plotdata)[,sizeq2])
ppc_dens_overlay(yGtest[sizeq3], as.matrix(plotdata)[,sizeq3])
ppc_dens_overlay(yGtest[sizeq4], as.matrix(plotdata)[,sizeq4])

## Subset posterior predicitve plot by ppt_norm
ppt_norm_q<-quantile(grow$ppt_norm)
ppt_normq1<-which(grow_test$ppt_norm<=ppt_norm_q[2])
ppt_normq2<-which(grow_test$ppt_norm>ppt_norm_q[2] & grow_test$ppt_norm<=ppt_norm_q[3])
ppt_normq3<-which(grow_test$ppt_norm>ppt_norm_q[3] &grow_test$ppt_norm<=ppt_norm_q[4])
ppt_normq4<-which(grow_test$ppt_norm>ppt_norm_q[4])

ppc_dens_overlay(yGtest[ppt_normq1], as.matrix(plotdata)[,ppt_normq1])
ppc_dens_overlay(yGtest[ppt_normq2], as.matrix(plotdata)[,ppt_normq2])
ppc_dens_overlay(yGtest[ppt_normq3], as.matrix(plotdata)[,ppt_normq3])
ppc_dens_overlay(yGtest[ppt_normq4], as.matrix(plotdata)[,ppt_normq4])

## Subset posterior predicitve plot by ppt_norm
tmp_norm_q<-quantile(grow$tmp_norm)
tmp_normq1<-which(grow_test$tmp_norm<=tmp_norm_q[2])
tmp_normq2<-which(grow_test$tmp_norm>tmp_norm_q[2] & grow_test$tmp_norm<=tmp_norm_q[3])
tmp_normq3<-which(grow_test$tmp_norm>tmp_norm_q[3] &grow_test$tmp_norm<=tmp_norm_q[4])
tmp_normq4<-which(grow_test$tmp_norm>tmp_norm_q[4])

ppc_dens_overlay(yGtest[tmp_normq1], as.matrix(plotdata)[,tmp_normq1])
ppc_dens_overlay(yGtest[tmp_normq2], as.matrix(plotdata)[,tmp_normq2])
ppc_dens_overlay(yGtest[tmp_normq3], as.matrix(plotdata)[,tmp_normq3])
ppc_dens_overlay(yGtest[tmp_normq4], as.matrix(plotdata)[,tmp_normq4])


##generating a plotting function
##ppcdensoverlay
make_plot <- function() {
  for (i in min(grow_test$year):max(grow_test$year)) {
    year<-which(grow_test$year == i)
    p = ppc_dens_overlay(yGtest[year], as.matrix(plotdata)[,year]) + 
      theme(
        plot.title = element_text(size = rel(2.5), legend.text = element_text(size = 16), 
                                  axis.text.x = element_text(size = 12),
                                  legend.key.size = unit(1.2, "lines")
        ) + xlim(-6.91, 3.96) +
          ggtitle(
            paste(i)
          ))
    print(p)
  }
}


if (!file.exists(here::here("images", "ppc_year-animation.gif"))) {
  
  gifski::save_gif(
    make_plot(),
    gif_file = here::here("images", "ppc_year-animation.gif"), 
    progress = FALSE,
    delay = 0.5, 
    height = 360, width = 640, units = "px"
  )
}

##density function for climate
make_ppt_plot <- function() {
  for (i in min(grow_test$year):max(grow_test$year)) {
    year <- i
    p = ggplot(grow_train[grow_train$year == year,], aes(x = ppt_yr )) + geom_density() +
      theme(
        plot.title = element_text(size = rel(2.5)),  legend.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.key.size = unit(1.2, "lines")
      ) + xlim(min(grow_train$ppt_yr), max(grow_train$ppt_yr)) +
      ggtitle(
        paste(i)
      )
    print(p)
  }
}

if (!file.exists(here::here("images", "ppt_year-animation.gif"))) {
  
  gifski::save_gif(
    make_ppt_plot(),
    gif_file = here::here("images", "ppt_year-animation.gif"), 
    progress = FALSE,
    delay = 0.5, 
    height = 360, width = 640, units = "px"
  )
}        


make_ppt_plot <- function() {
  for (i in min(grow_test$year):max(grow_test$year)) {
    year <- i
    p = ggplot(grow_train[grow_train$year == year,], aes(x = tmp_yr )) + geom_density() +
      theme(
        plot.title = element_text(size = rel(2.5)),  legend.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.key.size = unit(1.2, "lines")
      ) + xlim(min(grow_train$tmp_yr), max(grow_train$tmp_yr)) +
      ggtitle(
        paste(i)
      )
    print(p)
  }
}

if (!file.exists(here::here("images", "tmp_year-animation.gif"))) {
  
  gifski::save_gif(
    make_ppt_plot(),
    gif_file = here::here("images", "tmp_year-animation.gif"), 
    progress = FALSE,
    delay = 0.5, 
    height = 360, width = 640, units = "px"
  )
}        

#MCMC Intervals Plots
mcmc_intervals(plotdatainterval, prob = 0.5)
pdf("plotdatainterval_mcmc_intervals.pdf", height = 6, width = 10) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
mcmc_intervals(plotdatainterval, prob = 0.5) # put the code that makes one of the plots in here
dev.off() # "device off" tells R to stop printing stuff to the pdf


#MCMC Area Plot
pdf("plotdatainterval_mcmc_areas.pdf", height = 8, width = 20) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
color_scheme_set("purple")
mcmc_areas(plotdatainterval, prob = 0.8)
dev.off() # "device off" tells R to stop printing stuff to the pdf

#MCMC Traces
pdf("plotdatainterval_mcmc_traces.pdf", height = 6, width = 30) # tells R to save the following plots to a pdf named "filename.pdf" that is 6 inches wide and 6 inches width
color_scheme_set("mix-blue-red")
mcmc_trace(plotdatainterval,
           facet_args = list(nrow = 2))
#pars = c("alpha", "sigma")
dev.off() # "device off" tells R to stop printing stuff to the pdf

mcmcplot(As.mcmc.list(fit_grow))


png("plotdatainterval.png", height = 20, width = 24, units = "in", res = 200)
pairs.panels(as.matrix(plotdatainterval))
dev.off()



#Effects Plots
mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               legend.text=element_text(size=11),legend.title=element_text(size=12),
               legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
               axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
               axis.line.x = element_line(color="black", size = 0.3),
               axis.line.y = element_line(color="black", size = 0.3))

#effect of ppt_norm
ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for ppt_normrng
ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.25)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
growthpredictionpptnorm <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm))
pfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
  for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
    growthpredictionpptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  growthpredictionpptnorm
}
ppt_norm_prediction <- pfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, 
                                             Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
ppt_norm_prediction_tr <- exp(ppt_norm_prediction)
ci.ppt_norm <- apply(ppt_norm_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.ppt_norm.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_norm[2,], ci.low = ci.ppt_norm[1,], ci.high = ci.ppt_norm[3,])
ggplot() + 
  geom_ribbon(data = ci.ppt_norm.df, aes(x = ppt_norm, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.ppt_norm.df, aes(x = ppt_norm, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "tmp_norm")]), aes(x = ppt_norm))


#effect of tmp_norm
tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.35)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
growthpredictiontmpnorm <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm))
tfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
  for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
    growthpredictiontmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  growthpredictiontmpnorm
}
tmp_norm_prediction <- tfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, 
                                             Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
tmp_norm_prediction_tr <- exp(tmp_norm_prediction)
ci.tmp_norm <- apply(tmp_norm_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_norm.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_norm[2,], ci.low = ci.tmp_norm[1,], ci.high = ci.tmp_norm[3,])
ggplot() + 
  geom_ribbon(data = ci.tmp_norm.df, aes(x = tmp_norm, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.tmp_norm.df, aes(x = tmp_norm, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "tmp_norm")]), aes(x = tmp_norm))


#effect of Precip_JulAug
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
growthpredictionPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
pjafun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  growthpredictionPrecipJulAug
}
PrecipJulAug_prediction <- pjafun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, 
                                                   Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
PrecipJulAug_prediction_tr <- exp(PrecipJulAug_prediction)
ci.Precip_JulAug <- apply(PrecipJulAug_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAug.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAug[2,], ci.low = ci.Precip_JulAug[1,], ci.high = ci.Precip_JulAug[3,])
ggplot() + 
  geom_ribbon(data = ci.Precip_JulAug.df, aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.Precip_JulAug.df, aes(x = Precip_JulAug, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "Precip_JulAug")]), aes(x = Precip_JulAug))


#effect of Precip_NovDecJanFebMar
Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
growthpredictionPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar))
pndjfmfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  growthpredictionPrecipNovDecJanFebMar
}
PrecipNovDecJanFebMar_prediction <- pndjfmfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, 
                                                               Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
PrecipNovDecJanFebMar_prediction_tr <- exp(PrecipNovDecJanFebMar_prediction)
ci.Precip_NovDecJanFebMar <- apply(PrecipNovDecJanFebMar_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_NovDecJanFebMar.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_NovDecJanFebMar[3,])
ggplot() + 
  geom_ribbon(data = ci.Precip_NovDecJanFebMar.df, aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.Precip_NovDecJanFebMar.df, aes(x = Precip_NovDecJanFebMar, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "Precip_NovDecJanFebMar")]), aes(x = Precip_NovDecJanFebMar))

#effect of Tmean_AprMayJun
Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
growthpredictionTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun))
tamjfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  growthpredictionTmeanAprMayJun
}
TmeanAprMayJun_prediction <- tamjfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, 
                                                      Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
TmeanAprMayJun_prediction_tr <- exp(TmeanAprMayJun_prediction)
ci.Tmean_AprMayJun <- apply(TmeanAprMayJun_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Tmean_AprMayJun.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJun[2,], ci.low = ci.Tmean_AprMayJun[1,], ci.high = ci.Tmean_AprMayJun[3,])
ggplot() + 
  geom_ribbon(data = ci.Tmean_AprMayJun.df, aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.Tmean_AprMayJun.df, aes(x = Tmean_AprMayJun, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "Tmean_AprMayJun")]), aes(x = Tmean_AprMayJun))


#effect of Tmean_SepOct
Tmean_SepOctrng <- range(grow_train$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
growthpredictionTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct))
tsofun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,Tmean_AprMayJun,Tmean_SepOct){
  for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
    growthpredictionTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  growthpredictionTmeanSepOct
}
TmeanSepOct_prediction <- tsofun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, 
                                                  Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct)
TmeanSepOct_prediction_tr <- exp(TmeanSepOct_prediction)
ci.Tmean_SepOct <- apply(TmeanSepOct_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Tmean_SepOct.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_SepOct[2,], ci.low = ci.Tmean_SepOct[1,], ci.high = ci.Tmean_SepOct[3,])
ggplot() + 
  geom_ribbon(data = ci.Tmean_SepOct.df, aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.Tmean_SepOct.df, aes(x = Tmean_SepOct, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "Tmean_SepOct")]), aes(x = Tmean_SepOct))


#effect of DIA_prev
sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
size <- seq(sizerng[1], sizerng[2], by = 7.5)
Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
growthpredictionsize <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size))
sizefun_plotdatainterval<-function(size,tmp_norm,ppt_norm,Precip_JulAug,Precip_DecJanFeb,Tmean_JulAug,Tmean_DecJanFeb){
  for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
    growthpredictionsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
      plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
      plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
      plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
      plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
      plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
      plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
      plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
      plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
      plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
      plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
      plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
      plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
      plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  }
  growthpredictionsize
}
size_prediction <- sizefun_plotdatainterval(size = size, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_DecJanFeb = Precip_DecJanFeb, Tmean_JulAug = Tmean_JulAug, Tmean_DecJanFeb = Tmean_DecJanFeb)
size_prediction_tr <- exp(Tmean_DecJanFeb_prediction)
ci.size <- apply(size_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.size.df <- data.frame(size = size, median = ci.size[2,], ci.low = ci.size[1,], ci.high = ci.size[3,])
ggplot() + 
  geom_ribbon(data = ci.size.df, aes(x = size, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.size.df, aes(x = size, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "DIA_prev")]), aes(x = size, color = ppt_norm))


#Precip_JulAug and Tmean_AprMayJun
Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Precip_JulAug_range <- quantile(grow_train$Precip_JulAug, c(0.2, 0.8))
growthpredictionTmeanAprMayJun_highPrecipJulAug <- growthpredictionTmeanAprMayJun_lowPrecipJulAug <- growthpredictionTmeanAprMayJun_midPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
  growthpredictionTmeanAprMayJun_highPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug_range[2] +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug_range[2]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug_range[2]*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionTmeanAprMayJun_midPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionTmeanAprMayJun_lowPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug_range[1] +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug_range[1]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug_range[1]*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Tmean_AprMayJun_prediction_trlow <- exp(growthpredictionTmeanAprMayJun_lowPrecipJulAug)
Tmean_AprMayJun_prediction_trmid <- exp(growthpredictionTmeanAprMayJun_midPrecipJulAug)
Tmean_AprMayJun_prediction_trhigh <- exp(growthpredictionTmeanAprMayJun_highPrecipJulAug)
ci.Tmean_AprMayJunhigh <- apply(Tmean_AprMayJun_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Tmean_AprMayJunmid <- apply(Tmean_AprMayJun_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunlow <- apply(Tmean_AprMayJun_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highPrecip_JulAug")
ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midPrecip_JulAug")
ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowPrecip_JulAug")
Tmean_AprMayJun_Precip_JulAugint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
ggplot(data = Tmean_AprMayJun_Precip_JulAugint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = Tmean_AprMayJun_Precip_JulAugint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)

#Precip_JulAug and Tmean_SepOct
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Tmean_SepOct_range <- quantile(grow_train$Tmean_SepOct, c(0.2, 0.8))
growthpredictionPrecipJulAug_highTmeanSepOct <- growthpredictionPrecipJulAug_lowTmeanSepOct <- growthpredictionPrecipJulAug_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[2]
  
  growthpredictionPrecipJulAug_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[1]
}
Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowTmeanSepOct)
Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midTmeanSepOct)
Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highTmeanSepOct)
ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highTmean_SepOct")
ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midTmean_SepOct")
ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowTmean_SepOct")
Precip_JulAug_Tmean_SepOctint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
ggplot(data = Precip_JulAug_Tmean_SepOctint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = Precip_JulAug_Tmean_SepOctint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
  scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
  geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Precip_JulAug and Precip_NovDecJanFebMar
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Precip_NovDecJanFebMar_range <- quantile(grow_train$Precip_NovDecJanFebMar, c(0.2, 0.8))
growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar <- growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar <- growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[2] +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar_range[2] + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar_range[2] + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar_range[2] + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar_range[2]*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[1] +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar_range[1] + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar_range[1] + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar_range[1] + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar_range[1]*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Precip_JulAug_predictionpndjfm_trlow <- exp(growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar)
Precip_JulAug_predictionpndjfm_trmid <- exp(growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar)
Precip_JulAug_predictionpndjfm_trhigh <- exp(growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar)
ci.Precip_JulAughigh_Precip_NovDecJanFebMar <- apply(Precip_JulAug_predictionpndjfm_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugmid_Precip_NovDecJanFebMar <- apply(Precip_JulAug_predictionpndjfm_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAuglow_Precip_NovDecJanFebMar <- apply(Precip_JulAug_predictionpndjfm_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAughigh_Precip_NovDecJanFebMar.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh_Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_JulAughigh_Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_JulAughigh_Precip_NovDecJanFebMar[3,], ci.group = "highPrecip_NovDecJanFebMar")
ci.Precip_JulAugmid_Precip_NovDecJanFebMar.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid_Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_JulAugmid_Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_JulAugmid_Precip_NovDecJanFebMar[3,], ci.group = "midPrecip_NovDecJanFebMar")
ci.Precip_JulAuglow_Precip_NovDecJanFebMar.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow_Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_JulAuglow_Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_JulAuglow_Precip_NovDecJanFebMar[3,], ci.group = "lowPrecip_NovDecJanFebMar")
Precip_JulAug_Precip_NovDecJanFebMarint <- rbind(ci.Precip_JulAughigh_Precip_NovDecJanFebMar.df, ci.Precip_JulAugmid_Precip_NovDecJanFebMar.df, ci.Precip_JulAuglow_Precip_NovDecJanFebMar.df)
ggplot(data = Precip_JulAug_Precip_NovDecJanFebMarint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)

#Precip_JulAug and tmp_norm
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictionPrecipJulAug_hightnorm <- growthpredictionPrecipJulAug_lowtnorm <- growthpredictionPrecipJulAug_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowtnorm)
Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midtnorm)
Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_hightnorm)
ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "hightnorm")
ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midtnorm")
ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowtnorm")
Precip_JulAug_tnormint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
ggplot(data = Precip_JulAug_tnormint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Precip_JulAug and ppt_norm
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
growthpredictionPrecipJulAug_highpnorm <- growthpredictionPrecipJulAug_lowpnorm <- growthpredictionPrecipJulAug_midpnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_highpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[2]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[2]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_midpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_lowpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[1]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[1]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowpnorm)
Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midpnorm)
Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highpnorm)
ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highpnorm")
ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midpnorm")
ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowpnorm")
Precip_JulAug_pnormint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
ggplot(data = Precip_JulAug_pnormint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Precip_JulAug and size
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictionPrecipJulAug_highsize <- growthpredictionPrecipJulAug_lowsize <- growthpredictionPrecipJulAug_midsize <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipJulAug_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowsize)
Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midsize)
Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_highsize)
ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "highsize")
ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midsize")
ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowsize")
Precip_JulAug_sizeint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
ggplot(data = Precip_JulAug_sizeint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#ppt_norm and tmp_norm
ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictionpnorm_hightnorm <- growthpredictionpnorm_lowtnorm <- growthpredictionpnorm_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 

for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
  growthpredictionpnorm_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpnorm_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpnorm_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowtnorm)
ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midtnorm)
ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_hightnorm)
ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.ppt_normhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normhigh[2,], ci.low = ci.ppt_normhigh[1,], ci.high = ci.ppt_normhigh[3,], ci.group = "hightnorm")
ci.ppt_normmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normmid[2,], ci.low = ci.ppt_normmid[1,], ci.high = ci.ppt_normmid[3,], ci.group = "midtnorm")
ci.ppt_normlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normlow[2,], ci.low = ci.ppt_normlow[1,], ci.high = ci.ppt_normlow[3,], ci.group = "lowtnorm")
ppt_norm_pnormint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
ggplot(data = ppt_norm_pnormint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#ppt_norm and Precip_DecJanFeb interaction
Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
growthpredictionpptnorm_highPrecipNovDecJanFebMar <- growthpredictionpptnorm_lowPrecipNovDecJanFebMar <- growthpredictionpptnorm_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
  growthpredictionpptnorm_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[2]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[2]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpptnorm_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpptnorm_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[1]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[1]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Precip_NovDecJanFebMar_prediction_trlow <- exp(growthpredictionpptnorm_lowPrecipNovDecJanFebMar)
Precip_NovDecJanFebMar_prediction_trmid <- exp(growthpredictionpptnorm_midPrecipNovDecJanFebMar)
Precip_NovDecJanFebMar_prediction_trhigh <- exp(growthpredictionpptnorm_highPrecipNovDecJanFebMar)
ci.Precip_NovDecJanFebMarhigh <- apply(Precip_NovDecJanFebMar_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_NovDecJanFebMarmid <- apply(Precip_NovDecJanFebMar_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_NovDecJanFebMarlow <- apply(Precip_NovDecJanFebMar_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_NovDecJanFebMarhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "highppt_norm")
ci.Precip_NovDecJanFebMarmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midppt_norm")
ci.Precip_NovDecJanFebMarlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowppt_norm")
Precip_NovDecJanFebMar_ppt_normint <- rbind(ci.Precip_NovDecJanFebMarhigh.df, ci.Precip_NovDecJanFebMarmid.df, ci.Precip_NovDecJanFebMarlow.df)
ggplot(data = Precip_NovDecJanFebMar_ppt_normint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = Precip_NovDecJanFebMar_ppt_normint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#ppt_norm and Tmean_AprMayJun interaction
Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
ppt_norm <- mean(grow_train$ppt_norm)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
growthpredictionpnorm_highTmeanAprMayJun <- growthpredictionpnorm_lowTmeanAprMayJun <- growthpredictionpnorm_midTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
  growthpredictionpnorm_highTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[2]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[2]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpnorm_midTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpnorm_lowTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm_range[1]*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm_range[1]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowTmeanAprMayJun)
ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midTmeanAprMayJun)
ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highTmeanAprMayJun)
ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.ppt_normhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highppt_norm")
ci.ppt_normmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midppt_norm")
ci.ppt_normlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowppt_norm")
ppt_norm_Tmean_AprMayJunint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
ggplot(data = ppt_norm_Tmean_AprMayJunint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = ppt_norm_Tmean_AprMayJunint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#ppt_norm and Tmean_SepOct interaction
ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Tmean_SepOct_range <- quantile(grow_train$Tmean_SepOct, c(0.2, 0.8))
growthpredictionpnorm_highTmeanSepOct <- growthpredictionpnorm_lowTmeanSepOct <- growthpredictionpnorm_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 

for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
  growthpredictionpnorm_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[2]
  
  growthpredictionpnorm_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpnorm_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[1]
}
ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowTmeanSepOct)
ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midTmeanSepOct)
ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highTmeanSepOct)
ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.ppt_normhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normhigh[2,], ci.low = ci.ppt_normhigh[1,], ci.high = ci.ppt_normhigh[3,], ci.group = "highTmeanSepOct")
ci.ppt_normmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normmid[2,], ci.low = ci.ppt_normmid[1,], ci.high = ci.ppt_normmid[3,], ci.group = "midTmeanSepOct")
ci.ppt_normlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normlow[2,], ci.low = ci.ppt_normlow[1,], ci.high = ci.ppt_normlow[3,], ci.group = "lowTmeanSepOct")
ppt_norm_Tmean_SepOctint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
ggplot(data = ppt_norm_Tmean_SepOctint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#ppt_norm and size
ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) #setting range for tmp_normrng
ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictionpnorm_highsize <- growthpredictionpnorm_lowsize <- growthpredictionpnorm_midsize <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 

for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
  growthpredictionpnorm_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpnorm_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionpnorm_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowsize)
ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midsize)
ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highsize)
ci.ppt_normhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.ppt_normmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.ppt_normlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.ppt_normhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normhigh[2,], ci.low = ci.ppt_normhigh[1,], ci.high = ci.ppt_normhigh[3,], ci.group = "highsize")
ci.ppt_normmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normmid[2,], ci.low = ci.ppt_normmid[1,], ci.high = ci.ppt_normmid[3,], ci.group = "midTmeansize")
ci.ppt_normlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_normlow[2,], ci.low = ci.ppt_normlow[1,], ci.high = ci.ppt_normlow[3,], ci.group = "lowTmeansize")
ppt_norm_sizeint <- rbind(ci.ppt_normhigh.df, ci.ppt_normmid.df, ci.ppt_normlow.df)
ggplot(data = ppt_norm_sizeint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_norm and Precip_NovDecJanFebMar interaction
Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictiontnorm_highPrecipNovDecJanFebMar <- growthpredictiontnorm_lowPrecipNovDecJanFebMar <- growthpredictiontnorm_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
  growthpredictiontnorm_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictiontnorm_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictiontnorm_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowPrecipNovDecJanFebMar)
tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midPrecipNovDecJanFebMar)
tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highPrecipNovDecJanFebMar)
ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "hightmp_norm")
ci.tmp_normmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midtmp_norm")
ci.tmp_normlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowtmp_norm")
tmp_norm_Precip_NovDecJanFebMarint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
ggplot(data = tmp_norm_Precip_NovDecJanFebMarint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = tmp_norm_Precip_NovDecJanFebMarint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_norm and Precip_JulAug interaction
tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Precip_JulAug_range <- quantile(grow_train$Precip_JulAug, c(0.2, 0.8))
growthpredictiontnorm_highPrecipJulAug <- growthpredictiontnorm_lowPrecipJulAug <- growthpredictiontnorm_midPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
  growthpredictiontnorm_highPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug_range[2] +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug_range[2]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug_range[2]*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictiontnorm_midPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictiontnorm_lowPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug_range[1] +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug_range[1]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug_range[1]*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowPrecipJulAug)
tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midPrecipJulAug)
tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highPrecipJulAug)
ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highPrecipJulAug")
ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midPrecipJulAug")
ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowPrecipJulAug")
tmp_norm_Precip_JulAugint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
ggplot(data = tmp_norm_Precip_JulAugint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_norm and Tmean_AprMayJun interaction
tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Tmean_AprMayJun_range <- quantile(grow_train$Tmean_AprMayJun, c(0.2, 0.8))
growthpredictiontnorm_highTmeanAprMayJun <- growthpredictiontnorm_lowTmeanAprMayJun <- growthpredictiontnorm_midTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
  growthpredictiontnorm_highTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun_range[2] + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun_range[2] +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun_range[2] + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun_range[2] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun_range[2] +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun_range[2]*Tmean_SepOct
  
  growthpredictiontnorm_midTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictiontnorm_lowTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun_range[1] + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun_range[1] +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun_range[1] + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun_range[1] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun_range[1] +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun_range[1]*Tmean_SepOct
}
tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowTmeanAprMayJun)
tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midTmeanAprMayJun)
tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highTmeanAprMayJun)
ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highTmeanAprMayJun")
ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midTmeanAprMayJun")
ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowTmeanAprMayJun")
tmp_norm_Tmean_AprMayJunint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
ggplot(data = tmp_norm_Tmean_AprMayJunint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_norm and Tmean_SepOct interaction
tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Tmean_SepOct_range <- quantile(grow_train$Tmean_SepOct, c(0.2, 0.8))
growthpredictiontnorm_highTmeanSepOct <- growthpredictiontnorm_lowTmeanSepOct <- growthpredictiontnorm_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
  growthpredictiontnorm_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[2]
  
  growthpredictiontnorm_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictiontnorm_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[1]
}
tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowTmeanSepOct)
tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midTmeanSepOct)
tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highTmeanSepOct)
ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highTmeanSepOct")
ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midTmeanSepOct")
ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowTmeanSepOct")
tmp_norm_Tmean_SepOctint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
ggplot(data = tmp_norm_Tmean_SepOctint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_norm and size interaction
tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictiontnorm_highsize <- growthpredictiontnorm_lowsize <- growthpredictiontnorm_midsize <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
  growthpredictiontnorm_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictiontnorm_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictiontnorm_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowsize)
tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midsize)
tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highsize)
ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "highsize")
ci.tmp_normmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midsize")
ci.tmp_normlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowsize")
tmp_norm_sizeint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
ggplot(data = tmp_norm_sizeint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Precip_NovDecJanFebMar and Tmean_DecJanFeb interaction
Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Tmean_SepOct_range <- quantile(grow_train$Tmean_SepOct, c(0.2, 0.8))
growthpredictionPrecipNovDecJanFebMar_highTmeanSepOct <- growthpredictionPrecipNovDecJanFebMar_lowTmeanSepOct <- growthpredictionPrecipNovDecJanFebMar_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
  growthpredictionPrecipNovDecJanFebMar_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[2]
  
  growthpredictionPrecipNovDecJanFebMar_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipNovDecJanFebMar_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[1]
}
Precip_NovDecJanFebMar_prediction_trlow <- exp(growthpredictionPrecipNovDecJanFebMar_lowTmeanSepOct)
Precip_NovDecJanFebMar_prediction_trmid <- exp(growthpredictionPrecipNovDecJanFebMar_midTmeanSepOct)
Precip_NovDecJanFebMar_prediction_trhigh <- exp(growthpredictionPrecipNovDecJanFebMar_highTmeanSepOct)
ci.Precip_NovDecJanFebMarhigh <- apply(Precip_NovDecJanFebMar_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_NovDecJanFebMarmid <- apply(Precip_NovDecJanFebMar_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_NovDecJanFebMarlow <- apply(Precip_NovDecJanFebMar_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_NovDecJanFebMarhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "highTmeanSepOct")
ci.Precip_NovDecJanFebMarmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midTmeanSepOct")
ci.Precip_NovDecJanFebMarlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowTmeanSepOct")
Precip_NovDecJanFebMar_Tmean_SepOctint <- rbind(ci.Precip_NovDecJanFebMarhigh.df, ci.Precip_NovDecJanFebMarmid.df, ci.Precip_NovDecJanFebMarlow.df)
ggplot(data = Precip_NovDecJanFebMar_Tmean_SepOctint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_NovDecJanFebMar and Tmean_AprMayJun interaction
Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Precip_NovDecJanFebMar_range <- quantile(grow_train$Precip_NovDecJanFebMar, c(0.2, 0.8))
growthpredictionTmeanAprMayJun_highPrecipNovDecJanFebMar <- growthpredictionTmeanAprMayJun_lowPrecipNovDecJanFebMar <- growthpredictionTmeanAprMayJun_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
  growthpredictionTmeanAprMayJun_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[2] +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar_range[2] + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar_range[2] + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar_range[2] + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar_range[2]*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionTmeanAprMayJun_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionTmeanAprMayJun_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[1] +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar_range[1] + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar_range[1] + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar_range[1] + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar_range[1]*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Tmean_AprMayJun_prediction_trlow <- exp(growthpredictionTmeanAprMayJun_lowPrecipNovDecJanFebMar)
Tmean_AprMayJun_prediction_trmid <- exp(growthpredictionTmeanAprMayJun_midPrecipNovDecJanFebMar)
Tmean_AprMayJun_prediction_trhigh <- exp(growthpredictionTmeanAprMayJun_highPrecipNovDecJanFebMar)
ci.Tmean_AprMayJunhigh <- apply(Tmean_AprMayJun_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Tmean_AprMayJunmid <- apply(Tmean_AprMayJun_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunlow <- apply(Tmean_AprMayJun_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highPrecip_NovDecJanFebMar")
ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midPrecip_NovDecJanFebMar")
ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowPrecip_NovDecJanFebMar")
Tmean_AprMayJun_Precip_NovDecJanFebMarint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
ggplot(data = Tmean_AprMayJun_Precip_NovDecJanFebMarint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = Tmean_AprMayJun_Precip_NovDecJanFebMarint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Precip_NovDecJanFebMar and size interaction
Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictionPrecipNovDecJanFebMar_highsize <- growthpredictionPrecipNovDecJanFebMar_lowsize <- growthpredictionPrecipNovDecJanFebMar_midsize <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
  growthpredictionPrecipNovDecJanFebMar_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipNovDecJanFebMar_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionPrecipNovDecJanFebMar_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Precip_NovDecJanFebMar_prediction_trlow <- exp(growthpredictionPrecipNovDecJanFebMar_lowsize)
Precip_NovDecJanFebMar_prediction_trmid <- exp(growthpredictionPrecipNovDecJanFebMar_midsize)
Precip_NovDecJanFebMar_prediction_trhigh <- exp(growthpredictionPrecipNovDecJanFebMar_highsize)
ci.Precip_NovDecJanFebMarhigh <- apply(Precip_NovDecJanFebMar_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_NovDecJanFebMarmid <- apply(Precip_NovDecJanFebMar_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_NovDecJanFebMarlow <- apply(Precip_NovDecJanFebMar_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_NovDecJanFebMarhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "highsize")
ci.Precip_NovDecJanFebMarmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midsize")
ci.Precip_NovDecJanFebMarlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowsize")
Precip_NovDecJanFebMar_sizeint <- rbind(ci.Precip_NovDecJanFebMarhigh.df, ci.Precip_NovDecJanFebMarmid.df, ci.Precip_NovDecJanFebMarlow.df)
ggplot(data = Precip_NovDecJanFebMar_sizeint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Tmean_AprMayJun and Tmean_SepOct interaction
Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
Tmean_SepOct_range <- quantile(grow_train$Tmean_SepOct, c(0.2, 0.8))
growthpredictionTmeanAprMayJun_highTmeanSepOct <- growthpredictionTmeanAprMayJun_lowTmeanSepOct <- growthpredictionTmeanAprMayJun_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
  growthpredictionTmeanAprMayJun_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[2] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[2]
  
  growthpredictionTmeanAprMayJun_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionTmeanAprMayJun_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct_range[1] +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct_range[1]
}
Tmean_AprMayJun_prediction_trlow <- exp(growthpredictionTmeanAprMayJun_lowTmeanSepOct)
Tmean_AprMayJun_prediction_trmid <- exp(growthpredictionTmeanAprMayJun_midTmeanSepOct)
Tmean_AprMayJun_prediction_trhigh <- exp(growthpredictionTmeanAprMayJun_highTmeanSepOct)
ci.Tmean_AprMayJunhigh <- apply(Tmean_AprMayJun_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Tmean_AprMayJunmid <- apply(Tmean_AprMayJun_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunlow <- apply(Tmean_AprMayJun_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highTmeanSepOct")
ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midTmeanSepOct")
ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowTmeanSepOct")
Tmean_AprMayJun_Tmean_SepOctint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
ggplot(data = Tmean_AprMayJun_Tmean_SepOctint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Tmean_AprMayJun and size interaction
Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
x_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictionTmeanAprMayJun_highsize <- growthpredictionTmeanAprMayJun_lowsize <- growthpredictionTmeanAprMayJun_midsize <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
  growthpredictionTmeanAprMayJun_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionTmeanAprMayJun_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionTmeanAprMayJun_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
Tmean_AprMayJun_prediction_trlow <- exp(growthpredictionTmeanAprMayJun_lowsize)
Tmean_AprMayJun_prediction_trmid <- exp(growthpredictionTmeanAprMayJun_midsize)
Tmean_AprMayJun_prediction_trhigh <- exp(growthpredictionTmeanAprMayJun_highsize)
ci.Tmean_AprMayJunhigh <- apply(Tmean_AprMayJun_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Tmean_AprMayJunmid <- apply(Tmean_AprMayJun_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunlow <- apply(Tmean_AprMayJun_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highsize")
ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midsize")
ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowsize")
Tmean_AprMayJun_sizeint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
ggplot(data = Tmean_AprMayJun_sizeint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Tmean_DecJanFeb and size
sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
size <- seq(sizerng[1], sizerng[2], by = 7.5)
Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
Tmean_DecJanFeb_range <- quantile(grow_train$Tmean_DecJanFeb, c(0.2, 0.8))
growthpredictionsize_highTmeanDecJanFeb <- growthpredictionsize_lowTmeanDecJanFeb <- growthpredictionsize_midTmeanDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size)) 

for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
  growthpredictionsize_highTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[2] +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[2] + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[2] +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[2]*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[2]*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[2]*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[2]*size
  
  growthpredictionsize_midTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_lowTmeanDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb_range[1] +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb_range[1] + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb_range[1] +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb_range[1]*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb_range[1]*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb_range[1]*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb_range[1]*size
}
size_prediction_trlow <- exp(growthpredictionsize_lowTmeanDecJanFeb)
size_prediction_trmid <- exp(growthpredictionsize_midTmeanDecJanFeb)
size_prediction_trhigh <- exp(growthpredictionsize_highTmeanDecJanFeb)
ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highTmeanDecJanFeb")
ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midTmeanDecJanFeb")
ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowTmeanDecJanFeb")
size_Tmean_DecJanFebint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
ggplot(data = size_Tmean_DecJanFebint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#size and Tmean_JulAug
sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
size <- seq(sizerng[1], sizerng[2], by = 7.5)
Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
Tmean_JulAug_range <- quantile(grow_train$Tmean_JulAug, c(0.2, 0.8))
growthpredictionsize_highTmeanJulAug <- growthpredictionsize_lowTmeanJulAug <- growthpredictionsize_midTmeanJulAug <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size)) 

for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
  growthpredictionsize_highTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug_range[2] + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug_range[2] + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug_range[2]*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug_range[2]*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug_range[2]*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug_range[2]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug_range[2]*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_midTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_lowTmeanJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug_range[1] + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug_range[1] + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug_range[1]*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug_range[1]*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug_range[1]*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug_range[1]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug_range[1]*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
}
size_prediction_trlow <- exp(growthpredictionsize_lowTmeanJulAug)
size_prediction_trmid <- exp(growthpredictionsize_midTmeanJulAug)
size_prediction_trhigh <- exp(growthpredictionsize_highTmeanJulAug)
ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highTmeanJulAug")
ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midTmeanJulAug")
ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowTmeanJulAug")
size_Tmean_JulAugint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
ggplot(data = size_Tmean_JulAugint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#size and Precip_JulAug
sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
size <- seq(sizerng[1], sizerng[2], by = 7.5)
Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
Precip_JulAug_range <- quantile(grow_train$Precip_JulAug, c(0.2, 0.8))
growthpredictionsize_highPrecipJulAug <- growthpredictionsize_lowPrecipJulAug <- growthpredictionsize_midPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size)) 

for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
  growthpredictionsize_highPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[2] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug_range[2]*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug_range[2]*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug_range[2] +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_midPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_lowPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug_range[1] + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug_range[1]*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug_range[1]*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug_range[1] +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
}
size_prediction_trlow <- exp(growthpredictionsize_lowPrecipJulAug)
size_prediction_trmid <- exp(growthpredictionsize_midPrecipJulAug)
size_prediction_trhigh <- exp(growthpredictionsize_highPrecipJulAug)
ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highPrecipJulAug")
ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midPrecipJulAug")
ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowPrecipJulAug")
size_Precip_JulAugint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
ggplot(data = size_Precip_JulAugint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#size and ppt_norm
sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
size <- seq(sizerng[1], sizerng[2], by = 7.5)
Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
growthpredictionsize_highpnorm <- growthpredictionsize_lowpnorm <- growthpredictionsize_midpnorm <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size)) 

for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
  growthpredictionsize_highpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm_range[2] +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_midpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_lowpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm_range[1] +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
}
size_prediction_trlow <- exp(growthpredictionsize_lowpnorm)
size_prediction_trmid <- exp(growthpredictionsize_midpnorm)
size_prediction_trhigh <- exp(growthpredictionsize_highpnorm)
ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highpnorm")
ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midpnorm")
ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowpnorm")
size_pnormint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
ggplot(data = size_pnormint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#size and tmp_norm
sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
size <- seq(sizerng[1], sizerng[2], by = 7.5)
Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictionsize_hightnorm <- growthpredictionsize_lowtnorm <- growthpredictionsize_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size)) 

for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
  growthpredictionsize_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[2] + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm_range[2] + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm_range[2] + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm_range[2] + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[1] + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm_range[1] + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm_range[1] + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm_range[1] + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
}
size_prediction_trlow <- exp(growthpredictionsize_lowtnorm)
size_prediction_trmid <- exp(growthpredictionsize_midtnorm)
size_prediction_trhigh <- exp(growthpredictionsize_hightnorm)
ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "hightnorm")
ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midtnorm")
ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowtnorm")
size_tnormint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
ggplot(data = size_tnormint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#size and Precip_DecJanFeb
sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for tmp_normrng
size <- seq(sizerng[1], sizerng[2], by = 7.5)
Tmean_JulAug <- mean(grow_train$Tmean_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_DecJanFeb <- mean(grow_train$Precip_DecJanFeb)
Tmean_DecJanFeb <- mean(grow_train$Tmean_DecJanFeb)
Precip_DecJanFeb_range <- quantile(grow_train$Precip_DecJanFeb, c(0.2, 0.8))
growthpredictionsize_highPrecipDecJanFeb <- growthpredictionsize_lowPrecipDecJanFeb <- growthpredictionsize_midPrecipDecJanFeb <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size)) 

for(i in 1:length(plotdatainterval$u_beta_DIA_prev)){
  growthpredictionsize_highPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[2] +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[2]*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[2]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[2]*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[2]*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[2]*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_midPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
  
  growthpredictionsize_lowPrecipDecJanFeb[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_DecJanFeb"]*Precip_DecJanFeb_range[1] +
    plotdatainterval[i,"u_beta_Tmean_JulAug"]*Tmean_JulAug + plotdatainterval[i,"u_beta_Tmean_DecJanFeb"]*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size +
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size+
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_Precip_DecJanFeb_ppt_norm"]*Precip_DecJanFeb_range[1]*ppt_norm +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_tmp_norm"]*Precip_DecJanFeb_range[1]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Precip_JulAug"]*Precip_DecJanFeb_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_JulAug"]*Precip_DecJanFeb_range[1]*Tmean_JulAug + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_Tmean_DecJanFeb"]*Precip_DecJanFeb_range[1]*Tmean_DecJanFeb + 
    plotdatainterval[i,"u_beta_Precip_DecJanFeb_DIA_prev"]*Precip_DecJanFeb_range[1]*size + plotdatainterval[i,"u_beta_Tmean_JulAug_ppt_norm"]*Tmean_JulAug*ppt_norm +
    plotdatainterval[i,"u_beta_Tmean_JulAug_tmp_norm"]*Tmean_JulAug*tmp_norm + plotdatainterval[i,"u_beta_Tmean_JulAug_Precip_JulAug"]*Tmean_JulAug*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_JulAug_Tmean_DecJanFeb"]*Tmean_JulAug*Tmean_DecJanFeb +
    plotdatainterval[i,"u_beta_Tmean_JulAug_DIA_prev"]*Tmean_JulAug*size + plotdatainterval[i,"u_beta_Tmean_DecJanFeb_ppt_norm"]*Tmean_DecJanFeb*ppt_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_tmp_norm"]*Tmean_DecJanFeb*tmp_norm + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_Precip_JulAug"]*Tmean_DecJanFeb*Precip_JulAug + 
    plotdatainterval[i,"u_beta_Tmean_DecJanFeb_DIA_prev"]*Tmean_DecJanFeb*size
}
size_prediction_trlow <- exp(growthpredictionsize_lowPrecipDecJanFeb)
size_prediction_trmid <- exp(growthpredictionsize_midPrecipDecJanFeb)
size_prediction_trhigh <- exp(growthpredictionsize_highPrecipDecJanFeb)
ci.sizehigh <- apply(size_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.sizemid <- apply(size_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizelow <- apply(size_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.sizehigh.df <- data.frame(size = size, median = ci.sizehigh[2,], ci.low = ci.sizehigh[1,], ci.high = ci.sizehigh[3,], ci.group = "highPrecipDecJanFeb")
ci.sizemid.df <- data.frame(size = size, median = ci.sizemid[2,], ci.low = ci.sizemid[1,], ci.high = ci.sizemid[3,], ci.group = "midPrecipDecJanFeb")
ci.sizelow.df <- data.frame(size = size, median = ci.sizelow[2,], ci.low = ci.sizelow[1,], ci.high = ci.sizelow[3,], ci.group = "lowPrecipDecJanFeb")
size_PrecipDecJanFebint <- rbind(ci.sizehigh.df, ci.sizemid.df, ci.sizelow.df)
ggplot(data = size_PrecipDecJanFebint, aes(x = size, y = median, color = ci.group)) + geom_ribbon(aes(x = size, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)




#New interaction plots

biggrid <- expand_grid(ppt_norm = as.vector(quantile(grow_train$ppt_norm, c(0, .2, .6, .8, 1.0))), tmp_norm = as.vector(quantile(grow_train$tmp_norm, c(0, .2, .6, .8, 1.0))),
                       Precip_NovDecJanFebMar = as.vector(quantile(grow_train$Precip_NovDecJanFebMar, c(0, .2, .6, .8, 1.0))), Precip_JulAug = as.vector(quantile(grow_train$Precip_JulAug, c(0, .2, .6, .8, 1.0))),
                       Tmean_AprMayJun = as.vector(quantile(grow_train$Tmean_AprMayJun, c(0, .2, .6, .8, 1.0))), Tmean_SepOct = as.vector(quantile(grow_train$Tmean_SepOct, c(0, .2, .6, .8, 1.0))),
                       DIA_prev = as.vector(quantile(grow_train$DIA_prev, c(0, .2, .6, .8, 1.0))))

ppt_norm <- biggrid$ppt_norm
tmp_norm <- biggrid$tmp_norm
Precip_NovDecJanFebMar <- biggrid$Precip_NovDecJanFebMar
Precip_JulAug <- biggrid$Precip_JulAug
Tmean_AprMayJun <- biggrid$Tmean_AprMayJun
Tmean_SepOct <- biggrid$Tmean_SepOct
x <- biggrid$DIA_prev

growthprediction <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), nrow(biggrid)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
  growthprediction[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
}
growth_prediction <- exp(growthprediction)
ci.growthprediction <- apply(growth_prediction, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.growthprediction <- t(ci.growthprediction)
growthpredictioncbind <- cbind(ci.growthprediction, biggrid)
growthpredictionfilter2 <- growthpredictioncbind
growthpredictionfilter2$ppt_norm_character <- ifelse(growthpredictionfilter2$ppt_norm == quantile (grow_train$ppt_norm, c(0)), "lowest ppt_norm", 
                                                     ifelse(growthpredictionfilter2$ppt_norm == quantile(grow_train$ppt_norm, c(0.2)), "low ppt_norm", 
                                                            ifelse(growthpredictionfilter2$ppt_norm == quantile(grow_train$ppt_norm, c(0.6)), "mid ppt_norm", 
                                                                   ifelse(growthpredictionfilter2$ppt_norm == quantile(grow_train$ppt_norm, c(0.8)), "high ppt_norm", "highest ppt_norm"))))
growthpredictionfilter2$ppt_norm_character <- factor(growthpredictionfilter2$ppt_norm_character, levels = c("lowest ppt_norm", "low ppt_norm", "mid ppt_norm", "high ppt_norm", "highest ppt_norm"))
growthpredictionfilter2$tmp_norm_character <- ifelse(growthpredictionfilter2$tmp_norm == quantile (grow_train$tmp_norm, c(0)), "lowest tmp_norm", 
                                                     ifelse(growthpredictionfilter2$tmp_norm == quantile(grow_train$tmp_norm, c(0.2)), "low tmp_norm", 
                                                            ifelse(growthpredictionfilter2$tmp_norm == quantile(grow_train$tmp_norm, c(0.6)), "mid tmp_norm", 
                                                                   ifelse(growthpredictionfilter2$tmp_norm == quantile(grow_train$tmp_norm, c(0.8)), "high tmp_norm", "highest tmp_norm"))))
growthpredictionfilter2$tmp_norm_character <- factor(growthpredictionfilter2$tmp_norm_character, levels = c("lowest tmp_norm", "low tmp_norm", "mid tmp_norm", "high tmp_norm", "highest tmp_norm"))
growthpredictionfilter2$Precip_NovDecJanFebMar_character <- ifelse(growthpredictionfilter2$Precip_NovDecJanFebMar == quantile (grow_train$Precip_NovDecJanFebMar, c(0)), "lowest Precip_NovDecJanFebMar", 
                                                                   ifelse(growthpredictionfilter2$Precip_NovDecJanFebMar == quantile(grow_train$Precip_NovDecJanFebMar, c(0.2)), "low Precip_NovDecJanFebMar", 
                                                                          ifelse(growthpredictionfilter2$Precip_NovDecJanFebMar == quantile(grow_train$Precip_NovDecJanFebMar, c(0.6)), "mid Precip_NovDecJanFebMar", 
                                                                                 ifelse(growthpredictionfilter2$Precip_NovDecJanFebMar == quantile(grow_train$Precip_NovDecJanFebMar, c(0.8)), "high Precip_NovDecJanFebMar", "highest Precip_NovDecJanFebMar"))))
growthpredictionfilter2$Precip_NovDecJanFebMar_character <- factor(growthpredictionfilter2$Precip_NovDecJanFebMar_character, levels = c("lowest Precip_NovDecJanFebMar", "low Precip_NovDecJanFebMar", "mid Precip_NovDecJanFebMar", "high Precip_NovDecJanFebMar", "highest Precip_NovDecJanFebMar"))
growthpredictionfilter2$Precip_JulAug_character <- ifelse(growthpredictionfilter2$Precip_JulAug == quantile (grow_train$Precip_JulAug, c(0)), "lowest Precip_JulAug", 
                                                          ifelse(growthpredictionfilter2$Precip_JulAug == quantile(grow_train$Precip_JulAug, c(0.2)), "low Precip_JulAug", 
                                                                 ifelse(growthpredictionfilter2$Precip_JulAug == quantile(grow_train$Precip_JulAug, c(0.6)), "mid Precip_JulAug", 
                                                                        ifelse(growthpredictionfilter2$Precip_JulAug == quantile(grow_train$Precip_JulAug, c(0.8)), "high Precip_JulAug", "highest Precip_JulAug"))))
growthpredictionfilter2$Precip_JulAug_character <- factor(growthpredictionfilter2$Precip_JulAug_character, levels = c("lowest Precip_JulAug", "low Precip_JulAug", "mid Precip_JulAug", "high Precip_JulAug", "highest Precip_JulAug"))
growthpredictionfilter2$Tmean_AprMayJun_character <- ifelse(growthpredictionfilter2$Tmean_AprMayJun == quantile (grow_train$Tmean_AprMayJun, c(0)), "lowest Tmean_AprMayJun", 
                                                            ifelse(growthpredictionfilter2$Tmean_AprMayJun == quantile(grow_train$Tmean_AprMayJun, c(0.2)), "low Tmean_AprMayJun", 
                                                                   ifelse(growthpredictionfilter2$Tmean_AprMayJun == quantile(grow_train$Tmean_AprMayJun, c(0.6)), "mid Tmean_AprMayJun", 
                                                                          ifelse(growthpredictionfilter2$Tmean_AprMayJun == quantile(grow_train$Tmean_AprMayJun, c(0.8)), "high Tmean_AprMayJun", "highest Tmean_AprMayJun"))))
growthpredictionfilter2$Tmean_AprMayJun_character <- factor(growthpredictionfilter2$Tmean_AprMayJun_character, levels = c("lowest Tmean_AprMayJun", "low Tmean_AprMayJun", "mid Tmean_AprMayJun", "high Tmean_AprMayJun", "highest Tmean_AprMayJun"))
growthpredictionfilter2$Tmean_SepOct_character <- ifelse(growthpredictionfilter2$Tmean_SepOct == quantile (grow_train$Tmean_SepOct, c(0)), "lowest Tmean_SepOct", 
                                                         ifelse(growthpredictionfilter2$Tmean_SepOct == quantile(grow_train$Tmean_SepOct, c(0.2)), "low Tmean_SepOct", 
                                                                ifelse(growthpredictionfilter2$Tmean_SepOct == quantile(grow_train$Tmean_SepOct, c(0.6)), "mid Tmean_SepOct", 
                                                                       ifelse(growthpredictionfilter2$Tmean_SepOct == quantile(grow_train$Tmean_SepOct, c(0.8)), "high Tmean_SepOct", "highest Tmean_SepOct"))))
growthpredictionfilter2$Tmean_SepOct_character <- factor(growthpredictionfilter2$Tmean_SepOct_character, levels = c("lowest Tmean_SepOct", "low Tmean_SepOct", "mid Tmean_SepOct", "high Tmean_SepOct", "highest Tmean_SepOct"))
growthpredictionfilter2$DIA_prev_character <- ifelse(growthpredictionfilter2$DIA_prev == quantile (grow_train$DIA_prev, c(0)), "lowest DIA_prev", 
                                                     ifelse(growthpredictionfilter2$DIA_prev == quantile(grow_train$DIA_prev, c(0.2)), "low DIA_prev", 
                                                            ifelse(growthpredictionfilter2$DIA_prev == quantile(grow_train$DIA_prev, c(0.6)), "mid DIA_prev", 
                                                                   ifelse(growthpredictionfilter2$DIA_prev == quantile(grow_train$DIA_prev, c(0.8)), "high DIA_prev", "highest DIA_prev"))))
growthpredictionfilter2$DIA_prev_character <- factor(growthpredictionfilter2$DIA_prev_character, levels = c("lowest DIA_prev", "low DIA_prev", "mid DIA_prev", "high DIA_prev", "highest DIA_prev"))
growthpredictionfilter3 <- growthpredictionfilter2 %>% filter(Tmean_AprMayJun == quantile(grow_train$Tmean_AprMayJun, c(0.6)) & Tmean_SepOct == quantile(grow_train$Tmean_SepOct, c(0.6)),
                                                              DIA_prev == quantile(grow_train$DIA_prev, c(0.6)))
ggplot(data = growthpredictionfilter3, aes(x = Precip_NovDecJanFebMar, y = `50%`, color = ppt_norm_character)) + geom_ribbon(data = growthpredictionfilter3, aes(x = Precip_NovDecJanFebMar, ymin = `2.5%`, ymax = `97.5%`, fill = ppt_norm_character),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2) + facet_grid(tmp_norm_character ~ Precip_JulAug_character)








#plotting individual tree growth
#Tmean_AprMayJun and tmp_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_tnorm <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_tnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trtnorm <- exp(growthpredictionTmeanAprMayJun_tnorm)
  ci.Tmean_AprMayJuntnorm <- apply(Tmean_AprMayJun_prediction_trtnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJuntnorm.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, tmp_norm = tmp_norm, median = ci.Tmean_AprMayJuntnorm[2,], ci.low = ci.Tmean_AprMayJuntnorm[1,], ci.high = ci.Tmean_AprMayJuntnorm[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_tnormint <- rbind(ci.Tmean_AprMayJuntnorm.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_tnormint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = tmp_norm, group = ci.group)) + geom_line() + 
  mytheme + ylab("Predicted Growth") + ylim(0, 2) + scale_color_gradient2(low = "#4575b4", mid = "#fddbc7", high = "#b2182b")
#map of LATLONbin
all_states <- map_data("state")
states <- subset(all_states, region %in% c("arizona", "colorado", "utah"))
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata<-states
mapdata<-data.frame(mapdata)
ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
  geom_point(data = grow.monsoon, aes(x = LON, y = LAT, color = LATLONbin)) + theme_bw()


#Tmean_AprMayJun and ppt_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_pnorm <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_pnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trpnorm <- exp(growthpredictionTmeanAprMayJun_pnorm)
  ci.Tmean_AprMayJunpnorm <- apply(Tmean_AprMayJun_prediction_trpnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJunpnorm.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, ppt_norm = ppt_norm, median = ci.Tmean_AprMayJunpnorm[2,], ci.low = ci.Tmean_AprMayJunpnorm[1,], ci.high = ci.Tmean_AprMayJunpnorm[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_pnormint <- rbind(ci.Tmean_AprMayJunpnorm.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_pnormint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by ppt_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ppt_norm, group = ci.group)) + geom_line(alpha = 0.5) + 
  mytheme + ylab("Predicted Growth") + ylim(0, 2) + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")


#Tmean_AprMayJun and Tmean_SepOct
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Tmean_SepOct_range <- quantile(tree.grow$Tmean_SepOct, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_Tmean_SepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_Tmean_SepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trTmean_SepOct <- exp(growthpredictionTmeanAprMayJun_Tmean_SepOct)
  ci.Tmean_AprMayJunTmean_SepOct <- apply(Tmean_AprMayJun_prediction_trTmean_SepOct, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJunTmean_SepOct.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_AprMayJunTmean_SepOct[2,], ci.low = ci.Tmean_AprMayJunTmean_SepOct[1,], ci.high = ci.Tmean_AprMayJunTmean_SepOct[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_Tmean_SepOctint <- rbind(ci.Tmean_AprMayJunTmean_SepOct.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_Tmean_SepOctint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = Tmean_SepOct, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)





#Tmean_AprMayJun and Precip_JulAug
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Precip_JulAug_range <- quantile(tree.grow$Precip_JulAug, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_PrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_PrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trPrecipJulAug <- exp(growthpredictionTmeanAprMayJun_PrecipJulAug)
  ci.Tmean_AprMayJunPrecipJulAug <- apply(Tmean_AprMayJun_prediction_trPrecipJulAug, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJunPrecipJulAug.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, Precip_JulAug = Precip_JulAug, median = ci.Tmean_AprMayJunPrecipJulAug[2,], ci.low = ci.Tmean_AprMayJunPrecipJulAug[1,], ci.high = ci.Tmean_AprMayJunPrecipJulAug[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_PrecipJulAugint <- rbind(ci.Tmean_AprMayJunPrecipJulAug.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_PrecipJulAugint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = Precip_JulAug, group = ci.group)) + geom_line(alpha = 0.5) + 
  mytheme + ylab("Predicted Growth") + ylim(0, 2) + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")



#Tmean_AprMayJun and Precip_NovDecJanFebMar
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Precip_NovDecJanFebMar_range <- quantile(tree.grow$Precip_NovDecJanFebMar, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_PrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_PrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trPrecipNovDecJanFebMar <- exp(growthpredictionTmeanAprMayJun_PrecipNovDecJanFebMar)
  ci.Tmean_AprMayJunPrecipNovDecJanFebMar <- apply(Tmean_AprMayJun_prediction_trPrecipNovDecJanFebMar, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJunPrecipNovDecJanFebMar.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Tmean_AprMayJunPrecipNovDecJanFebMar[2,], ci.low = ci.Tmean_AprMayJunPrecipNovDecJanFebMar[1,], ci.high = ci.Tmean_AprMayJunPrecipNovDecJanFebMar[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_PrecipNovDecJanFebMarint <- rbind(ci.Tmean_AprMayJunPrecipNovDecJanFebMar.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_PrecipNovDecJanFebMarint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = Precip_NovDecJanFebMar, group = ci.group)) + geom_line(alpha = 0.5) + 
  mytheme + ylab("Predicted Growth") + ylim(0, 2) + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")



#Precip_NovDecJanFebMar and Precip_JulAug
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Precip_NovDecJanFebMar_range <- quantile(tree.grow$Precip_NovDecJanFebMar, c(0.2, 0.8))
  growthpredictionPrecipJulAug_PrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_PrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_JulAug_prediction_trPrecipNovDecJanFebMar <- exp(growthpredictionPrecipJulAug_PrecipNovDecJanFebMar)
  ci.Precip_JulAugPrecipNovDecJanFebMar <- apply(Precip_JulAug_prediction_trPrecipNovDecJanFebMar, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_JulAugPrecipNovDecJanFebMar.df <- data.frame(Precip_JulAug = Precip_JulAug, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_JulAugPrecipNovDecJanFebMar[2,], ci.low = ci.Precip_JulAugPrecipNovDecJanFebMar[1,], ci.high = ci.Precip_JulAugPrecipNovDecJanFebMar[3,], ci.group = tree.subset$treeCD)
  Precip_JulAug_PrecipNovDecJanFebMarint <- rbind(ci.Precip_JulAugPrecipNovDecJanFebMar.df)
  print(ind.samples[j,])
  Precip_JulAug_PrecipNovDecJanFebMarint  
}
#get.ind.tmp.response(i = 6)
Precip_JulAug_tree_response <- list()
Precip_JulAug_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_JulAug_tree_response.df <- do.call(rbind, Precip_JulAug_tree_response)
merged.response.samples <- merge(Precip_JulAug_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = Precip_NovDecJanFebMar, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_JulAug and ppt_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionPrecipJulAug_pptnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_pptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_JulAug_prediction_trpptnorm <- exp(growthpredictionPrecipJulAug_pptnorm)
  ci.Precip_JulAugpptnorm <- apply(Precip_JulAug_prediction_trpptnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_JulAugpptnorm.df <- data.frame(Precip_JulAug = Precip_JulAug, ppt_norm = ppt_norm, median = ci.Precip_JulAugpptnorm[2,], ci.low = ci.Precip_JulAugpptnorm[1,], ci.high = ci.Precip_JulAugpptnorm[3,], ci.group = tree.subset$treeCD)
  Precip_JulAug_pptnormint <- rbind(ci.Precip_JulAugpptnorm.df)
  print(ind.samples[j,])
  Precip_JulAug_pptnormint  
}
#get.ind.tmp.response(i = 6)
Precip_JulAug_tree_response <- list()
Precip_JulAug_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_JulAug_tree_response.df <- do.call(rbind, Precip_JulAug_tree_response)
merged.response.samples <- merge(Precip_JulAug_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ppt_norm, group = ci.group)) + geom_line() + 
  mytheme + ylab("Predicted Growth") + ylim(0, 2) + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")



#Precip_JulAug and tmp_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionPrecipJulAug_tmpnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_tmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
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
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = tmp_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_JulAug and Tmean_SepOct
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Tmean_SepOct_range <- quantile(tree.grow$Tmean_SepOct, c(0.2, 0.8))
  growthpredictionPrecipJulAug_TmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_TmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_JulAug_prediction_trTmeanSepOct <- exp(growthpredictionPrecipJulAug_TmeanSepOct)
  ci.Precip_JulAugTmeanSepOct <- apply(Precip_JulAug_prediction_trTmeanSepOct, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_JulAugTmeanSepOct.df <- data.frame(Precip_JulAug = Precip_JulAug, Tmean_SepOct = Tmean_SepOct, median = ci.Precip_JulAugTmeanSepOct[2,], ci.low = ci.Precip_JulAugTmeanSepOct[1,], ci.high = ci.Precip_JulAugTmeanSepOct[3,], ci.group = tree.subset$treeCD)
  Precip_JulAug_TmeanSepOctint <- rbind(ci.Precip_JulAugTmeanSepOct.df)
  print(ind.samples[j,])
  Precip_JulAug_TmeanSepOctint  
}
#get.ind.tmp.response(i = 6)
Precip_JulAug_tree_response <- list()
Precip_JulAug_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_JulAug_tree_response.df <- do.call(rbind, Precip_JulAug_tree_response)
merged.response.samples <- merge(Precip_JulAug_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = Tmean_SepOct, group = ci.group)) + geom_line(alpha = 0.5) + 
  mytheme + ylab("Predicted Growth") + ylim(0, 2) + scale_color_gradient2(low = "#4575b4", mid = "#fddbc7", high = "#b2182b")



#Precip_NovDecJanFebMar and tmp_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_NovDecJanFebMarrng <- range(tree.grow$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
  Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionPrecipNovDecJanFebMar_tmpnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionPrecipNovDecJanFebMar_tmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
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
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = tmp_norm, group = ci.group)) + geom_line(alpha = 0.5) + 
  mytheme + ylab("Predicted Growth") + ylim(0, 2) + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")



#Precip_NovDecJanFebMar and ppt_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)

growfilter <- filter(grow.monsoon, Precip_NovDecJanFebMar > 9)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_NovDecJanFebMarrng <- range(tree.grow$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
  Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionPrecipNovDecJanFebMar_pptnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionPrecipNovDecJanFebMar_pptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  } 
  Precip_NovDecJanFebMar_prediction_trpptnorm <- exp(growthpredictionPrecipNovDecJanFebMar_pptnorm)
  ci.Precip_NovDecJanFebMarpptnorm <- apply(Precip_NovDecJanFebMar_prediction_trpptnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_NovDecJanFebMarpptnorm.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, ppt_norm = ppt_norm, median = ci.Precip_NovDecJanFebMarpptnorm[2,], ci.low = ci.Precip_NovDecJanFebMarpptnorm[1,], ci.high = ci.Precip_NovDecJanFebMarpptnorm[3,], ci.group = tree.subset$treeCD)
  Precip_NovDecJanFebMar_pptnormint <- rbind(ci.Precip_NovDecJanFebMarpptnorm.df)
  print(ind.samples[j,])
  Precip_NovDecJanFebMar_pptnormint  
}
#get.ind.tmp.response(i = 6)
Precip_NovDecJanFebMar_tree_response <- list()
Precip_NovDecJanFebMar_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_NovDecJanFebMar_tree_response.df <- do.call(rbind, Precip_NovDecJanFebMar_tree_response)
merged.response.samples <- merge(Precip_NovDecJanFebMar_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ppt_norm, group = ci.group))  + geom_line(alpha = 0.5) + 
  mytheme + ylab("Predicted Growth") + ylim(0, 2) + scale_color_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#4575b4")



#Precip_NovDecJanFebMar and Tmean_SepOct
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_NovDecJanFebMarrng <- range(tree.grow$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
  Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Tmean_SepOct_range <- quantile(tree.grow$Tmean_SepOct, c(0.2, 0.8))
  growthpredictionPrecipNovDecJanFebMar_TmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionPrecipNovDecJanFebMar_TmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_NovDecJanFebMar_prediction_trTmeanSepOct <- exp(growthpredictionPrecipNovDecJanFebMar_TmeanSepOct)
  ci.Precip_NovDecJanFebMarTmeanSepOct <- apply(Precip_NovDecJanFebMar_prediction_trTmeanSepOct, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_NovDecJanFebMarTmeanSepOct.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_SepOct = Tmean_SepOct, median = ci.Precip_NovDecJanFebMarTmeanSepOct[2,], ci.low = ci.Precip_NovDecJanFebMarTmeanSepOct[1,], ci.high = ci.Precip_NovDecJanFebMarTmeanSepOct[3,], ci.group = tree.subset$treeCD)
  Precip_NovDecJanFebMar_TmeanSepOctint <- rbind(ci.Precip_NovDecJanFebMarTmeanSepOct.df)
  print(ind.samples[j,])
  Precip_NovDecJanFebMar_TmeanSepOctint  
}
#get.ind.tmp.response(i = 6)
Precip_NovDecJanFebMar_tree_response <- list()
Precip_NovDecJanFebMar_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_NovDecJanFebMar_tree_response.df <- do.call(rbind, Precip_NovDecJanFebMar_tree_response)
merged.response.samples <- merge(Precip_NovDecJanFebMar_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = Tmean_SepOct, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Tmean_SepOct and ppt_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Tmean_SepOctrng <- range(tree.grow$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionTmeanSepOct_pptnorm <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
    growthpredictionTmeanSepOct_pptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_SepOct_prediction_trpptnorm <- exp(growthpredictionTmeanSepOct_pptnorm)
  ci.Tmean_SepOctpptnorm <- apply(Tmean_SepOct_prediction_trpptnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_SepOctpptnorm.df <- data.frame(Tmean_SepOct = Tmean_SepOct, ppt_norm = ppt_norm, median = ci.Tmean_SepOctpptnorm[2,], ci.low = ci.Tmean_SepOctpptnorm[1,], ci.high = ci.Tmean_SepOctpptnorm[3,], ci.group = tree.subset$treeCD)
  Tmean_SepOct_pptnormint <- rbind(ci.Tmean_SepOctpptnorm.df)
  print(ind.samples[j,])
  Tmean_SepOct_pptnormint  
}
#get.ind.tmp.response(i = 6)
Tmean_SepOct_tree_response <- list()
Tmean_SepOct_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_SepOct_tree_response.df <- do.call(rbind, Tmean_SepOct_tree_response)
merged.response.samples <- merge(Tmean_SepOct_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = ppt_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Tmean_SepOct and ppt_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Tmean_SepOctrng <- range(tree.grow$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionTmeanSepOct_tmpnorm <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
    growthpredictionTmeanSepOct_tmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
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
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = tmp_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)





#TMPBIN monsoon
#Tmean_AprMayJun and tmp_norm
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_tnorm <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_tnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trtnorm <- exp(growthpredictionTmeanAprMayJun_tnorm)
  ci.Tmean_AprMayJuntnorm <- apply(Tmean_AprMayJun_prediction_trtnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJuntnorm.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, tmp_norm = tmp_norm, median = ci.Tmean_AprMayJuntnorm[2,], ci.low = ci.Tmean_AprMayJuntnorm[1,], ci.high = ci.Tmean_AprMayJuntnorm[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_tnormint <- rbind(ci.Tmean_AprMayJuntnorm.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_tnormint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = tmp_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#map of tmpbin
all_states <- map_data("state")
states <- subset(all_states, region %in% c("arizona", "utah", "colorado"))
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata<-states
mapdata<-data.frame(mapdata)
ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
  geom_point(data = grow.monsoon, aes(x = LON, y = LAT, color = tmpbin)) + theme_bw()



#Tmean_AprMayJun and Tmean_SepOct
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Tmean_SepOct_range <- quantile(tree.grow$Tmean_SepOct, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_TmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_TmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trTmeanSepOct <- exp(growthpredictionTmeanAprMayJun_TmeanSepOct)
  ci.Tmean_AprMayJunTmeanSepOct <- apply(Tmean_AprMayJun_prediction_trTmeanSepOct, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJunTmeanSepOct.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_AprMayJunTmeanSepOct[2,], ci.low = ci.Tmean_AprMayJunTmeanSepOct[1,], ci.high = ci.Tmean_AprMayJunTmeanSepOct[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_TmeanSepOctint <- rbind(ci.Tmean_AprMayJunTmeanSepOct.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_TmeanSepOctint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = Tmean_SepOct, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)





#Tmean_AprMayJun and Precip_JulAug
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Precip_JulAug_range <- quantile(tree.grow$Precip_JulAug, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_PrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_PrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trPrecipJulAug <- exp(growthpredictionTmeanAprMayJun_PrecipJulAug)
  ci.Tmean_AprMayJunPrecipJulAug <- apply(Tmean_AprMayJun_prediction_trPrecipJulAug, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJunPrecipJulAug.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, Precip_JulAug = Precip_JulAug, median = ci.Tmean_AprMayJunPrecipJulAug[2,], ci.low = ci.Tmean_AprMayJunPrecipJulAug[1,], ci.high = ci.Tmean_AprMayJunPrecipJulAug[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_PrecipJulAugint <- rbind(ci.Tmean_AprMayJunPrecipJulAug.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_PrecipJulAugint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = Precip_JulAug, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Tmean_AprMayJun and Precip_NovDecJanFebMar
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Precip_NovDecJanFebMar_range <- quantile(tree.grow$Precip_NovDecJanFebMar, c(0.2, 0.8))
  growthpredictionTmeanAprMayJun_PrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
    growthpredictionTmeanAprMayJun_PrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_AprMayJun_prediction_trPrecipNovDecJanFebMar <- exp(growthpredictionTmeanAprMayJun_PrecipNovDecJanFebMar)
  ci.Tmean_AprMayJunPrecipNovDecJanFebMar <- apply(Tmean_AprMayJun_prediction_trPrecipNovDecJanFebMar, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_AprMayJunPrecipNovDecJanFebMar.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Tmean_AprMayJunPrecipNovDecJanFebMar[2,], ci.low = ci.Tmean_AprMayJunPrecipNovDecJanFebMar[1,], ci.high = ci.Tmean_AprMayJunPrecipNovDecJanFebMar[3,], ci.group = tree.subset$treeCD)
  Tmean_AprMayJun_PrecipNovDecJanFebMarint <- rbind(ci.Tmean_AprMayJunPrecipNovDecJanFebMar.df)
  print(ind.samples[j,])
  Tmean_AprMayJun_PrecipNovDecJanFebMarint  
}
#get.ind.tmp.response(i = 6)
Tmean_AprMayJun_tree_response <- list()
Tmean_AprMayJun_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_AprMayJun_tree_response.df <- do.call(rbind, Tmean_AprMayJun_tree_response)
merged.response.samples <- merge(Tmean_AprMayJun_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = Precip_NovDecJanFebMar, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_NovDecJanFebMar and Precip_JulAug
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Precip_NovDecJanFebMarrng <- range(tree.grow$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
  Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Precip_JulAug_range <- quantile(tree.grow$Precip_JulAug, c(0.2, 0.8))
  growthpredictionPrecipNovDecJanFebMar_PrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionPrecipNovDecJanFebMar_PrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_NovDecJanFebMar_prediction_trPrecipJulAug <- exp(growthpredictionPrecipNovDecJanFebMar_PrecipJulAug)
  ci.Precip_NovDecJanFebMarPrecipJulAug <- apply(Precip_NovDecJanFebMar_prediction_trPrecipJulAug, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_NovDecJanFebMarPrecipJulAug.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Precip_JulAug = Precip_JulAug, median = ci.Precip_NovDecJanFebMarPrecipJulAug[2,], ci.low = ci.Precip_NovDecJanFebMarPrecipJulAug[1,], ci.high = ci.Precip_NovDecJanFebMarPrecipJulAug[3,], ci.group = tree.subset$treeCD)
  Precip_NovDecJanFebMar_PrecipJulAugint <- rbind(ci.Precip_NovDecJanFebMarPrecipJulAug.df)
  print(ind.samples[j,])
  Precip_NovDecJanFebMar_PrecipJulAugint  
}
#get.ind.tmp.response(i = 6)
Precip_NovDecJanFebMar_tree_response <- list()
Precip_NovDecJanFebMar_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_NovDecJanFebMar_tree_response.df <- do.call(rbind, Precip_NovDecJanFebMar_tree_response)
merged.response.samples <- merge(Precip_NovDecJanFebMar_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = Precip_JulAug, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_NovDecJanFebMar and ppt_norm
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Precip_NovDecJanFebMarrng <- range(tree.grow$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
  Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionPrecipNovDecJanFebMar_pptnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionPrecipNovDecJanFebMar_pptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_NovDecJanFebMar_prediction_trpptnorm <- exp(growthpredictionPrecipNovDecJanFebMar_pptnorm)
  ci.Precip_NovDecJanFebMarpptnorm <- apply(Precip_NovDecJanFebMar_prediction_trpptnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_NovDecJanFebMarpptnorm.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, ppt_norm = ppt_norm, median = ci.Precip_NovDecJanFebMarpptnorm[2,], ci.low = ci.Precip_NovDecJanFebMarpptnorm[1,], ci.high = ci.Precip_NovDecJanFebMarpptnorm[3,], ci.group = tree.subset$treeCD)
  Precip_NovDecJanFebMar_pptnormint <- rbind(ci.Precip_NovDecJanFebMarpptnorm.df)
  print(ind.samples[j,])
  Precip_NovDecJanFebMar_pptnormint  
}
#get.ind.tmp.response(i = 6)
Precip_NovDecJanFebMar_tree_response <- list()
Precip_NovDecJanFebMar_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_NovDecJanFebMar_tree_response.df <- do.call(rbind, Precip_NovDecJanFebMar_tree_response)
merged.response.samples <- merge(Precip_NovDecJanFebMar_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ppt_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_NovDecJanFebMar and tmp_norm
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Precip_NovDecJanFebMarrng <- range(tree.grow$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
  Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionPrecipNovDecJanFebMar_tmpnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionPrecipNovDecJanFebMar_tmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
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
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = tmp_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_NovDecJanFebMar and Tmean_SepOct
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Precip_NovDecJanFebMarrng <- range(tree.grow$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
  Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Tmean_SepOct_range <- quantile(tree.grow$Tmean_SepOct, c(0.2, 0.8))
  growthpredictionPrecipNovDecJanFebMar_TmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionPrecipNovDecJanFebMar_TmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_NovDecJanFebMar_prediction_trTmeanSepOct <- exp(growthpredictionPrecipNovDecJanFebMar_TmeanSepOct)
  ci.Precip_NovDecJanFebMarTmeanSepOct <- apply(Precip_NovDecJanFebMar_prediction_trTmeanSepOct, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_NovDecJanFebMarTmeanSepOct.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, Tmean_SepOct = Tmean_SepOct, median = ci.Precip_NovDecJanFebMarTmeanSepOct[2,], ci.low = ci.Precip_NovDecJanFebMarTmeanSepOct[1,], ci.high = ci.Precip_NovDecJanFebMarTmeanSepOct[3,], ci.group = tree.subset$treeCD)
  Precip_NovDecJanFebMar_TmeanSepOctint <- rbind(ci.Precip_NovDecJanFebMarTmeanSepOct.df)
  print(ind.samples[j,])
  Precip_NovDecJanFebMar_TmeanSepOctint  
}
#get.ind.tmp.response(i = 6)
Precip_NovDecJanFebMar_tree_response <- list()
Precip_NovDecJanFebMar_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_NovDecJanFebMar_tree_response.df <- do.call(rbind, Precip_NovDecJanFebMar_tree_response)
merged.response.samples <- merge(Precip_NovDecJanFebMar_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_NovDecJanFebMar, y = median, color = Tmean_SepOct, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_JulAug and ppt_norm
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionPrecipJulAug_pptnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_pptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_JulAug_prediction_trpptnorm <- exp(growthpredictionPrecipJulAug_pptnorm)
  ci.Precip_JulAugpptnorm <- apply(Precip_JulAug_prediction_trpptnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_JulAugpptnorm.df <- data.frame(Precip_JulAug = Precip_JulAug, ppt_norm = ppt_norm, median = ci.Precip_JulAugpptnorm[2,], ci.low = ci.Precip_JulAugpptnorm[1,], ci.high = ci.Precip_JulAugpptnorm[3,], ci.group = tree.subset$treeCD)
  Precip_JulAug_pptnormint <- rbind(ci.Precip_JulAugpptnorm.df)
  print(ind.samples[j,])
  Precip_JulAug_pptnormint  
}
#get.ind.tmp.response(i = 6)
Precip_JulAug_tree_response <- list()
Precip_JulAug_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_JulAug_tree_response.df <- do.call(rbind, Precip_JulAug_tree_response)
merged.response.samples <- merge(Precip_JulAug_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ppt_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_JulAug and tmp_norm
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionPrecipJulAug_tmpnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_tmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
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
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = tmp_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Precip_JulAug and Tmean_SepOct
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Tmean_SepOct <- mean(tree.grow$Tmean_SepOct)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  Tmean_SepOct_range <- quantile(tree.grow$Tmean_SepOct, c(0.2, 0.8))
  growthpredictionPrecipJulAug_TmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_TmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Precip_JulAug_prediction_trTmeanSepOct <- exp(growthpredictionPrecipJulAug_TmeanSepOct)
  ci.Precip_JulAugTmeanSepOct <- apply(Precip_JulAug_prediction_trTmeanSepOct, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Precip_JulAugTmeanSepOct.df <- data.frame(Precip_JulAug = Precip_JulAug, Tmean_SepOct = Tmean_SepOct, median = ci.Precip_JulAugTmeanSepOct[2,], ci.low = ci.Precip_JulAugTmeanSepOct[1,], ci.high = ci.Precip_JulAugTmeanSepOct[3,], ci.group = tree.subset$treeCD)
  Precip_JulAug_TmeanSepOctint <- rbind(ci.Precip_JulAugTmeanSepOct.df)
  print(ind.samples[j,])
  Precip_JulAug_TmeanSepOctint  
}
#get.ind.tmp.response(i = 6)
Precip_JulAug_tree_response <- list()
Precip_JulAug_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Precip_JulAug_tree_response.df <- do.call(rbind, Precip_JulAug_tree_response)
merged.response.samples <- merge(Precip_JulAug_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = tmp_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Tmean_SepOct and ppt_norm
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Tmean_SepOctrng <- range(tree.grow$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionTmeanSepOct_pptnorm <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
    growthpredictionTmeanSepOct_pptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
  Tmean_SepOct_prediction_trpptnorm <- exp(growthpredictionTmeanSepOct_pptnorm)
  ci.Tmean_SepOctpptnorm <- apply(Tmean_SepOct_prediction_trpptnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.Tmean_SepOctpptnorm.df <- data.frame(Tmean_SepOct = Tmean_SepOct, ppt_norm = ppt_norm, median = ci.Tmean_SepOctpptnorm[2,], ci.low = ci.Tmean_SepOctpptnorm[1,], ci.high = ci.Tmean_SepOctpptnorm[3,], ci.group = tree.subset$treeCD)
  Tmean_SepOct_pptnormint <- rbind(ci.Tmean_SepOctpptnorm.df)
  print(ind.samples[j,])
  Tmean_SepOct_pptnormint  
}
#get.ind.tmp.response(i = 6)
Tmean_SepOct_tree_response <- list()
Tmean_SepOct_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
Tmean_SepOct_tree_response.df <- do.call(rbind, Tmean_SepOct_tree_response)
merged.response.samples <- merge(Tmean_SepOct_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = ppt_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Tmean_SepOct and tmp_norm
grow.monsoon$tmpbin <- ifelse(grow.monsoon$tmp_yr > 10 & grow.monsoon$tmp_yr <= 14, "10 to 14", ifelse(grow.monsoon$tmp_yr > 14 &  grow.monsoon$tmp_yr <= 19, 
                                                                                                       "14 to 19", ifelse(grow.monsoon$tmp_yr > 6 &  grow.monsoon$tmp_yr <= 10, "6 to 10", "2 to 6")))
grow.monsoon$tmpbin <- as.vector(grow.monsoon$tmpbin)
ind.samples <- unique(grow.monsoon[,c("tmpbin", "treeCD")]) %>% group_by(tmpbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(tmpbin == tree.subset$tmpbin & treeCD == tree.subset$treeCD)
  
  Tmean_SepOctrng <- range(tree.grow$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  Precip_NovDecJanFebMar <- mean(tree.grow$Precip_NovDecJanFebMar)
  Tmean_AprMayJun <- mean(tree.grow$Tmean_AprMayJun)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionTmeanSepOct_tmpnorm <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
    growthpredictionTmeanSepOct_tmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x +
      plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*x*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*x*Precip_NovDecJanFebMar + 
      plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*x*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*x*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
      plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
      plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  }
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
merged.response.samples <- merge(Precip_JulAug_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = ci.group)) + geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = tmpbin, group = ci.group)) + geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = tmpbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Tmean_SepOct, y = median, color = tmp_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



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
all_states <- map_data("state")
states <- subset(all_states, region %in% c("arizona", "colorado", "utah", "new mexico"))
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata<-states
mapdata<-data.frame(mapdata)
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





line_type<-c(NA,"solid","dashed","dotted","dotdash")
quartiles<-c(NA,"1st quartile","2nd quartile","3rd quartile","4th quartile")

q=quantile(climdata.scaled$ppt_norm)
norms=c(0)
for (i in 2:5){
  
  n=mean(c(q[i-1],q[i]))
  norms[i]<-n
}

climdata.scaled$ppt_norm_cat[climdata.scaled$ppt_norm<=q[2]]<-quartiles[2]
climdata.scaled$ppt_norm_cat[climdata.scaled$ppt_norm>q[2] & climdata.scaled$ppt_norm<=q[3]]<-quartiles[3]
climdata.scaled$ppt_norm_cat[climdata.scaled$ppt_norm>q[3] & climdata.scaled$ppt_norm<=q[4]]<-quartiles[4]
climdata.scaled$ppt_norm_cat[climdata.scaled$ppt_norm>q[4] & climdata.scaled$ppt_norm<=q[5]]<-quartiles[5]

pptyr_plot <- ggplot(data = climdata.scaled, aes(x = ppt_yr, y = growth, col = ppt_norm_cat)) + 
  geom_point()+
  stat_function(fun=pfun_interall,args=list(temp=mean(climdata.scaled$tmp_yr),norm=norms[2]),
                aes(linetype=quartiles[2]),size=1,colour="black")+
  stat_function(fun=pfun_interall,args=list(temp=mean(climdata.scaled$tmp_yr),norm=norms[3]),
                aes(linetype=quartiles[3]),size=1,colour="black")+
  stat_function(fun=pfun_interall,args=list(temp=mean(climdata.scaled$tmp_yr),norm=norms[4]),
                aes(linetype=quartiles[4]),size=1,colour="black")+
  stat_function(fun=pfun_interall,args=list(temp=mean(climdata.scaled$tmp_yr),norm=norms[5]),
                aes(linetype=quartiles[5]),size=1,colour="black")+
  scale_colour_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                      values=c("1st quartile"="#41ae76","2nd quartile"="#238b45",
                               "3rd quartile"="#006d2c","4th quartile"="#00441b"),
                      labels=quartiles[2:5])+
  scale_linetype_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                        values=c("1st quartile"="solid","2nd quartile"="dashed",
                                 "3rd quartile"="dotdash","4th quartile"="dotted"),labels=quartiles[2:5])+
  mytheme+labs(x="Annual precipitation",y="Growth")
pptyr_plot  

tmpyr_plot <- ggplot(data = climdata.scaled, aes(x = tmp_yr, y = growth, col = ppt_norm_cat)) + 
  geom_point()+
  stat_function(fun=tfun_interall,args=list(ppt=mean(climdata.scaled$ppt_yr),norm=norms[2]),
                aes(linetype=quartiles[2]),size=0.75,colour="black")+
  stat_function(fun=tfun_interall,args=list(ppt=mean(climdata.scaled$ppt_yr),norm=norms[3]),
                aes(linetype=quartiles[3]),size=0.75,colour="black")+
  stat_function(fun=tfun_interall,args=list(ppt=mean(climdata.scaled$ppt_yr),norm=norms[4]),
                aes(linetype=quartiles[4]),size=0.75,colour="black")+
  stat_function(fun=tfun_interall,args=list(ppt=mean(climdata.scaled$ppt_yr),norm=norms[5]),
                aes(linetype=quartiles[5]),size=0.75,colour="black")+
  scale_colour_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                      values=c("1st quartile"="#41ae76","2nd quartile"="#238b45",
                               "3rd quartile"="#006d2c","4th quartile"="#00441b"),
                      labels=quartiles[2:5])+
  scale_linetype_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                        values=c("1st quartile"="solid","2nd quartile"="dashed",
                                 "3rd quartile"="dotdash","4th quartile"="dotted"),labels=quartiles[2:5])+
  mytheme+labs(x="Annual temperature",y="Growth")
tmpyr_plot  

ppt_violin <- ggplot(data = climdata.scaled, aes(x = ppt_norm_cat, y = ppt_yr)) + 
  geom_violin()+
  mytheme+labs(x="PPT norm",y="Annual precipitation")
ppt_violin

ppt_hist <- ggplot(data = climdata.scaled, aes(x = ppt_yr,fill = ppt_norm_cat)) + 
  geom_histogram()+
  scale_fill_manual("PPT Norm",breaks=c("1st quartile","2nd quartile","3rd quartile","4th quartile"),
                    values=c("1st quartile"="#41ae76","2nd quartile"="#238b45",
                             "3rd quartile"="#006d2c","4th quartile"="#00441b"),
                    labels=quartiles[2:5])+
  mytheme+labs(x="Annual precipitation",y="Count")
ppt_hist

ggsave("pptyr_plot.png",pptyr_plot)
ggsave("tmpyr_plot.png",tmpyr_plot)
ggsave("ppt_violin.png",ppt_violin)
ggsave("ppt_hist.png",ppt_hist)




write.csv(summary$summary,"./piedsummarylong.csv")

pied_grow_coef2<-summary$summary
pied_grow_coef2
save(pied_grow_coef2, grow, file="./pied_grow_coef2.rda")

save(pied_grow_coef, grow, file="./pied_grow_coef1.rda")

pial_post_grow<-as.data.frame(fit_grow)
write.csv(pial_post_grow,"./Ogle/pial_post_grow.csv",row.names=F)

pial_grow_coef<-read.csv("pial_grow_ogle_coef1.csv",row.names=1)

#Check residuals
pred<-grow_coef[1,2]+grow_coef[2,2]*grow$Size_t+grow_coef[28,2]*(grow$scovar%*%grow_coef[3:27,2])+grow_coef[29,2]*(grow$scovar%*%grow_coef[3:27,2])^2
resid<-log(grow$Growth+0.001)*10-pial_grow_coef[83:1533,1]   #pred
hist(resid)
qqnorm(resid)
qqline(resid)

#Plot of temperature coefficients
mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               legend.text=element_text(size=8),axis.text=element_text(size=8),
               axis.title.x=element_text(size=9),axis.title.y=element_text(size=9),
               axis.line.x = element_line(color="black", size = 0.3),
               axis.line.y = element_line(color="black", size = 0.3))

##Growth
grow$Growth_1000<-grow$Size_t1-grow$Size_t
q<-quantile(grow$Size_t,seq(0,1,length=4))
Size=(q[1:3]+q[2:4])/2
Dens=median(grow$Density)
grow1<-subset(grow,grow$Size_t<q[2])
grow2<-subset(grow,grow$Size_t>=q[2] & grow$Size_t<q[3])
grow3<-subset(grow,grow$Size_t>q[3])

spei_min<-apply(grow$scovar,2,min)
spei_max<-apply(grow$scovar,2,max)
spei<-matrix(NA,50,25)
for(i in 1:25){
  spei[,i]<-seq(spei_min[i],spei_max[i],length=50)
}

ncuts=25

Temperature1<-as.matrix(grow1$scovar*pial_grow_coef[9:33,1])
Temperature1<-as.data.frame(.rowSums(Temperature1,nrow(Temperature1),25))
Temperature1_2<-as.matrix((grow1$scovar^2)*pial_grow_coef[9:33,1])
Temperature1_2<-as.data.frame(.rowSums(Temperature1_2,nrow(Temperature1),25))

names(Temperature1)<-"Temps"
names(Temperature1_2)<-"Temps"
chopsize1<-cut(Temperature1$Temps,ncuts)
growbinned1<-as.vector(sapply(split(grow1$Growth_1000,chopsize1),mean))
tempbinned1<-as.vector(sapply(split(Temperature1$Temps,chopsize1),mean))
tempbinned1_2<-as.vector(sapply(split(Temperature1_2$Temps,chopsize1),mean))

g_binned1<-as.data.frame(cbind(tempbinned1,tempbinned1_2,growbinned1))
names(g_binned1)<-c("temperature","temperature2","grow")
g_binned1$size<-Size[1]

grow1<-as.data.frame(cbind(Temperature1,Temperature1_2,grow1$Growth_1000))
names(grow1)<-c("temperature","temperature2","grow")
grow1$size<-Size[1]

Temperature2<-as.matrix(grow2$scovar*pial_grow_coef[9:33,1])
Temperature2<-as.data.frame(.rowSums(Temperature2,nrow(Temperature2),25))
Temperature2_2<-as.matrix((grow2$scovar^2)*pial_grow_coef[9:33,1])
Temperature2_2<-as.data.frame(.rowSums(Temperature2_2,nrow(Temperature2),25))

names(Temperature2)<-"Temps"
names(Temperature2_2)<-"Temps"
chopsize2<-cut(Temperature2$Temps,ncuts)
growbinned2<-as.vector(sapply(split(grow2$Growth_1000,chopsize2),mean))
tempbinned2<-as.vector(sapply(split(Temperature2$Temps,chopsize2),mean))
tempbinned2_2<-as.vector(sapply(split(Temperature2_2$Temps,chopsize2),mean))

g_binned2<-as.data.frame(cbind(tempbinned2,tempbinned2_2,growbinned2))
names(g_binned2)<-c("temperature","temperature2","grow")
g_binned2$size<-Size[2]

grow2<-as.data.frame(cbind(Temperature2,Temperature2_2,grow2$Growth_1000))
names(grow2)<-c("temperature","temperature2","grow")
grow2$size<-Size[2]

Temperature3<-as.matrix(grow3$scovar*pial_grow_coef[9:33,1])
Temperature3<-as.data.frame(.rowSums(Temperature3,nrow(Temperature3),25))
Temperature3_2<-as.matrix((grow3$scovar^2)*pial_grow_coef[9:33,1])
Temperature3_2<-as.data.frame(.rowSums(Temperature3_2,nrow(Temperature3),25))

names(Temperature3)<-"Temps"
names(Temperature3_2)<-"Temps"
chopsize3<-cut(Temperature3$Temps,ncuts)
growbinned3<-as.vector(sapply(split(grow3$Growth_1000,chopsize3),mean))
tempbinned3<-as.vector(sapply(split(Temperature3$Temps,chopsize3),mean))
tempbinned3_2<-as.vector(sapply(split(Temperature3_2$Temps,chopsize3),mean))

g_binned3<-as.data.frame(cbind(tempbinned3,tempbinned3_2,growbinned3))
names(g_binned3)<-c("temperature","temperature2","grow")
g_binned3$size<-Size[3]

grow3<-as.data.frame(cbind(Temperature3,Temperature3_2,grow3$Growth_1000))
names(grow3)<-c("temperature","temperature2","grow")
grow3$size<-Size[3]

grow_plot_data<-rbind(grow1,grow2,grow3)
grow_binned_plot_data<-rbind(g_binned1,g_binned2,g_binned3)

g_fun<-function(spei,size){
  spei1<-(spei/100)%*%pial_grow_coef[9:33,1]
  spei2<-((spei/100)^2)%*%pial_grow_coef[9:33,1]
  size=size
  grow<-(pial_grow_coef[1,1]+pial_grow_coef[8,1]*size+
           pial_grow_coef[34,1]*spei1+pial_grow_coef[35,1]*spei2+
           pial_grow_coef[36,1]*spei1*size+pial_grow_coef[37,1]*spei2*size)-size
  return(data.frame(spei=spei1*100,grow=grow,size=size,row.names=NULL))
}

Growth1<-g_fun(spei=spei,size=Size[1])
Growth2<-g_fun(spei=spei,size=Size[2])
Growth3<-g_fun(spei=spei,size=Size[3])

Growth<-rbind(Growth1,Growth2,Growth3)

grow_plot<-ggplot(grow_binned_plot_data,aes(x=temperature,y=grow,col=size))+
  labs(x="SPEI",y="Growth (cm)")+
  geom_point()+
  geom_line(data=Growth,aes(x=spei,y=grow,group=size))+
  scale_colour_viridis("Size")+
  coord_cartesian(xlim=c(min(grow_plot_data$temperature,na.rm=T),max(grow_plot_data$temperature,na.rm=T)))+
  #ylim=c((min(grow_plot_data$size,na.rm=T))+5,(max(grow_plot_data$size,na.rm=T)+1)))+
  mytheme

png(file="./Output/grow_spei.png",400,360,type="cairo")
grow_plot
dev.off()
