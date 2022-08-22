### Stan models for PIAL growth using technique in Ogle et al., Ecology Letters to test for climate effects
## Sharmila Dey
# 22 June 2020
setwd("/home/rstudio")
load("/home/rstudio/data/pied_grow_coef2.rda")
fit_grow <- readRDS(url("https://de.cyverse.org/dl/d/888FD6F6-DAEA-46AE-BDF4-036A708990CC/log_normal_monsoonoos_pptoos.RDS"))
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


PIED.all <- read.csv("data/pied_all_growth_v3.csv")
full.ppt.tmean.norms <- read.csv("data/pied_all_tmean_ppt_v3.csv")
grow.new <- merge(PIED.all, full.ppt.tmean.norms, by.x = c("name", "year", "LON", "LAT"), by.y = c("name", "year", "lon", "lat"))

#grow.new <- read.csv("data/pied_all_tmean_ppt_v3.csv")

#grow.new <- merge(grow, newclimate, by.x = c("LON", "LAT", "name", "year"), by.y = c("lon", "lat", "name", "year"))

grow.monsoon<-na.omit(grow.new) %>% 
  mutate_at(scale, .vars = vars(Precip_JulAug, Precip_NovDecJanFebMar, tmp_norm, ppt_norm, tmp_yr)) %>%
  arrange(PLOT,SUBP,name) %>%
  mutate(PlotCD=as.numeric(factor(PLOT, levels = unique(PLOT))),treeCD=as.numeric(factor(name,levels=unique(name))),
         growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2))

split=0.20
trainIndex <- createDataPartition(grow.monsoon$name, p=split, list=FALSE)
grow_test <- grow.monsoon[trainIndex,]
grow_train <- grow.monsoon[-trainIndex,]



xG<-as.matrix(cbind(grow_train$ppt_norm, grow_train$tmp_norm, grow_train$Precip_JulAug, grow_train$Precip_NovDecJanFebMar,
                    grow_train$tmp_yr, grow_train$DIA_prev,
                    grow_train$ppt_norm*grow_train$tmp_norm, grow_train$ppt_norm*grow_train$tmp_yr, 
                    grow_train$ppt_norm*grow_train$Precip_JulAug, grow_train$ppt_norm*grow_train$DIA_prev,
                    grow_train$Precip_JulAug*grow_train$tmp_yr, grow_train$Precip_JulAug*grow_train$tmp_norm, grow_train$Precip_JulAug*grow_train$DIA_prev, 
                    grow_train$tmp_norm*grow_train$tmp_yr, grow_train$tmp_norm*grow_train$DIA_prev, 
                    grow_train$tmp_yr*grow_train$DIA_prev, grow_train$Precip_NovDecJanFebMar*grow_train$ppt_norm,
                    grow_train$Precip_NovDecJanFebMar*grow_train$tmp_norm, grow_train$Precip_NovDecJanFebMar*grow_train$Precip_JulAug,
                    grow_train$Precip_NovDecJanFebMar*grow_train$tmp_yr, grow_train$Precip_NovDecJanFebMar*grow_train$DIA_prev))
xGtest<-as.matrix(cbind(grow_test$ppt_norm, grow_test$tmp_norm, grow_test$Precip_JulAug, grow_test$Precip_NovDecJanFebMar,
                        grow_test$tmp_yr, grow_test$DIA_prev,
                        grow_test$ppt_norm*grow_test$tmp_norm, grow_test$ppt_norm*grow_test$tmp_yr, 
                        grow_test$ppt_norm*grow_test$Precip_JulAug, grow_test$ppt_norm*grow_test$DIA_prev,
                        grow_test$Precip_JulAug*grow_test$tmp_yr, grow_test$Precip_JulAug*grow_test$tmp_norm, grow_test$Precip_JulAug*grow_test$DIA_prev, 
                        grow_test$tmp_norm*grow_test$tmp_yr, grow_test$tmp_norm*grow_test$DIA_prev, 
                        grow_test$tmp_yr*grow_test$DIA_prev, grow_test$Precip_NovDecJanFebMar*grow_test$ppt_norm,
                        grow_test$Precip_NovDecJanFebMar*grow_test$tmp_norm, grow_test$Precip_NovDecJanFebMar*grow_test$Precip_JulAug,
                        grow_test$Precip_NovDecJanFebMar*grow_test$tmp_yr, grow_test$Precip_NovDecJanFebMar*grow_test$DIA_prev))
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


sink("pied_grow.stan")
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



fit_grow <- stan(file = 'pied_grow.stan', data = pied_dat, 
                 iter = 5000, warmup = 1000, chains = 3)
saveRDS(fit_grow, file = "log_normal_monsoonoos_pptoos.RDS")
summary<-summary(fit_grow)
summary

plotdata<-select(as.data.frame(fit_grow),"yrep[1]":"yrep[8741]")
plotdatainterval<-select(as.data.frame(fit_grow), "u_beta[1]":"u_beta[21]")
colnames(plotdatainterval) <- c("u_beta_ppt_norm", "u_beta_tmp_norm", "u_beta_Precip_JulAug", "u_beta_Precip_NovDecJanFebMar",
                                "u_beta_tmp_yr", "u_beta_DIA_prev",
                                "u_beta_ppt_norm_tmp_norm", "u_beta_ppt_norm_tmp_yr", 
                                "u_beta_ppt_norm_Precip_JulAug", "u_beta_ppt_norm_DIA_prev",
                                "u_beta_Precip_JulAug_tmp_yr", "u_beta_Precip_JulAug_tmp_norm", "u_beta_Precip_JulAug_DIA_prev", 
                                "u_beta_tmp_norm_tmp_yr", "u_beta_tmp_norm_DIA_prev", 
                                "u_beta_tmp_yr_DIA_prev", "u_beta_Precip_NovDecJanFebMar_ppt_norm",
                                "u_beta_Precip_NovDecJanFebMar_tmp_norm", "u_beta_Precip_NovDecJanFebMar_Precip_JulAug",
                                "u_beta_Precip_NovDecJanFebMar_tmp_yr", "u_beta_Precip_NovDecJanFebMar_DIA_prev")
ppc_dens_overlay(yGtest, as.matrix(plotdata))

ext_fit <- rstan::extract(fit_grow)
yrep <- ext_fit$yrep
#yrep <- exp(yrep)
mean.pred <- apply(ext_fit$yrep, 2, median)
p.o.df <- data.frame(predicted = exp(mean.pred), observed = exp(grow_test$loggrowth), error = (exp(mean.pred) - exp(grow_test$loggrowth)))
meansqrd <- (mean(p.o.df$error))^2
ggplot(p.o.df, aes(predicted, observed)) + geom_point(alpha = 0.1) + geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dotted") +
  ylim(0, 10) + xlim(0,10)


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

## Subset posterior predicitve plot by tmp_norm
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
tmp_yr <- mean(grow_train$tmp_yr)
growthpredictionpptnorm <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm))
pfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,tmp_yr){
  for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
    growthpredictionpptnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  }
  growthpredictionpptnorm
}
ppt_norm_prediction <- pfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, tmp_yr = tmp_yr)
ppt_norm_prediction_tr <- exp(ppt_norm_prediction)
ci.ppt_norm <- apply(ppt_norm_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.ppt_norm.df <- data.frame(ppt_norm = ppt_norm, median = ci.ppt_norm[2,], ci.low = ci.ppt_norm[1,], ci.high = ci.ppt_norm[3,])
ggplot() + 
  geom_ribbon(data = ci.ppt_norm.df, aes(x = ppt_norm, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.ppt_norm.df, aes(x = ppt_norm, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "tmp_norm")]), aes(x = ppt_norm, color = tmp_norm))


#effect of tmp_norm
tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.35)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
tmp_yr <- mean(grow_train$tmp_yr)
growthpredictiontmpnorm <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm))
tfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,tmp_yr){
  for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
    growthpredictiontmpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  }
  growthpredictiontmpnorm
}
tmp_norm_prediction <- tfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, tmp_yr = tmp_yr)
tmp_norm_prediction_tr <- exp(tmp_norm_prediction)
ci.tmp_norm <- apply(tmp_norm_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_norm.df <- data.frame(tmp_norm = tmp_norm, median = ci.tmp_norm[2,], ci.low = ci.tmp_norm[1,], ci.high = ci.tmp_norm[3,])
ggplot() + 
  geom_ribbon(data = ci.tmp_norm.df, aes(x = tmp_norm, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.tmp_norm.df, aes(x = tmp_norm, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "ppt_norm", "tmp_norm")]), aes(x = tmp_norm, color = ppt_norm))


#effect of Precip_JulAug
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
tmp_yr <- mean(grow_train$tmp_yr)
growthpredictionpptyr <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug))
pyfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,tmp_yr){
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionpptyr[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  }
  growthpredictionpptyr
}
Precip_JulAug_prediction <- pyfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, tmp_yr = tmp_yr)
Precip_JulAug_prediction_tr <- exp(Precip_JulAug_prediction)
ci.Precip_JulAug <- apply(Precip_JulAug_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAug.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAug[2,], ci.low = ci.Precip_JulAug[1,], ci.high = ci.Precip_JulAug[3,])
ggplot() + 
  geom_ribbon(data = ci.Precip_JulAug.df, aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.Precip_JulAug.df, aes(x = Precip_JulAug, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "Precip_JulAug", "ppt_norm")]), aes(x = Precip_JulAug, color = ppt_norm))



#effect of Precip_NovDecJanFebMar
Precip_NovDecJanFebMarrng <- range(grow_train$Precip_NovDecJanFebMar,na.rm = TRUE) #setting range for tmp_normrng
Precip_NovDecJanFebMar <- seq(Precip_NovDecJanFebMarrng[1], Precip_NovDecJanFebMarrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
tmp_yr <- mean(grow_train$tmp_yr)
growthpredictionpndjfm <- matrix(NA, length(plotdatainterval$u_beta_Precip_NovDecJanFebMar), length(Precip_NovDecJanFebMar))
pndjfmfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,tmp_yr){
  for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
    growthpredictionpndjfm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  }
  growthpredictionpndjfm
}
Precip_NovDecJanFebMar_prediction <- pndjfmfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, tmp_yr = tmp_yr)
Precip_NovDecJanFebMar_prediction_tr <- exp(Precip_NovDecJanFebMar_prediction)
ci.Precip_NovDecJanFebMar <- apply(Precip_NovDecJanFebMar_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_NovDecJanFebMar.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMar[2,], ci.low = ci.Precip_NovDecJanFebMar[1,], ci.high = ci.Precip_NovDecJanFebMar[3,])
ggplot() + 
  geom_ribbon(data = ci.Precip_NovDecJanFebMar.df, aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.Precip_NovDecJanFebMar.df, aes(x = Precip_NovDecJanFebMar, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "Precip_NovDecJanFebMar", "tmp_norm")]), aes(x = Precip_NovDecJanFebMar, color = tmp_norm))



#effect of tmp_yr
tmp_yrrng <- range(grow_train$tmp_yr,na.rm = TRUE) #setting range for tmp_normrng
tmp_yr <- seq(tmp_yrrng[1], tmp_yrrng[2], by = 0.3)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
tmp_norm <- mean(grow_train$tmp_norm)
growthpredictiontmpyr <- matrix(NA, length(plotdatainterval$u_beta_tmp_yr), length(tmp_yr))
tyfun_plotdatainterval<-function(x,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebMar,tmp_yr){
  for(i in 1:length(plotdatainterval$u_beta_tmp_yr)){
    growthpredictiontmpyr[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  }
  growthpredictiontmpyr
}
tmp_yr_prediction <- tyfun_plotdatainterval(x = x, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, tmp_yr = tmp_yr)
tmp_yr_prediction_tr <- exp(tmp_yr_prediction)
ci.tmp_yr <- apply(tmp_yr_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_yr.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yr[2,], ci.low = ci.tmp_yr[1,], ci.high = ci.tmp_yr[3,])
ggplot() + 
  geom_ribbon(data = ci.tmp_yr.df, aes(x = tmp_yr, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.tmp_yr.df, aes(x = tmp_yr, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "tmp_yr", "tmp_norm")]), aes(x = tmp_yr, color = tmp_norm))


#effect of DIA_prev
sizerng <- range(grow_train$DIA_prev,na.rm = TRUE) #setting range for sizerng
size <- seq(sizerng[1], sizerng[2], by = 9)
tmp_yr <- mean(grow_train$tmp_yr)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
tmp_norm <- mean(grow_train$tmp_norm)
growthpredictionsize <- matrix(NA, length(plotdatainterval$u_beta_DIA_prev), length(size))
sizefun_plotdatainterval<-function(size,tmp_norm,ppt_norm,Precip_JulAug,Precip_NovDecJanFebmar,tmp_yr){
  for(i in 1:length(plotdatainterval$u_beta_tmp_yr)){
    growthpredictionsize[i,] <-  plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  }
  growthpredictionsize
}
size_prediction <- sizefun_plotdatainterval(size = size, tmp_norm = tmp_norm, ppt_norm = ppt_norm, Precip_JulAug = Precip_JulAug, Precip_NovDecJanFebmar = Precip_NovDecJanFebmar, tmp_yr = tmp_yr)
size_prediction_tr <- exp(size_prediction)
ci.size <- apply(size_prediction_tr, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.size.df <- data.frame(size = size, median = ci.size[2,], ci.low = ci.size[1,], ci.high = ci.size[3,])
ggplot() + 
  geom_ribbon(data = ci.size.df, aes(x = size, ymin = ci.low, ymax = ci.high), alpha = 0.75, fill = "cadetblue2") + 
  geom_line(data = ci.size.df, aes(x = size, y = median), color = "indianred2") + mytheme + ylab("Predicted Growth") + ylim(0, 2) + 
  geom_rug(data = unique(grow_train[,c("LAT", "LON", "DIA_prev", "tmp_norm")]), aes(x = size, color = tmp_norm))


#Precip_JulAug and tmp_yr interaction
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
tmp_yr <- mean(grow_train$tmp_yr)
tmp_yr_range <- quantile(grow_train$tmp_yr, c(0.2, 0.8))
growthpredictionPrecipJulAug_hightmpyr <- growthpredictionPrecipJulAug_lowtmpyr <- growthpredictionPrecipJulAug_midtmpyr <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_hightmpyr[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr_range[2] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr_range[2] + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr_range[2] + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr_range[2] + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr_range[2]*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictionPrecipJulAug_midtmpyr[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictionPrecipJulAug_lowtmpyr[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr_range[1] +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr_range[1] + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr_range[1] + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr_range[1] + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr_range[1]*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
}
Precip_JulAug_prediction_trlow <- exp(growthpredictionPrecipJulAug_lowtmpyr)
Precip_JulAug_prediction_trmid <- exp(growthpredictionPrecipJulAug_midtmpyr)
Precip_JulAug_prediction_trhigh <- exp(growthpredictionPrecipJulAug_hightmpyr)
ci.Precip_JulAughigh <- apply(Precip_JulAug_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugmid <- apply(Precip_JulAug_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAuglow <- apply(Precip_JulAug_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAughigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAughigh[2,], ci.low = ci.Precip_JulAughigh[1,], ci.high = ci.Precip_JulAughigh[3,], ci.group = "hightmp_yr")
ci.Precip_JulAugmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugmid[2,], ci.low = ci.Precip_JulAugmid[1,], ci.high = ci.Precip_JulAugmid[3,], ci.group = "midtmp_yr")
ci.Precip_JulAuglow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAuglow[2,], ci.low = ci.Precip_JulAuglow[1,], ci.high = ci.Precip_JulAuglow[3,], ci.group = "lowtmp_yr")
Precip_JulAug_tmp_yrint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
ggplot(data = Precip_JulAug_tmp_yrint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_JulAug and ppt_norm interaction
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
tmp_yr <- mean(grow_train$tmp_yr)
ppt_norm_range <- quantile(grow_train$tmp_yr, c(0.2, 0.8))
growthpredictionPrecipJulAug_highpnorm <- growthpredictionPrecipJulAug_lowpnorm <- growthpredictionPrecipJulAug_midpnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_highpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm_range[2]*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictionPrecipJulAug_midpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictionPrecipJulAug_lowpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm_range[1]*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
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
Precip_JulAug_ppt_normint <- rbind(ci.Precip_JulAughigh.df, ci.Precip_JulAugmid.df, ci.Precip_JulAuglow.df)
ggplot(data = Precip_JulAug_ppt_normint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Precip_JulAug and tmp_norm interaction
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
tmp_yr <- mean(grow_train$tmp_yr)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictionPrecipJulAug_hightnorm <- growthpredictionPrecipJulAug_lowtnorm <- growthpredictionPrecipJulAug_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm_range[2]*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictionPrecipJulAug_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictionPrecipJulAug_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm_range[1]*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
}
Precip_JulAug_prediction_trtlow <- exp(growthpredictionPrecipJulAug_lowtnorm)
Precip_JulAug_prediction_trtmid <- exp(growthpredictionPrecipJulAug_midtnorm)
Precip_JulAug_prediction_trthigh <- exp(growthpredictionPrecipJulAug_hightnorm)
ci.Precip_JulAugthigh <- apply(Precip_JulAug_prediction_trthigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugtmid <- apply(Precip_JulAug_prediction_trtmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAugtlow <- apply(Precip_JulAug_prediction_trtlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAugthigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugthigh[2,], ci.low = ci.Precip_JulAugthigh[1,], ci.high = ci.Precip_JulAugthigh[3,], ci.group = "hightmp_norm")
ci.Precip_JulAugtmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugtmid[2,], ci.low = ci.Precip_JulAugtmid[1,], ci.high = ci.Precip_JulAugtmid[3,], ci.group = "midtmp_norm")
ci.Precip_JulAugtlow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugtlow[2,], ci.low = ci.Precip_JulAugtlow[1,], ci.high = ci.Precip_JulAugtlow[3,], ci.group = "lowtmp_norm")
Precip_JulAug_tmp_normint <- rbind(ci.Precip_JulAugthigh.df, ci.Precip_JulAugtmid.df, ci.Precip_JulAugtlow.df)
ggplot(data = Precip_JulAug_tmp_normint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_JulAug and Precip_NovDecJanFebMar interaction
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
tmp_yr <- mean(grow_train$tmp_yr)
Precip_NovDecJanFebMar_range <- quantile(grow_train$Precip_NovDecJanFebMar, c(0.2, 0.8))
growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar <- growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar <- growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[2] +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar_range[2]*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar_range[2]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar_range[2]*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar_range[2]*x
  
  growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[1] +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar_range[1]*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar_range[1]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar_range[1]*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar_range[1]*x
}
Precip_JulAug_prediction_trtlow <- exp(growthpredictionPrecipJulAug_lowPrecipNovDecJanFebMar)
Precip_JulAug_prediction_trtmid <- exp(growthpredictionPrecipJulAug_midPrecipNovDecJanFebMar)
Precip_JulAug_prediction_trthigh <- exp(growthpredictionPrecipJulAug_highPrecipNovDecJanFebMar)
ci.Precip_JulAugthigh <- apply(Precip_JulAug_prediction_trthigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugtmid <- apply(Precip_JulAug_prediction_trtmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAugtlow <- apply(Precip_JulAug_prediction_trtlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAugthigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugthigh[2,], ci.low = ci.Precip_JulAugthigh[1,], ci.high = ci.Precip_JulAugthigh[3,], ci.group = "highPrecip_NovDecJanFebMar")
ci.Precip_JulAugtmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugtmid[2,], ci.low = ci.Precip_JulAugtmid[1,], ci.high = ci.Precip_JulAugtmid[3,], ci.group = "midPrecip_NovDecJanFebMar")
ci.Precip_JulAugtlow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugtlow[2,], ci.low = ci.Precip_JulAugtlow[1,], ci.high = ci.Precip_JulAugtlow[3,], ci.group = "lowPrecip_NovDecJanFebMar")
Precip_JulAug_PrecipNovDecJanFebMarint <- rbind(ci.Precip_JulAugthigh.df, ci.Precip_JulAugtmid.df, ci.Precip_JulAugtlow.df)
ggplot(data = Precip_JulAug_PrecipNovDecJanFebMarint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#Precip_JulAug and size interaction
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
size <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
tmp_yr <- mean(grow_train$tmp_yr)
size_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictionPrecipJulAug_highsize <- growthpredictionPrecipJulAug_lowsize <- growthpredictionPrecipJulAug_midsize <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictionPrecipJulAug_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[2]+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size_range[2] + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[2] + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size_range[2] +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size_range[2]
  
  growthpredictionPrecipJulAug_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictionPrecipJulAug_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[1]+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size_range[1] + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[1] + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size_range[1] +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size_range[1]
}
Precip_JulAug_prediction_trsizelow <- exp(growthpredictionPrecipJulAug_lowsize)
Precip_JulAug_prediction_trsizemid <- exp(growthpredictionPrecipJulAug_midsize)
Precip_JulAug_prediction_trsizehigh <- exp(growthpredictionPrecipJulAug_highsize)
ci.Precip_JulAugsizehigh <- apply(Precip_JulAug_prediction_trsizehigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_JulAugsizemid <- apply(Precip_JulAug_prediction_trsizemid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAugsizelow <- apply(Precip_JulAug_prediction_trsizelow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_JulAugsizehigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugsizehigh[2,], ci.low = ci.Precip_JulAugsizehigh[1,], ci.high = ci.Precip_JulAugsizehigh[3,], ci.group = "highsize")
ci.Precip_JulAugsizemid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugsizemid[2,], ci.low = ci.Precip_JulAugsizemid[1,], ci.high = ci.Precip_JulAugsizemid[3,], ci.group = "midsize")
ci.Precip_JulAugsizelow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.Precip_JulAugsizelow[2,], ci.low = ci.Precip_JulAugsizelow[1,], ci.high = ci.Precip_JulAugsizelow[3,], ci.group = "lowsize")
Precip_JulAug_sizeint <- rbind(ci.Precip_JulAugsizehigh.df, ci.Precip_JulAugsizemid.df, ci.Precip_JulAugsizelow.df)
ggplot(data = Precip_JulAug_sizeint, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#tmp_yr and ppt_norm interaction
tmp_yrrng <- range(grow_train$tmp_yr,na.rm = TRUE) #setting range for tmp_normrng
tmp_yr <- seq(tmp_yrrng[1], tmp_yrrng[2], by = 0.50)
size <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
growthpredictiontmpyr_highpnorm <- growthpredictiontmpyr_lowpnorm <- growthpredictiontmpyr_midpnorm <- matrix(NA, length(plotdatainterval$u_beta_tmp_yr), length(tmp_yr)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_yr)){
  growthpredictiontmpyr_highpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[2]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm_range[2]*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[2]*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictiontmpyr_midpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictiontmpyr_lowpnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm_range[1]*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm_range[1]*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm_range[1]*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
}
tmp_yr_prediction_trpnormlow <- exp(growthpredictiontmpyr_lowpnorm)
tmp_yr_prediction_trpnormmid <- exp(growthpredictiontmpyr_midpnorm)
tmp_yr_prediction_trpnormhigh <- exp(growthpredictiontmpyr_highpnorm)
ci.tmp_yrpnormhigh <- apply(tmp_yr_prediction_trpnormhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_yrpnormmid <- apply(tmp_yr_prediction_trpnormmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_yrpnormlow <- apply(tmp_yr_prediction_trpnormlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_yrpnormhigh.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrpnormhigh[2,], ci.low = ci.tmp_yrpnormhigh[1,], ci.high = ci.tmp_yrpnormhigh[3,], ci.group = "highpnorm")
ci.tmp_yrpnormmid.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrpnormmid[2,], ci.low = ci.tmp_yrpnormmid[1,], ci.high = ci.tmp_yrpnormmid[3,], ci.group = "midpnorm")
ci.tmp_yrpnormlow.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrpnormlow[2,], ci.low = ci.tmp_yrpnormlow[1,], ci.high = ci.tmp_yrpnormlow[3,], ci.group = "lowpnorm")
tmp_yr_pnormint <- rbind(ci.tmp_yrpnormhigh.df, ci.tmp_yrpnormmid.df, ci.tmp_yrpnormlow.df)
ggplot(data = tmp_yr_pnormint, aes(x = tmp_yr, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)




#tmp_yr and tmp_norm interaction
tmp_yrrng <- range(grow_train$tmp_yr,na.rm = TRUE) #setting range for tmp_normrng
tmp_yr <- seq(tmp_yrrng[1], tmp_yrrng[2], by = 0.50)
size <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictiontmpyr_hightnorm <- growthpredictiontmpyr_lowtnorm <- growthpredictiontmpyr_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_tmp_yr), length(tmp_yr)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_yr)){
  growthpredictiontmpyr_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm_range[2]*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictiontmpyr_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictiontmpyr_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm_range[1]*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
}
tmp_yr_prediction_trtnormlow <- exp(growthpredictiontmpyr_lowtnorm)
tmp_yr_prediction_trtnormmid <- exp(growthpredictiontmpyr_midtnorm)
tmp_yr_prediction_trtnormhigh <- exp(growthpredictiontmpyr_hightnorm)
ci.tmp_yrtnormhigh <- apply(tmp_yr_prediction_trtnormhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_yrtnormmid <- apply(tmp_yr_prediction_trtnormmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_yrtnormlow <- apply(tmp_yr_prediction_trtnormlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_yrtnormhigh.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrtnormhigh[2,], ci.low = ci.tmp_yrtnormhigh[1,], ci.high = ci.tmp_yrtnormhigh[3,], ci.group = "hightnorm")
ci.tmp_yrtnormmid.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrtnormmid[2,], ci.low = ci.tmp_yrtnormmid[1,], ci.high = ci.tmp_yrtnormmid[3,], ci.group = "midtnorm")
ci.tmp_yrtnormlow.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrtnormlow[2,], ci.low = ci.tmp_yrtnormlow[1,], ci.high = ci.tmp_yrtnormlow[3,], ci.group = "lowtnorm")
tmp_yr_tnormint <- rbind(ci.tmp_yrtnormhigh.df, ci.tmp_yrtnormmid.df, ci.tmp_yrtnormlow.df)
ggplot(data = tmp_yr_tnormint, aes(x = tmp_yr, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#tmp_yr and Precip_NovDecJanFebMar interaction
tmp_yrrng <- range(grow_train$tmp_yr,na.rm = TRUE) #setting range for tmp_normrng
tmp_yr <- seq(tmp_yrrng[1], tmp_yrrng[2], by = 0.50)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
x <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar_range <- quantile(grow_train$Precip_NovDecJanFebMar, c(0.2, 0.8))
growthpredictiontmpyr_highPrecipNovDecJanFebMar <- growthpredictiontmpyr_lowPrecipNovDecJanFebMar <- growthpredictiontmpyr_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_tmp_yr), length(tmp_yr)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_yr)){
  growthpredictiontmpyr_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[2] +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar_range[2]*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar_range[2]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar_range[2]*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar_range[2]*x
  
  growthpredictiontmpyr_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*x
  
  growthpredictiontmpyr_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[1] +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*x + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*x + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*x + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*x +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar_range[1]*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar_range[1]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar_range[1]*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar_range[1]*x
}
tmp_yr_prediction_trPrecipNovDecJanFebMarlow <- exp(growthpredictiontmpyr_lowPrecipNovDecJanFebMar)
tmp_yr_prediction_trPrecipNovDecJanFebMarmid <- exp(growthpredictiontmpyr_midPrecipNovDecJanFebMar)
tmp_yr_prediction_trPrecipNovDecJanFebMarhigh <- exp(growthpredictiontmpyr_highPrecipNovDecJanFebMar)
ci.tmp_yrPrecipNovDecJanFebMarhigh <- apply(tmp_yr_prediction_trPrecipNovDecJanFebMarhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_yrPrecipNovDecJanFebMarmid <- apply(tmp_yr_prediction_trPrecipNovDecJanFebMarmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_yrPrecipNovDecJanFebMarlow <- apply(tmp_yr_prediction_trPrecipNovDecJanFebMarlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_yrPrecipNovDecJanFebMarhigh.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrtnormhigh[2,], ci.low = ci.tmp_yrtnormhigh[1,], ci.high = ci.tmp_yrtnormhigh[3,], ci.group = "highPrecipNovDecJanFebMar")
ci.tmp_yrPrecipNovDecJanFebMarmid.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrtnormmid[2,], ci.low = ci.tmp_yrtnormmid[1,], ci.high = ci.tmp_yrtnormmid[3,], ci.group = "midPrecipNovDecJanFebMar")
ci.tmp_yrPrecipNovDecJanFebMarlow.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrtnormlow[2,], ci.low = ci.tmp_yrtnormlow[1,], ci.high = ci.tmp_yrtnormlow[3,], ci.group = "lowPrecipNovDecJanFebMar")
tmp_yr_PrecipNovDecJanFebMarint <- rbind(ci.tmp_yrPrecipNovDecJanFebMarhigh.df, ci.tmp_yrPrecipNovDecJanFebMarmid.df, ci.tmp_yrPrecipNovDecJanFebMarlow.df)
ggplot(data = tmp_yr_PrecipNovDecJanFebMarint, aes(x = tmp_yr, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_yr and size interaction
tmp_yrrng <- range(grow_train$tmp_yr,na.rm = TRUE) #setting range for tmp_normrng
tmp_yr <- seq(tmp_yrrng[1], tmp_yrrng[2], by = 0.50)
size <- mean(grow_train$DIA_prev)
ppt_norm <- mean(grow_train$ppt_norm)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
size_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictiontmpyr_highsize <- growthpredictiontmpyr_lowsize <- growthpredictiontmpyr_midsize <- matrix(NA, length(plotdatainterval$u_beta_tmp_yr), length(tmp_yr)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_yr)){
  growthpredictiontmpyr_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[2]+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size_range[2] + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[2] + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size_range[2] +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size_range[2]
  
  growthpredictiontmpyr_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictiontmpyr_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[1]+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size_range[1] + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[1] + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size_range[1] +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size_range[1]
}
tmp_yr_prediction_trsizelow <- exp(growthpredictiontmpyr_lowsize)
tmp_yr_prediction_trsizemid <- exp(growthpredictiontmpyr_midsize)
tmp_yr_prediction_trsizehigh <- exp(growthpredictiontmpyr_highsize)
ci.tmp_yrsizehigh <- apply(tmp_yr_prediction_trsizehigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_yrsizemid <- apply(tmp_yr_prediction_trsizemid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_yrsizelow <- apply(tmp_yr_prediction_trsizelow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_yrsizehigh.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrsizehigh[2,], ci.low = ci.tmp_yrsizehigh[1,], ci.high = ci.tmp_yrsizehigh[3,], ci.group = "highsize")
ci.tmp_yrsizemid.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrsizemid[2,], ci.low = ci.tmp_yrsizemid[1,], ci.high = ci.tmp_yrsizemid[3,], ci.group = "midsize")
ci.tmp_yrsizelow.df <- data.frame(tmp_yr = tmp_yr, median = ci.tmp_yrsizelow[2,], ci.low = ci.tmp_yrsizelow[1,], ci.high = ci.tmp_yrsizelow[3,], ci.group = "lowsize")
tmp_yr_sizeint <- rbind(ci.tmp_yrsizehigh.df, ci.tmp_yrsizemid.df, ci.tmp_yrsizelow.df)
ggplot(data = tmp_yr_sizeint, aes(x = tmp_yr, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#ppt_norm and tmp_norm interaction
ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) 
ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
size <- mean(grow_train$DIA_prev)
tmp_yr <- mean(grow_train$tmp_yr)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictionpnorm_hightnorm <- growthpredictionpnorm_lowtnorm <- growthpredictionpnorm_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 

for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
  growthpredictionpnorm_hightnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[2] +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm_range[2]*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[2]*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm_range[2] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictionpnorm_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictionpnorm_lowtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm_range[1] +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm_range[1]*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm_range[1]*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm_range[1] + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
}
pnorm_prediction_trtnormlow <- exp(growthpredictionpnorm_lowtnorm)
pnorm_prediction_trtnormmid <- exp(growthpredictionpnorm_midtnorm)
pnorm_prediction_trtnormhigh <- exp(growthpredictionpnorm_hightnorm)
ci.pnorm_tnormhigh <- apply(pnorm_prediction_trtnormhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.pnorm_tnormmid <- apply(pnorm_prediction_trtnormmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.pnorm_tnormlow <- apply(pnorm_prediction_trtnormlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.pnorm_tnormhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_tnormhigh[2,], ci.low = ci.pnorm_tnormhigh[1,], ci.high = ci.pnorm_tnormhigh[3,], ci.group = "hightmp_norm")
ci.pnorm_tnormmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_tnormmid[2,], ci.low = ci.pnorm_tnormmid[1,], ci.high = ci.pnorm_tnormmid[3,], ci.group = "midtmp_norm")
ci.pnorm_tnormlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_tnormlow[2,], ci.low = ci.pnorm_tnormlow[1,], ci.high = ci.pnorm_tnormlow[3,], ci.group = "lowtmp_norm")
pnorm_tnormint <- rbind(ci.pnorm_tnormhigh.df, ci.pnorm_tnormmid.df, ci.pnorm_tnormlow.df)
ggplot(data = pnorm_tnormint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#ppt_norm and size interaction
ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) 
ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
size <- mean(grow_train$DIA_prev)
tmp_yr <- mean(grow_train$tmp_yr)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
size_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictionpnorm_highsize <- growthpredictionpnorm_lowsize <- growthpredictionpnorm_midsize <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 

for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
  growthpredictionpnorm_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[2]+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size_range[2] + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[2] + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size_range[2] +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size_range[2]
  
  growthpredictionpnorm_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictionpnorm_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[1]+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size_range[1] + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[1] + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size_range[1] +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size_range[1]
}
pnorm_prediction_trsizelow <- exp(growthpredictionpnorm_lowsize)
pnorm_prediction_trsizemid <- exp(growthpredictionpnorm_midsize)
pnorm_prediction_trsizehigh <- exp(growthpredictionpnorm_highsize)
ci.pnorm_sizehigh <- apply(pnorm_prediction_trsizehigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.pnorm_sizemid <- apply(pnorm_prediction_trsizemid, 2, quantile, c(0.025, 0.5, 0.975))
ci.pnorm_sizelow <- apply(pnorm_prediction_trsizelow, 2, quantile, c(0.025, 0.5, 0.975))
ci.pnorm_sizehigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_sizehigh[2,], ci.low = ci.pnorm_sizehigh[1,], ci.high = ci.pnorm_sizehigh[3,], ci.group = "highsize")
ci.pnorm_sizemid.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_sizemid[2,], ci.low = ci.pnorm_sizemid[1,], ci.high = ci.pnorm_sizemid[3,], ci.group = "midsize")
ci.pnorm_sizelow.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_sizelow[2,], ci.low = ci.pnorm_sizelow[1,], ci.high = ci.pnorm_sizelow[3,], ci.group = "lowsize")
pnorm_sizeint <- rbind(ci.pnorm_sizehigh.df, ci.pnorm_sizemid.df, ci.pnorm_sizelow.df)
ggplot(data = pnorm_sizeint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#ppt_norm and Precip_NovDecJanFebMar
ppt_normrng <- range(grow_train$ppt_norm,na.rm = TRUE) 
ppt_norm <- seq(ppt_normrng[1], ppt_normrng[2], by = 0.50)
size <- mean(grow_train$DIA_prev)
tmp_yr <- mean(grow_train$tmp_yr)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Precip_NovDecJanFebMar_range <- quantile(grow_train$Precip_NovDecJanFebMar, c(0.2, 0.8))
growthpredictionpnorm_highPrecipNovDecJanFebMar <- growthpredictionpnorm_lowPrecipNovDecJanFebMar <- growthpredictionpnorm_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_ppt_norm), length(ppt_norm)) 

for(i in 1:length(plotdatainterval$u_beta_ppt_norm)){
  growthpredictionpnorm_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[2] +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar_range[2]*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar_range[2]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar_range[2]*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar_range[2]*size
  
  growthpredictionpnorm_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictionpnorm_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[1] +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar_range[1]*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar_range[1]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar_range[1]*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar_range[1]*size
}
pnorm_prediction_trPrecipNovDecJanFebMarlow <- exp(growthpredictionpnorm_lowPrecipNovDecJanFebMar)
pnorm_prediction_trPrecipNovDecJanFebMarmid <- exp(growthpredictionpnorm_midPrecipNovDecJanFebMar)
pnorm_prediction_trPrecipNovDecJanFebMarhigh <- exp(growthpredictionpnorm_highPrecipNovDecJanFebMar)
ci.pnorm_PrecipNovDecJanFebMarhigh <- apply(pnorm_prediction_trPrecipNovDecJanFebMarhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.pnorm_PrecipNovDecJanFebMarmid <- apply(pnorm_prediction_trPrecipNovDecJanFebMarmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.pnorm_PrecipNovDecJanFebMarlow <- apply(pnorm_prediction_trPrecipNovDecJanFebMarlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.pnorm_PrecipNovDecJanFebMarhigh.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_PrecipNovDecJanFebMarhigh[2,], ci.low = ci.pnorm_PrecipNovDecJanFebMarhigh[1,], ci.high = ci.pnorm_PrecipNovDecJanFebMarhigh[3,], ci.group = "highPrecipNovDecJanFebMar")
ci.pnorm_PrecipNovDecJanFebMarmid.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_PrecipNovDecJanFebMarmid[2,], ci.low = ci.pnorm_PrecipNovDecJanFebMarmid[1,], ci.high = ci.pnorm_PrecipNovDecJanFebMarmid[3,], ci.group = "midPrecipNovDecJanFebMar")
ci.pnorm_PrecipNovDecJanFebMarlow.df <- data.frame(ppt_norm = ppt_norm, median = ci.pnorm_PrecipNovDecJanFebMarlow[2,], ci.low = ci.pnorm_PrecipNovDecJanFebMarlow[1,], ci.high = ci.pnorm_PrecipNovDecJanFebMarlow[3,], ci.group = "lowPrecipNovDecJanFebMar")
pnorm_PrecipNovDecJanFebMarint <- rbind(ci.pnorm_PrecipNovDecJanFebMarhigh.df, ci.pnorm_PrecipNovDecJanFebMarmid.df, ci.pnorm_PrecipNovDecJanFebMarlow.df)
ggplot(data = pnorm_PrecipNovDecJanFebMarint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#tmp_norm and Precip_NovDecJanFebMar
tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) 
tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
size <- mean(grow_train$DIA_prev)
tmp_yr <- mean(grow_train$tmp_yr)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Precip_NovDecJanFebMar_range <- quantile(grow_train$Precip_NovDecJanFebMar, c(0.2, 0.8))
growthpredictiontnorm_highPrecipNovDecJanFebMar <- growthpredictiontnorm_lowPrecipNovDecJanFebMar <- growthpredictiontnorm_midPrecipNovDecJanFebMar <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
  growthpredictiontnorm_highPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[2] +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar_range[2]*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar_range[2]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar_range[2]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar_range[2]*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar_range[2]*size
  
  growthpredictiontnorm_midPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictiontnorm_lowPrecipNovDecJanFebMar[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar_range[1] +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar_range[1]*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar_range[1]*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar_range[1]*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar_range[1]*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar_range[1]*size
}
tnorm_prediction_trPrecipNovDecJanFebMarlow <- exp(growthpredictiontnorm_lowPrecipNovDecJanFebMar)
tnorm_prediction_trPrecipNovDecJanFebMarmid <- exp(growthpredictiontnorm_midPrecipNovDecJanFebMar)
tnorm_prediction_trPrecipNovDecJanFebMarhigh <- exp(growthpredictiontnorm_highPrecipNovDecJanFebMar)
ci.tnorm_PrecipNovDecJanFebMarhigh <- apply(tnorm_prediction_trPrecipNovDecJanFebMarhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tnorm_PrecipNovDecJanFebMarmid <- apply(tnorm_prediction_trPrecipNovDecJanFebMarmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tnorm_PrecipNovDecJanFebMarlow <- apply(tnorm_prediction_trPrecipNovDecJanFebMarlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tnorm_PrecipNovDecJanFebMarhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tnorm_PrecipNovDecJanFebMarhigh[2,], ci.low = ci.tnorm_PrecipNovDecJanFebMarhigh[1,], ci.high = ci.tnorm_PrecipNovDecJanFebMarhigh[3,], ci.group = "highPrecipNovDecJanFebMar")
ci.tnorm_PrecipNovDecJanFebMarmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tnorm_PrecipNovDecJanFebMarmid[2,], ci.low = ci.tnorm_PrecipNovDecJanFebMarmid[1,], ci.high = ci.tnorm_PrecipNovDecJanFebMarmid[3,], ci.group = "midPrecipNovDecJanFebMar")
ci.tnorm_PrecipNovDecJanFebMarlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tnorm_PrecipNovDecJanFebMarlow[2,], ci.low = ci.tnorm_PrecipNovDecJanFebMarlow[1,], ci.high = ci.tnorm_PrecipNovDecJanFebMarlow[3,], ci.group = "lowPrecipNovDecJanFebMar")
tnorm_PrecipNovDecJanFebMarint <- rbind(ci.tnorm_PrecipNovDecJanFebMarhigh.df, ci.tnorm_PrecipNovDecJanFebMarmid.df, ci.tnorm_PrecipNovDecJanFebMarlow.df)
ggplot(data = tnorm_PrecipNovDecJanFebMarint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#tmp_norm and size interaction
tmp_normrng <- range(grow_train$tmp_norm,na.rm = TRUE) #setting range for tmp_normrng
tmp_norm <- seq(tmp_normrng[1], tmp_normrng[2], by = 0.50)
size <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar
ppt_norm <- mean(grow_train$ppt_norm)
tmp_yr <- mean(grow_train$tmp_yr)
size_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictiontnorm_highsize <- growthpredictiontnorm_lowsize <- growthpredictiontnorm_midsize <- matrix(NA, length(plotdatainterval$u_beta_tmp_norm), length(tmp_norm)) 

for(i in 1:length(plotdatainterval$u_beta_tmp_norm)){
  growthpredictiontnorm_highsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[2]+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size_range[2] + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[2] + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size_range[2] +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size_range[2]
  
  growthpredictiontnorm_midsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  
  growthpredictiontnorm_lowsize[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[1]+
    plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size_range[1] + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[1] + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size_range[1] +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size_range[1]
}
tnorm_prediction_trlowsize <- exp(growthpredictiontnorm_lowsize)
tnorm_prediction_trmidsize <- exp(growthpredictiontnorm_midsize)
tnorm_prediction_trhighsize <- exp(growthpredictiontnorm_highsize)
ci.tnormhigh <- apply(tnorm_prediction_trhighsize, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tnormmid <- apply(tnorm_prediction_trmidsize, 2, quantile, c(0.025, 0.5, 0.975))
ci.tnormlow <- apply(tnorm_prediction_trlowsize, 2, quantile, c(0.025, 0.5, 0.975))
ci.tnormhigh.df <- data.frame(tmp_norm = tmp_norm, median = ci.tnormhigh[2,], ci.low = ci.tnormhigh[1,], ci.high = ci.tnormhigh[3,], ci.group = "highsize")
ci.tnormmid.df <- data.frame(tmp_norm = tmp_norm, median = ci.tnormmid[2,], ci.low = ci.tnormmid[1,], ci.high = ci.tnormmid[3,], ci.group = "midsize")
ci.tnormlow.df <- data.frame(tmp_norm = tmp_norm, median = ci.tnormlow[2,], ci.low = ci.tnormlow[1,], ci.high = ci.tnormlow[3,], ci.group = "lowsize")
tnorm_sizeint <- rbind(ci.tnormhigh.df, ci.tnormmid.df, ci.tnormlow.df)
ggplot(data = tnorm_sizeint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#plotting individual tree growth
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -111, "-111 to -109", "-114 to -111")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 35, "35 to 37", "32 to 35")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(16)

#tmp_yr and tmp_norm
get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  tmp_yrrng <- range(tree.grow$tmp_yr,na.rm = TRUE) #setting range for tmp_normrng
  tmp_yr <- seq(tmp_yrrng[1], tmp_yrrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictiontmpyr_hightnorm <- growthpredictiontmpyr_lowtnorm <- growthpredictiontmpyr_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_tmp_yr), length(tmp_yr)) 
  
  for(i in 1:length(plotdatainterval$u_beta_tmp_yr)){
    growthpredictiontmpyr_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
      plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_ppt_norm"]*Precip_NovDecJanFebMar*ppt_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_norm"]*Precip_NovDecJanFebMar*tmp_norm + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Precip_JulAug"]*Precip_NovDecJanFebMar*Precip_JulAug +
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_tmp_yr"]*Precip_NovDecJanFebMar*tmp_yr + 
      plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_DIA_prev"]*Precip_NovDecJanFebMar*size
  }
  tmp_yr_prediction_trtnormmid <- exp(growthpredictiontmpyr_midtnorm)
  ci.tmp_yrtnormmid <- apply(tmp_yr_prediction_trtnormmid, 2, quantile, c(0.025, 0.5, 0.975))
  ci.tmp_yrtnormmid.df <- data.frame(tmp_yr = tmp_yr, tmp_norm = tmp_norm, median = ci.tmp_yrtnormmid[2,], ci.low = ci.tmp_yrtnormmid[1,], ci.high = ci.tmp_yrtnormmid[3,], ci.group = tree.subset$treeCD)
  tmp_yr_tnormint <- rbind(ci.tmp_yrtnormmid.df)
  print(ind.samples[j,])
  tmp_yr_tnormint  
}
#get.ind.tmp.response(i = 6)
tmp_yr_tree_response <- list()
tmp_yr_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
tmp_yr_tree_response.df <- do.call(rbind, tmp_yr_tree_response)
merged.response.samples <- merge(tmp_yr_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = tmp_yr, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
 #color by LATLONbin
ggplot(data = merged.response.samples, aes(x = tmp_yr, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = tmp_yr, y = median, color = tmp_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#map of LATLONbin
all_states <- map_data("state")
states <- subset(all_states, region %in% c("arizona"))
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata<-states
mapdata<-data.frame(mapdata)
ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
  geom_point(data = grow.monsoon, aes(x = LON, y = LAT, color = LATLONbin)) + theme_bw()




#Precip_JulAug and tmp_yr
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -111, "-111 to -109", "-114 to -111")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 35, "35 to 37", "32 to 35")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for Precip_JulAugrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  tmp_yr <- mean(tree.grow$tmp_yr)
  tmp_yr_range <- quantile(tree.grow$tmp_yr, c(0.2, 0.8))
  growthpredictionPrecipJulAug_midtmpyr <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_midtmpyr[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size
  }
  PrecipJulAug_prediction_trtmpyr <- exp(growthpredictionPrecipJulAug_midtmpyr)
  ci.PrecipJulAug_tmpyrmid <- apply(PrecipJulAug_prediction_trtmpyr, 2, quantile, c(0.025, 0.5, 0.975))
  ci.PrecipJulAug_tmpyrmid.df <- data.frame(Precip_JulAug = Precip_JulAug, tmp_yr = tmp_yr, median = ci.PrecipJulAug_tmpyrmid[2,], ci.low = ci.PrecipJulAug_tmpyrmid[1,], ci.high = ci.PrecipJulAug_tmpyrmid[3,], ci.group = tree.subset$treeCD)
  PrecipJulAug_tmpyrint <- rbind(ci.PrecipJulAug_tmpyrmid.df)
  print(ind.samples[j,])
  PrecipJulAug_tmpyrint  
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
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = tmp_yr, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)

  

#Precip_JulAug and tmp_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -111, "-111 to -109", "-114 to -111")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 35, "35 to 37", "32 to 35")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for Precip_JulAugrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  tmp_yr <- mean(tree.grow$tmp_yr)
  tmp_norm_range <- quantile(tree.grow$tmp_norm, c(0.2, 0.8))
  growthpredictionPrecipJulAug_midtnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_midtnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size
  }
  PrecipJulAug_prediction_trtnorm <- exp(growthpredictionPrecipJulAug_midtnorm)
  ci.PrecipJulAug_tnormmid <- apply(PrecipJulAug_prediction_trtnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.PrecipJulAug_tnormmid.df <- data.frame(Precip_JulAug = Precip_JulAug, tmp_norm = tmp_norm, median = ci.PrecipJulAug_tnormmid[2,], ci.low = ci.PrecipJulAug_tnormmid[1,], ci.high = ci.PrecipJulAug_tnormmid[3,], ci.group = tree.subset$treeCD)
  PrecipJulAug_tnormint <- rbind(ci.PrecipJulAug_tnormmid.df)
  print(ind.samples[j,])
  PrecipJulAug_tnormint  
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



#ppt_norm and tmp_yr
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -111, "-111 to -109", "-114 to -111")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 35, "35 to 37", "32 to 35")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  tmp_yrrng <- range(tree.grow$tmp_yr,na.rm = TRUE) #setting range for Precip_JulAugrng
  tmp_yr <- seq(tmp_yrrng[1], tmp_yrrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  Precip_JulAug <- mean(tree.grow$Precip_JulAug)
  tmp_norm <- mean(tree.grow$tmp_norm)
  ppt_norm <- mean(tree.grow$ppt_norm)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionpnorm_tmpyr <- matrix(NA, length(plotdatainterval$u_beta_tmp_yr), length(tmp_yr)) 
  
  for(i in 1:length(plotdatainterval$u_beta_tmp_yr)){
    growthpredictionpnorm_tmpyr[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size
  }
  pnorm_prediction_trtmpyr <- exp(growthpredictionpnorm_tmpyr)
  ci.pnorm_tmpyr <- apply(pnorm_prediction_trtmpyr, 2, quantile, c(0.025, 0.5, 0.975))
  ci.pnorm_tmpyr.df <- data.frame(ppt_norm = ppt_norm, tmp_yr = tmp_yr, median = ci.pnorm_tmpyr[2,], ci.low = ci.pnorm_tmpyr[1,], ci.high = ci.pnorm_tmpyr[3,], ci.group = tree.subset$treeCD)
  pnorm_tmpyrint <- rbind(ci.pnorm_tmpyr.df)
  print(ind.samples[j,])
  pnorm_tmpyrint  
}
#get.ind.tmp.response(i = 6)
ppt_norm_tree_response <- list()
ppt_norm_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
ppt_norm_tree_response.df <- do.call(rbind, ppt_norm_tree_response)
merged.response.samples <- merge(ppt_norm_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = tmp_yr, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = tmp_yr, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_yr
ggplot(data = merged.response.samples, aes(x = tmp_yr, y = median, color = ppt_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Precip_JulAug and ppt_norm
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -111, "-111 to -109", "-114 to -111")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 35, "35 to 37", "32 to 35")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin) %>% sample_n(7)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Precip_JulAugrng <- range(tree.grow$Precip_JulAug,na.rm = TRUE) #setting range for Precip_JulAugrng
  Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.1)
  size <- mean(tree.grow$DIA_prev)
  ppt_norm <- mean(tree.grow$ppt_norm)
  tmp_norm <- mean(tree.grow$tmp_norm)
  tmp_yr <- mean(tree.grow$tmp_yr)
  ppt_norm_range <- quantile(tree.grow$ppt_norm, c(0.2, 0.8))
  growthpredictionPrecipJulAug_pnorm <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 
  
  for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
    growthpredictionPrecipJulAug_pnorm[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_tmp_yr"]*tmp_yr +
      plotdatainterval[i,"u_beta_DIA_prev"]*size + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
      plotdatainterval[i,"u_beta_ppt_norm_tmp_yr"]*ppt_norm*tmp_yr + plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug +
      plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size+
      plotdatainterval[i,"u_beta_Precip_JulAug_tmp_yr"]*Precip_JulAug*tmp_yr + plotdatainterval[i,"u_beta_Precip_JulAug_tmp_norm"]*Precip_JulAug*tmp_norm +
      plotdatainterval[i,"u_beta_Precip_JulAug_DIA_prev"]*Precip_JulAug*size + plotdatainterval[i,"u_beta_tmp_norm_tmp_yr"]*tmp_norm*tmp_yr + 
      plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size + plotdatainterval[i,"u_beta_tmp_yr_DIA_prev"]*tmp_yr*size
  }
  PrecipJulAug_prediction_trpnorm <- exp(growthpredictionPrecipJulAug_pnorm)
  ci.PrecipJulAug_pnorm <- apply(PrecipJulAug_prediction_trpnorm, 2, quantile, c(0.025, 0.5, 0.975))
  ci.PrecipJulAug_pnorm.df <- data.frame(Precip_JulAug = Precip_JulAug, ppt_norm = ppt_norm, median = ci.PrecipJulAug_pnorm[2,], ci.low = ci.PrecipJulAug_pnorm[1,], ci.high = ci.PrecipJulAug_pnorm[3,], ci.group = tree.subset$treeCD)
  PrecipJulAug_pnormint <- rbind(ci.PrecipJulAug_pnorm.df)
  print(ind.samples[j,])
  PrecipJulAug_pnormint  
}
#get.ind.tmp.response(i = 6)
PrecipJulAug_tree_response <- list()
PrecipJulAug_tree_response <- lapply(1:length(ind.samples$treeCD), FUN = get.ind.tmp.response)
PrecipJulAug_tree_response.df <- do.call(rbind, PrecipJulAug_tree_response)
merged.response.samples <- merge(PrecipJulAug_tree_response.df, ind.samples, by.x = "ci.group", by.y = "treeCD")
merged.response.samples$ci.group <- as.character(merged.response.samples$ci.group)
#color by group
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by LATLONbin
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = LATLONbin, group = ci.group)) + geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = LATLONbin, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
#color by tmp_norm
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = ppt_norm, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)








#ggplot map of tmp_norm
library(hexbin)
all_states <- map_data("state")
states <- subset(all_states, region %in% c("arizona"))
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata<-states
mapdata<-data.frame(mapdata)
ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
  geom_point(data = grow, aes(x = LON, y = LAT, color = tmp_norm))+ scale_color_gradient2(low = "blue", mid = "white", high = "red")+
  theme_bw()


#ggplot map of ppt_norm
all_states <- map_data("state")
states <- subset(all_states, region %in% c("arizona"))
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata<-states
mapdata<-data.frame(mapdata)
ggplot() + geom_polygon(data=mapdata, aes(x=long, y=lat, group = group), color ="darkgray", fill = "darkgray")+
  geom_point(data = grow, aes(x = LON, y = LAT, color = ppt_norm))+ scale_color_gradient2(low = "red", mid = "white", high = "blue")+
  theme_bw()




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
