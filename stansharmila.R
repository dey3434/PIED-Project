### Stan models for PIAL growth using technique in Ogle et al., Ecology Letters to test for climate effects
## Sharmila Dey
# 22 June 2020
setwd("~/Desktop/tree ring lab")
load("./pied_grow_coef2.rda")
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

data<-read.csv("data.clim copy.csv") #Bring in ring-widths and climate data

#make data object for stan
grow<-na.omit(data) %>% 
  mutate_at(scale, .vars = vars(ppt_norm,ppt_yr,tmp_norm,tmp_yr,solrad_an)) %>%  #standardizing covariates
  arrange(PLOT.x,SUBP,name) %>% #ordering by plot, subplot, tree
  mutate(PlotCD=as.numeric(factor(PLOT.x, levels = unique(PLOT.x))),treeCD=as.numeric(factor(name,levels=unique(name))),
         growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2))   #assigning plot and tree indicies, log transformation of ring-widths

split=0.20  #setting aside testing data
trainIndex <- createDataPartition(grow$name, p=split, list=FALSE)  #20% of trees
grow_test <- grow[trainIndex,]
grow_train <- grow[-trainIndex,]


#create training design matrix
xG<-as.matrix(cbind(grow_train$ppt_norm, grow_train$tmp_norm, grow_train$ppt_yr, grow_train$tmp_yr, grow_train$DIA_prev,
                    grow_train$ppt_norm*grow_train$tmp_norm, grow_train$ppt_norm*grow_train$tmp_yr, 
                    grow_train$ppt_norm*grow_train$ppt_yr, grow_train$ppt_norm*grow_train$DIA_prev,
                    grow_train$ppt_yr*grow_train$tmp_yr, grow_train$ppt_yr*grow_train$tmp_norm, grow_train$ppt_yr*grow_train$DIA_prev, 
                    grow_train$tmp_norm*grow_train$tmp_yr, grow_train$tmp_norm*grow_train$DIA_prev, 
                    grow_train$tmp_yr*grow_train$DIA_prev))
#create testing design matrix
xGtest<-as.matrix(cbind(grow_test$ppt_norm, grow_test$tmp_norm, grow_test$ppt_yr, grow_test$tmp_yr, grow_test$DIA_prev,
                                  grow_test$ppt_norm*grow_test$tmp_norm, grow_test$ppt_norm*grow_test$tmp_yr, 
                                  grow_test$ppt_norm*grow_test$ppt_yr, grow_test$ppt_norm*grow_test$DIA_prev,
                                  grow_test$ppt_yr*grow_test$tmp_yr, grow_test$ppt_yr*grow_test$tmp_norm, grow_test$ppt_yr*grow_test$DIA_prev, 
                                  grow_test$tmp_norm*grow_test$tmp_yr, grow_test$tmp_norm*grow_test$DIA_prev, 
                                  grow_test$tmp_yr*grow_test$DIA_prev))
yG<-as.vector(grow_train$loggrowth)   #create training response vector
yGtest<-as.vector(grow_test$loggrowth)   #create testing response vector
nG<-nrow(grow_train)   #number of training data points
nGtest<-nrow(grow_test)   #number of testing data points
plot<-grow_train$PlotCD   #create plot identifier
nplot<-length(unique(grow_train$PlotCD))   #how many plots
K<-ncol(xG)   #number of predictors
tree<-grow_train$treeCD   #create tree identifier in train dataset
treetest<-grow_test$treeCD   #create tree identifier in test dataset
ntree<-length(unique(grow_train$treeCD))   #number of trees in training dataset
plotfortree<-grow_train %>%   #plot code for each tree
  group_by(treeCD) %>%
  summarize(Plot=mean(PlotCD))
plotfortree<-plotfortree$Plot   #create vector of plot code


sink("pied_grow.stan")   #create stan file
#put code in stan file
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
    
    //for(p in 1:nplot){
    //beta0_p[p] ~ normal(u_beta0, s_beta0_p);
    //  for(t in 1:ntree){
    //  beta0_t[t] ~ normal(beta0_p[plotfortree[t]], s_beta0_t);
    //  }
    //}
    
    
    // GROWTH MODEL
    
    for(n in 1:nG){
    mG[n] = beta0_t[tree[n]]+xG[n]*u_beta;
    }
    
    //yG ~ normal(mG,sigma_y);
    yG ~ gamma(mG,sigma_y);

    }
    
    generated quantities{
    vector[nGtest] yrep;
    //for(n in 1:nGtest){
    //yrep[n] = normal_rng(beta0_t[treetest[n]]+xGtest[n]*u_beta,sigma_y);
    //}

    for(n in 1:nGtest){
    yrep[n] = gamma_rng(beta0_t[treetest[n]]+xGtest[n]*u_beta,sigma_y);
    }

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

summary<-summary(fit_grow)
summary

plotdata<-select(as.data.frame(fit_grow),"yrep[1]":"yrep[8780]")
plotdatainterval<-select(as.data.frame(fit_grow), "u_beta[1]":"u_beta[15]")
colnames(plotdatainterval) <- c("u_beta_ppt_norm", "u_beta_tmp_norm", "u_beta_ppt_yr", "u_beta_tmp_yr", "u_beta_DIA_prev",
 "u_beta_ppt_norm_tmp_norm", "u_beta_ppt_norm_tmp_yr", "u_beta_ppt_norm_ppt_yr", "u_beta_ppt_norm_DIA_prev",
 "u_beta_ppt_yr_tmp_yr", "u_beta_ppt_yr_tmp_norm", "u_beta_ppt_yr_DIA_prev", "u_beta_tmp_norm_tmp_yr", "u_beta_tmp_norm_DIA_prev", 
 "u_beta_tmp_yr_DIA_prev")
ppc_dens_overlay(yGtest, as.matrix(plotdata))

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

#MCMC Area Plot
color_scheme_set("purple")
mcmc_areas(plotdatainterval, prob = 0.8)

#MCMC Traces
color_scheme_set("mix-blue-red")
mcmc_trace(plotdatainterval,
           facet_args = list(nrow = 2))
#pars = c("alpha", "sigma")

mcmcplot(As.mcmc.list(fit_grow))

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
