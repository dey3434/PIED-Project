### Stan models for PIAL growth using technique in Ogle et al., Ecology Letters to test for climate effects
## Sharmila Dey
# 22 June 2020
setwd("~/Desktop/tree ring lab")
library(rstan)
options(mc.cores = parallel::detectCores())
library(parallel) 
library(mcmcplots) ; library(lattice) ; library(MASS)
library(lme4) ; library(nlme) ; library(splines); library(MCMCpack)
library(ggplot2)
library(caret) ; library(tidyverse)
library(bayesplot)

data<-read.csv("data.clim copy.csv")

grow<-na.omit(data) %>% 
  mutate_at(scale, .vars = vars(ppt_norm,ppt_yr,tmp_norm,tmp_yr,growth,solrad_an)) %>%
  arrange(PLOT.x,SUBP,name) %>%
  mutate(PlotCD=as.numeric(factor(PLOT.x, levels = unique(PLOT.x))),treeCD=as.numeric(factor(name,levels=unique(name))))

split=0.20
trainIndex <- createDataPartition(grow$name, p=split, list=FALSE)
grow_test <- grow[trainIndex,]
grow_train <- grow[-trainIndex,]



xG<-as.matrix(cbind(grow_train$ppt_norm, grow_train$tmp_norm, grow_train$ppt_yr, grow_train$tmp_yr, grow_train$DIA_prev,
                    grow_train$ppt_norm*grow_train$tmp_norm, grow_train$ppt_norm*grow_train$tmp_yr, 
                    grow_train$ppt_norm*grow_train$ppt_yr, grow_train$ppt_norm*grow_train$DIA_prev,
                    grow_train$ppt_yr*grow_train$tmp_yr, grow_train$ppt_yr*grow_train$tmp_norm, grow_train$ppt_yr*grow_train$DIA_prev, 
                    grow_train$tmp_norm*grow_train$tmp_yr, grow_train$tmp_norm*grow_train$DIA_prev, 
                    grow_train$tmp_yr*grow_train$DIA_prev))
xGtest<-as.matrix(cbind(grow_test$ppt_norm, grow_test$tmp_norm, grow_test$ppt_yr, grow_test$tmp_yr, grow_test$DIA_prev,
                                  grow_test$ppt_norm*grow_test$tmp_norm, grow_test$ppt_norm*grow_test$tmp_yr, 
                                  grow_test$ppt_norm*grow_test$ppt_yr, grow_test$ppt_norm*grow_test$DIA_prev,
                                  grow_test$ppt_yr*grow_test$tmp_yr, grow_test$ppt_yr*grow_test$tmp_norm, grow_test$ppt_yr*grow_test$DIA_prev, 
                                  grow_test$tmp_norm*grow_test$tmp_yr, grow_test$tmp_norm*grow_test$DIA_prev, 
                                  grow_test$tmp_yr*grow_test$DIA_prev))
yG<-as.vector(grow_train$growth)
yGtest<-as.vector(grow_test$growth)
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
    
    yG ~ normal(mG,sigma_y);
    
    }
    
    generated quantities{
    vector[nGtest] yrep;
    for(n in 1:nGtest){
    yrep[n] = normal_rng(beta0_t[treetest[n]]+xGtest[n]*u_beta,sigma_y);
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

plotdata<-as.data.frame(fit_grow)
ppc_dens_overlay(yG, as.matrix(select(plotdata,"yrep[1]":"yrep[8780]")))

mcmcplot(As.mcmc.list(fit_grow))

write.csv(summary$summary,"./piedsummarylong.csv")

pied_grow_coef<-summary$summary
pied_grow_coef
save(pied_grow_coef, grow, file="./pied_grow_coef.rda")

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