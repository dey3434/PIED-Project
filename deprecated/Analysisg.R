library(tidyverse)
library(lme4)
library(bbmle)

setwd("~/Research/FIAcores")

climdata<-read.csv("data.clim.csv")

climdata2 <- subset(climdata, climdata$growth < 10)
climdata2$loggrowth <- log(climdata2$growth)
hist(climdata2$loggrowth)
hist(climdata2$growth)
hist(climdata2$ppt_norm)
hist(climdata2$ppt_yr)
hist(climdata2$tmp_norm)
hist(climdata2$tmp_yr)
hist(climdata2$solrad_an)

#visualization ggplot
growthxpptnorm <- ggplot(data = climdata2, aes(x = ppt_norm, y = loggrowth)) + 
  geom_point()
growthxpptnorm

growthxpptyr <- ggplot(data = climdata2, aes(x = ppt_yr, y = loggrowth)) + 
  geom_point()
growthxpptyr

growthxtmpnorm <- ggplot(data = climdata2, aes(x = tmp_norm, y = loggrowth)) + 
  geom_point()
growthxtmpnorm

growthxtmpyr <- ggplot(data = climdata2, aes(x = tmp_yr, y = loggrowth)) + 
  geom_point()
growthxtmpyr

growthxrad <- ggplot(data = climdata2, aes(x = solrad_an, y = loggrowth)) +
  geom_point()
growthxrad

growthxdiam <- ggplot(data = climdata2, aes(x = DIA, y = loggrowth)) +
  geom_point()
growthxdiam

#correlation of predictors
tmpnormxpptnorm <- ggplot(data = climdata2, aes(x = tmp_norm, y = ppt_norm)) +
  geom_point()
tmpnormxpptnorm

tmpyrxpptyr <- ggplot(data = climdata2, aes(x = tmp_yr, y = ppt_yr)) +
  geom_point()
tmpyrxpptyr

tmpnormxtmpyr <- ggplot(data = climdata2, aes(x = tmp_norm, y = tmp_yr)) +
  geom_point()
tmpnormxtmpyr

pptnormxpptyr <- ggplot(data = climdata2, aes(x = ppt_norm, y = ppt_yr)) +
  geom_point()
pptnormxpptyr

radxpptyr <- ggplot(data = climdata2, aes(x = solrad_an, y = ppt_yr)) +
  geom_point()
radxpptyr

radxtmpyr <- ggplot(data = climdata2, aes(x = solrad_an, y = tmp_yr)) +
  geom_point()
radxtmpyr

#scale predictor variables
climdata.scaled <- na.omit(climdata2) %>% 
  mutate_at(scale, .vars = vars(ppt_norm,ppt_yr,tmp_norm,tmp_yr,growth,solrad_an))

hist(climdata.scaled$growth)
hist(climdata.scaled$solrad_an)
str(climdata.scaled)

#models
#nomodel
null <- lm(growth~1, climdata.scaled)
null
mean(climdata.scaled$growth)

#temperature
tempnorm <- lm(growth~tmp_norm, climdata.scaled)
tempnorm
summary(tempnorm)

tempyr <- lm(growth~tmp_yr, climdata.scaled)
tempyr

pptnorm <- lm(growth~ppt_norm, climdata.scaled)
pptnorm

pptyr <- lm(growth~ppt_yr, climdata.scaled)
pptyr

pptyrnorm <- lm(growth~ppt_yr+ppt_norm, climdata.scaled)
pptyrnorm

tempyrpptnorm <- lm(growth~tmp_yr+ppt_norm, climdata.scaled)
tempyrpptnorm

pptyrnormtempyr <- lm(growth~tmp_yr+ppt_norm+ppt_yr, climdata.scaled)
pptyrnormtempyr

radpptyrtempyr <- lm(growth~tmp_yr+solrad_an+ppt_yr, climdata.scaled)
radpptyrtempyr

#interactions
interannual <- lm(growth~tmp_yr*ppt_yr+ppt_norm, climdata.scaled)
interannual

interall <- lm(growth~tmp_yr*ppt_yr*ppt_norm, climdata.scaled)
interall

interannual_rad <- lm(growth~tmp_yr*ppt_yr+solrad_an, climdata.scaled)
interannual_rad

interall_rad <- lm(growth~tmp_yr*ppt_yr*solrad_an, climdata.scaled)
interall_rad

#compare models
AICtab(null, tempnorm, pptnorm, tempyr, pptyr)
AICtab(pptnorm, pptyrnorm, tempyrpptnorm, pptyrnormtempyr,radpptyrtempyr)
AICtab(pptyrnormtempyr, interannual, interall, interannual_rad, interall_rad)

#Plots
mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               legend.text=element_text(size=11),legend.title=element_text(size=12),
               legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
               axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
               axis.line.x = element_line(color="black", size = 0.3),
               axis.line.y = element_line(color="black", size = 0.3))

pfun_interall<-function(x,temp,norm){
  interall$coef[1] + interall$coef[2]*temp +
    interall$coef[3]*x + interall$coef[4]*norm +
    interall$coef[5]*temp*x + interall$coef[6]*temp*norm +
    interall$coef[7]*x*norm + interall$coef[8]*temp*x*norm
}

tfun_interall<-function(x,ppt,norm){
  interall$coef[1] + interall$coef[2]*x +
    interall$coef[3]*ppt + interall$coef[4]*norm +
    interall$coef[5]*x*ppt + interall$coef[6]*x*norm +
    interall$coef[7]*ppt*norm + interall$coef[8]*ppt*x*norm
}

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
