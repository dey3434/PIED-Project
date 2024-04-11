#### Hi Sharmila, I hope this helps! I included code for running the logistic regressions; setting up a dataframe with binned values of a 
 # predictor (lambda, in my case) and occurence, as well as predicted occurrence from the models; and created a plot with points for the binned
 # values and a line with the model prediction. I did not include the code for setting up the original dataframe (FIA_lambda_noex) because I didn't
 # think it would be helpful, but the basic structure is simple: it could include any number of predictors (mine had four different precictions
 # for lambda) and a column with presence/absence data.
library(raster)

PApied <- raster("data/presenceAbsenceRasterNA.tif")
# getting a warning: Warning message:
#In showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj = prefer_proj) :
#  Discarded datum unknown in Proj4 definition
plot(PApied)

FIA_PApied <- as.data.frame(PApied,xy=TRUE)
names(FIA_PApied )<-c("lon","lat","PApied")

summary(FIA_PApied)

# save this as a data frame so that you could read it in without needing the raster package:
write.csv(FIA_PApied, "data/PresenceAbsernce_PIEDFIA.csv")

## Logistic models (lambda vs occurrence)
k=5
pa_c <- gam(PApied~s(lambda_c,k=k),family=binomial(link = cloglog),
            data=FIA_lambda_noex)
pa_ci <- gam(PApied~s(lambda_ci,k=k),family=binomial(link = cloglog),
             data=FIA_lambda)
pa_cc <- gam(PApied~s(lambda_cc,k=k),family=binomial(link = cloglog),
             data=FIA_lambda_noex)
pa_i <- gam(PApied~s(lambda_i),family=binomial(link = cloglog),
            data=FIA_lambda_noex)

## Set up data frame with binned values
ncuts=50 # number of cut to create bins

# divide into intervals based on number of cuts
chopsize_lam_c<-cut(FIA_lambda_noex$lambda_c,ncuts)
chopsize_lam_ci<-cut(FIA_lambda_noex$lambda_ci,ncuts)
chopsize_lam_cc<-cut(FIA_lambda_noex$lambda_cc,ncuts)
chopsize_lam_i<-cut(FIA_lambda_noex$lambda_i,ncuts)

count_binned_lam_c<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_c),length)) # calculate number of data points in each bin
lam_c_binned<-as.vector(sapply(split(FIA_lambda_noex$lambda_c,chopsize_lam_c),mean,na.rm=T)) # calculate mean lambda in each bin
pres_c_binned<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_c),mean,na.rm=T)) # calculate mean prob of occurrence in each bin

count_binned_lam_ci<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_ci),length))
lam_ci_binned<-as.vector(sapply(split(FIA_lambda_noex$lambda_ci,chopsize_lam_ci),mean,na.rm=T))
pres_ci_binned<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_ci),mean,na.rm=T))

count_binned_lam_cc<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_cc),length))
lam_cc_binned<-as.vector(sapply(split(FIA_lambda_noex$lambda_cc,chopsize_lam_cc),mean,na.rm=T))
pres_cc_binned<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_cc),mean,na.rm=T))

count_binned_lam_i<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_i),length))
lam_i_binned<-as.vector(sapply(split(FIA_lambda_noex$lambda_i,chopsize_lam_i),mean,na.rm=T))
pres_i_binned<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_i),mean,na.rm=T))

# combine into dataframe
pres_binned<-data.frame(count_lam=c(count_binned_lam_c,count_binned_lam_cc,count_binned_lam_i),
                        lam=c(lam_c_binned,lam_cc_binned,lam_i_binned),
                        pres=c(pres_c_binned,pres_cc_binned,pres_i_binned),
                        pred=c(invlogit(predict(pa_c,newdata=data.frame(lambda_c=lam_c_binned))),
                               invlogit(predict(pa_cc,newdata=data.frame(lambda_cc=lam_cc_binned))),
                               invlogit(predict(pa_i,newdata=data.frame(lambda_i=lam_i_binned)))),
                        model=c(rep("c",ncuts),rep("cc",ncuts),rep("i",ncuts)))

# Plot
pres_plot_c<-ggplot(data=subset(pres_binned,model=="c"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  annotate("label", x = 0.89, y = 0.5, 
           label = paste("Deviance=",round(pa_c$deviance,2),
                         "\nAIC =",round(pa_c$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = expression(paste(lambda)), y = "Probability of occurrence")
