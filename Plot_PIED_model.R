#----------------------------------------------------------------------------------------------
#  Code to make plots of posterior estimates from PIED stan models
#----------------------------------------------------------------------------------------------
#MCMC PLOTS
#MCMC Intervals Plots
mcmc_intervals(plotdatainterval, prob = 0.5, )
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
color_scheme_set("viridis")
mcmc_trace(plotdatainterval,
           facet_args = list(nrow = 2), 
           pars = c("u_beta_ppt_norm", "u_beta_tmp_norm"))
dev.off() # "device off" tells R to stop printing stuff to the pdf

mcmcplot(As.mcmc.list(fit_grow))


#PAIRS PANELS PLOTS
library(psych)

#raw data
xGpanels <- data.frame(xG[,1:7])
colnames(xGpanels) <- c("ppt_norm", "tmp_norm", "Precip_JulAug", "Precip_NovDecJanFebMar", "Tmean_AprMayJun", "Tmean_SepOct", "DIA_prev")
png("pairs_panels.png", width = 8, height = 8, units = "in", res = 200) 
pairs.panels(xGpanels)
dev.off()

#posterior data
colnames(plotdataintervalpanels) <- c("ppt_norm", "tmp_norm", "Precip_JA", "Precip_NDJFM", "Tmean_AMJ", "Tmean_SO", "size",
                                      "pn_tn", "pn_PJA", "pn_size", "pn_PNDJFM", "pn_TAMJ", "pn_TSO", "tn_size", "tn_PJA",
                                      "tn_PNDJFM", "tn_TAMJ", "tn_TSO", "size_PJA", "size_PNDJFM", "size_TAMJ",
                                      "size_TSO", "PJA_PNDJFM", "PJA_TAMJ", "PJA_TSO", "PNDJFM_TAMJ", "PNDJFM_TSO", "TAMJ_TSO")
png("ppc_pairs_panels.png", width = 15, height = 15, units = "in", res = 200) 
pairs.panels(plotdatainterval)
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


#INTERACTION PLOTS
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
ggplot(data = Precip_JulAug_Precip_NovDecJanFebMarint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
  scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
  geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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
ggplot(data = Precip_JulAug_tnormint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
  scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
  geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



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
ggplot(data = Precip_JulAug_pnormint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
  scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
  geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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
ggplot(data = Precip_JulAug_sizeint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
  scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
  geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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
ggplot(data = ppt_norm_pnormint, aes(x = ppt_norm, y = median, color = ci.group)) +
  scale_color_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) +
  geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#B2182B", "#FDDBC7", "#4575B4")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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
ci.Tmean_AprMayJunhigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Tmean_AprMayJunmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_AprMayJunhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunhigh[2,], ci.low = ci.Tmean_AprMayJunhigh[1,], ci.high = ci.Tmean_AprMayJunhigh[3,], ci.group = "highppt_norm")
ci.Tmean_AprMayJunmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunmid[2,], ci.low = ci.Tmean_AprMayJunmid[1,], ci.high = ci.Tmean_AprMayJunmid[3,], ci.group = "midppt_norm")
ci.Tmean_AprMayJunlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.Tmean_AprMayJunlow[2,], ci.low = ci.Tmean_AprMayJunlow[1,], ci.high = ci.Tmean_AprMayJunlow[3,], ci.group = "lowppt_norm")
ppt_norm_Tmean_AprMayJunint <- rbind(ci.Tmean_AprMayJunhigh.df, ci.Tmean_AprMayJunmid.df, ci.Tmean_AprMayJunlow.df)
ggplot(data = ppt_norm_Tmean_AprMayJunint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = ppt_norm_Tmean_AprMayJunint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#ppt_norm and Tmean_SepOct interaction
Tmean_SepOctrng <- range(grow_train$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
tmp_norm <- mean(grow_train$tmp_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
ppt_norm <- mean(grow_train$ppt_norm)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
ppt_norm_range <- quantile(grow_train$ppt_norm, c(0.2, 0.8))
growthpredictionpnorm_highTmeanSepOct <- growthpredictionpnorm_lowTmeanSepOct <- growthpredictionpnorm_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
  growthpredictionpnorm_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[2] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
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
  
  growthpredictionpnorm_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm_range[1] + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
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
ppt_norm_prediction_trlow <- exp(growthpredictionpnorm_lowTmeanSepOct)
ppt_norm_prediction_trmid <- exp(growthpredictionpnorm_midTmeanSepOct)
ppt_norm_prediction_trhigh <- exp(growthpredictionpnorm_highTmeanSepOct)
ci.Tmean_SepOcthigh <- apply(ppt_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Tmean_SepOctmid <- apply(ppt_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_SepOctlow <- apply(ppt_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Tmean_SepOcthigh.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_SepOcthigh[2,], ci.low = ci.Tmean_SepOcthigh[1,], ci.high = ci.Tmean_SepOcthigh[3,], ci.group = "highppt_norm")
ci.Tmean_SepOctmid.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_SepOctmid[2,], ci.low = ci.Tmean_SepOctmid[1,], ci.high = ci.Tmean_SepOctmid[3,], ci.group = "midppt_norm")
ci.Tmean_SepOctlow.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.Tmean_SepOctlow[2,], ci.low = ci.Tmean_SepOctlow[1,], ci.high = ci.Tmean_SepOctlow[3,], ci.group = "lowppt_norm")
ppt_norm_Tmean_SepOctint <- rbind(ci.Tmean_SepOcthigh.df, ci.Tmean_SepOctmid.df, ci.Tmean_SepOctlow.df)
ggplot(data = ppt_norm_Tmean_AprMayJunint, aes(x = ppt_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = ppt_norm_Tmean_SepOctint, aes(x = Tmean_SepOct, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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
ggplot(data = ppt_norm_sizeint, aes(x = ppt_norm, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = ppt_norm, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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

for(i in 1:length(plotdatainterval$u_beta_Precip_NovDecJanFebMar)){
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
ci.Precip_NovDecJanFebMarhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.Precip_NovDecJanFebMarmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_NovDecJanFebMarlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.Precip_NovDecJanFebMarhigh.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarhigh[2,], ci.low = ci.Precip_NovDecJanFebMarhigh[1,], ci.high = ci.Precip_NovDecJanFebMarhigh[3,], ci.group = "hightmp_norm")
ci.Precip_NovDecJanFebMarmid.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarmid[2,], ci.low = ci.Precip_NovDecJanFebMarmid[1,], ci.high = ci.Precip_NovDecJanFebMarmid[3,], ci.group = "midtmp_norm")
ci.Precip_NovDecJanFebMarlow.df <- data.frame(Precip_NovDecJanFebMar = Precip_NovDecJanFebMar, median = ci.Precip_NovDecJanFebMarlow[2,], ci.low = ci.Precip_NovDecJanFebMarlow[1,], ci.high = ci.Precip_NovDecJanFebMarlow[3,], ci.group = "lowtmp_norm")
tmp_norm_Precip_NovDecJanFebMarint <- rbind(ci.Precip_NovDecJanFebMarhigh.df, ci.Precip_NovDecJanFebMarmid.df, ci.Precip_NovDecJanFebMarlow.df)
ggplot(data = tmp_norm_Precip_NovDecJanFebMarint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = tmp_norm_Precip_NovDecJanFebMarint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_norm and Precip_JulAug interaction
Precip_JulAugrng <- range(grow_train$Precip_JulAug,na.rm = TRUE) #setting range for tmp_normrng
Precip_JulAug <- seq(Precip_JulAugrng[1], Precip_JulAugrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictiontnorm_highPrecipJulAug <- growthpredictiontnorm_lowPrecipJulAug <- growthpredictiontnorm_midPrecipJulAug <- matrix(NA, length(plotdatainterval$u_beta_Precip_JulAug), length(Precip_JulAug)) 

for(i in 1:length(plotdatainterval$u_beta_Precip_JulAug)){
  growthpredictiontnorm_highPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
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
  
  growthpredictiontnorm_lowPrecipJulAug[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
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
tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowPrecipJulAug)
tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midPrecipJulAug)
tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highPrecipJulAug)
ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normhigh.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "hightnorm")
ci.tmp_normmid.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midtnorm")
ci.tmp_normlow.df <- data.frame(Precip_JulAug = Precip_JulAug, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowtnorm")
tmp_norm_Precip_JulAugint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
ggplot(data = tmp_norm_Precip_JulAugint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = tmp_norm_Precip_JulAugint, aes(x = Precip_JulAug, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Precip_JulAug, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_norm and Tmean_AprMayJun interaction
Tmean_AprMayJun <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictiontnorm_highTmeanAprMayJun <- growthpredictiontnorm_lowTmeanAprMayJun <- growthpredictiontnorm_midTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
  growthpredictiontnorm_highTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
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
  
  growthpredictiontnorm_lowTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
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
tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowTmeanAprMayJun)
tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midTmeanAprMayJun)
tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highTmeanAprMayJun)
ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normhigh.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "hightnorm")
ci.tmp_normmid.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midtnorm")
ci.tmp_normlow.df <- data.frame(Tmean_AprMayJun = Tmean_AprMayJun, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowtnorm")
tmp_norm_Precip_JulAugint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
ggplot(data = tmp_norm_Precip_JulAugint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = tmp_norm_Precip_JulAugint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


#tmp_norm and Tmean_SepOct interaction
Tmean_SepOct <- range(grow_train$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
Tmean_SepOct <- seq(Tmean_SepOctrng[1], Tmean_SepOctrng[2], by = 0.50)
x <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
tmp_norm_range <- quantile(grow_train$tmp_norm, c(0.2, 0.8))
growthpredictiontnorm_highTmeanSepOct <- growthpredictiontnorm_lowTmeanSepOct <- growthpredictiontnorm_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
  growthpredictiontnorm_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[2] +
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
  
  growthpredictiontnorm_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm_range[1] +
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
tmp_norm_prediction_trlow <- exp(growthpredictiontnorm_lowTmeanSepOct)
tmp_norm_prediction_trmid <- exp(growthpredictiontnorm_midTmeanSepOct)
tmp_norm_prediction_trhigh <- exp(growthpredictiontnorm_highTmeanSepOct)
ci.tmp_normhigh <- apply(tmp_norm_prediction_trhigh, 2, quantile, c(0.025,0.5,0.975)) #confidence intervals
ci.tmp_normmid <- apply(tmp_norm_prediction_trmid, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normlow <- apply(tmp_norm_prediction_trlow, 2, quantile, c(0.025, 0.5, 0.975))
ci.tmp_normhigh.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.tmp_normhigh[2,], ci.low = ci.tmp_normhigh[1,], ci.high = ci.tmp_normhigh[3,], ci.group = "hightnorm")
ci.tmp_normmid.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.tmp_normmid[2,], ci.low = ci.tmp_normmid[1,], ci.high = ci.tmp_normmid[3,], ci.group = "midtnorm")
ci.tmp_normlow.df <- data.frame(Tmean_SepOct = Tmean_SepOct, median = ci.tmp_normlow[2,], ci.low = ci.tmp_normlow[1,], ci.high = ci.tmp_normlow[3,], ci.group = "lowtnorm")
tmp_norm_Tmean_SepOctint <- rbind(ci.tmp_normhigh.df, ci.tmp_normmid.df, ci.tmp_normlow.df)
ggplot(data = tmp_norm_Tmean_SepOctint, aes(x = tmp_norm, y = median, color = ci.group)) + geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)
ggplot(data = tmp_norm_Tmean_SepOctint, aes(x = Tmean_SepOct, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_SepOct, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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
ggplot(data = tmp_norm_sizeint, aes(x = tmp_norm, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = tmp_norm, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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
ggplot(data = Precip_NovDecJanFebMar_Tmean_SepOctint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



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
ggplot(data = Precip_NovDecJanFebMar_sizeint, aes(x = Precip_NovDecJanFebMar, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Precip_NovDecJanFebMar, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



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
ggplot(data = Tmean_AprMayJun_Tmean_SepOctint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)


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
ggplot(data = Tmean_AprMayJun_sizeint, aes(x = Tmean_AprMayJun, y = median, color = ci.group)) +
  scale_color_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) +
  geom_ribbon(aes(x = Tmean_AprMayJun, ymin = ci.low, ymax = ci.high, fill = ci.group), color = NA, alpha = 0.5) +
  scale_fill_manual(values=c("#4575B4", "#FDDBC7", "#B2182B")) + geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



#Tmean_AprMayJun and size
Tmean_AprMayJunrng <- range(grow_train$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 7.5)
Tmean_SepOct <- mean(grow_train$Tmean_SepOct)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
size <- mean(grow_train$DIA_prev)
size_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictionsize_highTmeanAprMayJun <- growthpredictionsize_lowTmeanAprMayJun <- growthpredictionsize_midTmeanAprMayJun <- matrix(NA, length(plotdatainterval$u_beta_Tmean_AprMayJun), length(Tmean_AprMayJun)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_AprMayJun)){
  growthpredictionsize_highTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionsize_midTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
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
  
  growthpredictionsize_lowTmeanAprMayJun[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*x_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
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


#size and Tmean_SepOct
Tmean_SepOctrng <- range(grow_train$Tmean_SepOct,na.rm = TRUE) #setting range for tmp_normrng
Tmean_SepOct <- seq(sizerng[1], sizerng[2], by = 7.5)
size <- mean(grow_train$DIA_prev)
tmp_norm <- mean(grow_train$tmp_norm)
ppt_norm <- mean(grow_train$ppt_norm)
Precip_JulAug <- mean(grow_train$Precip_JulAug)
Precip_NovDecJanFebMar <- mean(grow_train$Precip_NovDecJanFebMar)
Tmean_AprMayJun <- mean(grow_train$Tmean_AprMayJun)
size_range <- quantile(grow_train$DIA_prev, c(0.2, 0.8))
growthpredictionsize_highTmeanSepOct <- growthpredictionsize_lowTmeanSepOct <- growthpredictionsize_midTmeanSepOct <- matrix(NA, length(plotdatainterval$u_beta_Tmean_SepOct), length(Tmean_SepOct)) 

for(i in 1:length(plotdatainterval$u_beta_Tmean_SepOct)){
  growthpredictionsize_highTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[2] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[2] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[2] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size_range[2]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size_range[2]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size_range[2]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size_range[2]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
  
  growthpredictionsize_midTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
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
  
  growthpredictionsize_lowTmeanSepOct[i,] <- plotdatainterval[i,"u_beta_ppt_norm"]*ppt_norm + plotdatainterval[i,"u_beta_tmp_norm"]*tmp_norm +
    plotdatainterval[i,"u_beta_Precip_JulAug"]*Precip_JulAug + plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar"]*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Tmean_AprMayJun"]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_Tmean_SepOct"]*Tmean_SepOct +
    plotdatainterval[i,"u_beta_DIA_prev"]*size_range[1] + plotdatainterval[i,"u_beta_ppt_norm_tmp_norm"]*ppt_norm*tmp_norm +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_JulAug"]*ppt_norm*Precip_JulAug + plotdatainterval[i,"u_beta_ppt_norm_DIA_prev"]*ppt_norm*size_range[1] +
    plotdatainterval[i,"u_beta_ppt_norm_Precip_NovDecJanFebMar"]*ppt_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_ppt_norm_Tmean_AprMayJun"]*ppt_norm*Tmean_AprMayJun +  plotdatainterval[i,"u_beta_ppt_norm_Tmean_SepOct"]*ppt_norm*Tmean_SepOct +
    plotdatainterval[i,"u_beta_tmp_norm_DIA_prev"]*tmp_norm*size_range[1] + plotdatainterval[i,"u_beta_tmp_norm_Precip_JulAug"]*tmp_norm*Precip_JulAug +
    plotdatainterval[i,"u_beta_tmp_norm_Precip_NovDecJanFebMar"]*tmp_norm*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_tmp_norm_Tmean_AprMayJun"]*tmp_norm*Tmean_AprMayJun + plotdatainterval[i,"u_beta_tmp_norm_Tmean_SepOct"]*tmp_norm*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_DIA_prev_Precip_JulAug"]*size_range[1]*Precip_JulAug + plotdatainterval[i,"u_beta_DIA_prev_Precip_NovDecJanFebMar"]*size_range[1]*Precip_NovDecJanFebMar + 
    plotdatainterval[i,"u_beta_DIA_prev_Tmean_AprMayJun"]*size_range[1]*Tmean_AprMayJun + plotdatainterval[i,"u_beta_DIA_prev_Tmean_SepOct"]*size_range[1]*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Precip_NovDecJanFebMar"]*Precip_JulAug*Precip_NovDecJanFebMar +
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_AprMayJun"]*Precip_JulAug*Tmean_AprMayJun + 
    plotdatainterval[i,"u_beta_Precip_JulAug_Tmean_SepOct"]*Precip_JulAug*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_AprMayJun"]*Precip_NovDecJanFebMar*Tmean_AprMayJun +  
    plotdatainterval[i,"u_beta_Precip_NovDecJanFebMar_Tmean_SepOct"]*Precip_NovDecJanFebMar*Tmean_SepOct + 
    plotdatainterval[i,"u_beta_Tmean_AprMayJun_Tmean_SepOct"]*Tmean_AprMayJun*Tmean_SepOct
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
