#BIG GRID PLOT (note: code is altered to make one plot at a time)
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
growthpredictionfilter2$ppt_norm_character <- factor(growthpredictionfilter2$ppt_norm_character, levels = c("low ppt_norm", "mid ppt_norm", "high ppt_norm"))
growthpredictionfilter2$tmp_norm_character <- ifelse(growthpredictionfilter2$tmp_norm == quantile(grow_train$tmp_norm, c(0.2)), "low tmp_norm", 
                                                     ifelse(growthpredictionfilter2$tmp_norm == quantile(grow_train$tmp_norm, c(0.6)), "mid tmp_norm",
                                                            "high tmp_norm"))
growthpredictionfilter2$tmp_norm_character <- factor(growthpredictionfilter2$tmp_norm_character, levels = c("low tmp_norm", "mid tmp_norm", "high tmp_norm"))
growthpredictionfilter2$Precip_NovDecJanFebMar_character <- ifelse(growthpredictionfilter2$Precip_NovDecJanFebMar == quantile(grow_train$Precip_NovDecJanFebMar, c(0.2)), "low Precip_NovDecJanFebMar", 
                                                                   ifelse(growthpredictionfilter2$Precip_NovDecJanFebMar == quantile(grow_train$Precip_NovDecJanFebMar, c(0.6)), "mid Precip_NovDecJanFebMar",
                                                                          "high Precip_NovDecJanFebMar"))
growthpredictionfilter2$Precip_NovDecJanFebMar_character <- factor(growthpredictionfilter2$Precip_NovDecJanFebMar_character, levels = c("low Precip_NovDecJanFebMar", "mid Precip_NovDecJanFebMar", "high Precip_NovDecJanFebMar"))
growthpredictionfilter2$Precip_JulAug_character <- ifelse(growthpredictionfilter2$Precip_JulAug == quantile(grow_train$Precip_JulAug, c(0.2)), "low Precip_JulAug", 
                                                          ifelse(growthpredictionfilter2$Precip_JulAug == quantile(grow_train$Precip_JulAug, c(0.6)), "mid Precip_JulAug",
                                                                 "high Precip_JulAug"))
growthpredictionfilter2$Precip_JulAug_character <- factor(growthpredictionfilter2$Precip_JulAug_character, levels = c("low Precip_JulAug", "mid Precip_JulAug", "high Precip_JulAug"))
growthpredictionfilter2$Tmean_AprMayJun_character <- ifelse(growthpredictionfilter2$Tmean_AprMayJun == quantile(grow_train$Tmean_AprMayJun, c(0.2)), "low Tmean_AprMayJun", 
                                                            ifelse(growthpredictionfilter2$Tmean_AprMayJun == quantile(grow_train$Tmean_AprMayJun, c(0.6)), "mid Tmean_AprMayJun",
                                                                   "high Tmean_AprMayJun"))
growthpredictionfilter2$Tmean_AprMayJun_character <- factor(growthpredictionfilter2$Tmean_AprMayJun_character, levels = c("low Tmean_AprMayJun", "mid Tmean_AprMayJun", "high Tmean_AprMayJun"))
growthpredictionfilter2$Tmean_SepOct_character <- ifelse(growthpredictionfilter2$Tmean_SepOct == quantile(grow_train$Tmean_SepOct, c(0.2)), "low Tmean_SepOct", 
                                                         ifelse(growthpredictionfilter2$Tmean_SepOct == quantile(grow_train$Tmean_SepOct, c(0.6)), "mid Tmean_SepOct",
                                                                "high Tmean_SepOct"))
growthpredictionfilter2$Tmean_SepOct_character <- factor(growthpredictionfilter2$Tmean_SepOct_character, levels = c("low Tmean_SepOct", "mid Tmean_SepOct", "high Tmean_SepOct"))
growthpredictionfilter2$DIA_prev_character <- ifelse(growthpredictionfilter2$DIA_prev == quantile(grow_train$DIA_prev, c(0.2)), "low DIA_prev", 
                                                     ifelse(growthpredictionfilter2$DIA_prev == quantile(grow_train$DIA_prev, c(0.6)), "mid DIA_prev",
                                                            "high DIA_prev"))
growthpredictionfilter2$DIA_prev_character <- factor(growthpredictionfilter2$DIA_prev_character, levels = c("low DIA_prev", "mid DIA_prev", "high DIA_prev"))
growthpredictionfilter3 <- growthpredictioncbind %>% filter(Tmean_AprMayJun == quantile(grow_train$Tmean_AprMayJun, c(0.6)) & Tmean_SepOct == quantile(grow_train$Tmean_SepOct, c(0.6)), 
                                                            DIA_prev == quantile(grow_train$DIA_prev, c(0.6)))
ggplot(data = growthpredictionfilter3, aes(x = Precip_NovDecJanFebMar, y = `50%`, color = ppt_norm_character)) + geom_ribbon(data = growthpredictionfilter3, aes(x = Precip_NovDecJanFebMar, ymin = `2.5%`, ymax = `97.5%`, fill = ppt_norm_character),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2) + facet_grid(tmp_norm_character ~ Precip_JulAug_character)

