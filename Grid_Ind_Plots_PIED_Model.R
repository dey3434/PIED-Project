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



#INDIVIDUAL RESPONSE PLOTS
grow.monsoon$LONbin <- ifelse(grow.monsoon$LON > -109, "-109 to -104", "-114 to -109")
grow.monsoon$LATbin <- ifelse(grow.monsoon$LAT > 37, "37 to 41", "32 to 37")
grow.monsoon$LATLONbin <- paste(grow.monsoon$LONbin, grow.monsoon$LATbin)
ind.samples <- unique(grow.monsoon[,c("LATLONbin", "treeCD")]) %>% group_by(LATLONbin)

get.ind.tmp.response<- function(j){
  tree.subset <- ind.samples[j,]
  tree.grow <- grow.monsoon %>% filter(LATLONbin == tree.subset$LATLONbin & treeCD == tree.subset$treeCD)
  
  Tmean_AprMayJunrng <- range(tree.grow$Tmean_AprMayJun,na.rm = TRUE) #setting range for tmp_normrng
  Tmean_AprMayJun <- seq(Tmean_AprMayJunrng[1], Tmean_AprMayJunrng[2], by = 0.1)
  x <- mean(tree.grow$DIA_prev)
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
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = tmp_norm, group = ci.group)) + geom_line(alpha = 0.5) + 
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
  x <- mean(tree.grow$DIA_prev)
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
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = Precip_JulAug, group = ci.group)) + geom_line() + 
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
ggplot(data = merged.response.samples, aes(x = Tmean_AprMayJun, y = median, color = Precip_NovDecJanFebMar, group = ci.group)) + geom_line() + 
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
ggplot(data = merged.response.samples, aes(x = Precip_JulAug, y = median, color = Tmean_SepOct, group = ci.group)) + #geom_ribbon(aes(x = tmp_yr, ymin = ci.low, ymax = ci.high, fill = tmp_norm, group = ci.group),color = NA, alpha = 0.5) + 
  geom_line() + mytheme + ylab("Predicted Growth") + ylim(0, 2)



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


