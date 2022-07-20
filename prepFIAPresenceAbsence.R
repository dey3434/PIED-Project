# prep code to recalculate presence abserce of PIED FIA data from Emily Schultz's code
library(raster)
library(tidyverse)
library(dplyr)

fiadb <- readRDS("data/InWeUS_FIAdb.rds") # this was pulled using the fiadb package
TREE <- fiadb$TREE
PLOT <- fiadb$PLOT
unique(PLOT$STATECD)

fourcorners <- c(4, 8, 35, 49)

PLOT_COMBINED <- PLOT %>% filter(STATECD %in% fourcorners)
TREE_COMBINED <- TREE %>% filter(STATECD %in% fourcorners)

write.csv(TREE_COMBINED, "data/TREE_COMBINED.csv")
write.csv(PLOT_COMBINED, "data/PLOT_COMBINED.csv")

# Read in data
all <- read.csv("data/TREE_COMBINED.csv", header = T, stringsAsFactors = F)

# Keep only plot and species code
all <- all[, c("PLT_CN", "SPCD")]

# Read in plots
plots <- read.csv("data/PLOT_COMBINED.csv")
plots <- plots[, c("CN", "LAT", "LON")]

# Processing
# KH notes....this is taking awhile...trying to speed up
all$LAT <- NA
all$LON <- NA
# progress <- 1
# no <- length(plots$CN)
# for (i in plots$CN) {
#   if (nrow(all[all$PLT_CN == i,]) > 0) {
#     all[all$PLT_CN == i,]$LAT <- subset(plots, CN == i)$LAT
#     all[all$PLT_CN == i,]$LON <- subset(plots, CN == i)$LON
#   }
#   progress <- progress + 1
#   print(progress / no * 100)
# }

# I think we just want lat and lon: so we can match it up
all$LAT <- plots$LAT[match(all$PLT_CN, plots$CN)]
all$LON <- plots$LON[match(all$PLT_CN, plots$CN)]



# Export
#write.csv(all, "Processed/validationRecords.csv", row.names = F)

# Read in a raster for alignment
pres <- raster("./Output/tifs/PIED.clim_lambda_gam.tif")# I dont have this (KH), but is shoudl be the same size as the file emily sent
pres <- raster("presenceAbsenceRaster2.tif")
pres <- setValues(pres, NA)

# Make records spatial
proj <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0"
allSp <- SpatialPointsDataFrame(all[, c("LON", "LAT")], as.data.frame(all[, "SPCD"]), proj4string = CRS("+proj=longlat +datum=NAD83"))

# For each tree, find in which raster cell it falls
allSp$cellNumber <- NA
# progress <- 1
# no <- nrow(allSp)
# for (i in 1:nrow(allSp)) {
#  tmp <- as.integer(extract(pres, allSp[i,], cellnumbers = T)[[1]])
#  allSp[i, "cellNumber"] <- tmp
#  progress <- progress + 1
#  print(progress / no * 100)
# } 


# i dont think we need to individually do the extraction...kh...extract should allow us to do this all at once
tmp <- raster::extract(pres, allSp, cellnumbers = T)
allSp$cellNumber <- tmp[,"cells"]
# trying to convert to a raster for each tree

# For each raster cell, find whether it contains PIED (SPCD = 106)
names(allSp) <- c("SPCD", "cellNumber")
PIED <- 106
system.time(
for (i in 1:ncell(pres)) {
  tmp <- subset(allSp, cellNumber == i)
  if (nrow(tmp) > 0) {
    if (PIED %in% unique(tmp$SPCD)) pres[i] <- 1
    else pres[i] <- 0
  } else pres[i] <- NA # If there are no trees/FIA plots in allSp in thot fall into this grid cell 
}
)

# takes about 40 seconds...so I wont bother turning this into an apply function:
# user  system elapsed 
# 68.563  40.426 113.896 

# Export presence/absence raster
writeRaster(pres, "data/presenceAbsenceRasterNA.tif", overwrite = T)
raster::plot(pres) # so now NAs mean there is no FIA data, zeros mean there is no PIED presernce and 1's mean there is PIED presence


# not sure we needed this code for Sharmilas project so commenting out (KH)
# climate_wna<-plots %>%
#   group_by(PLOT) %>%
#   summarise(lat=mean(LAT,na.rm=T),long=mean(LON,na.rm=T),el=(mean(ELEV,na.rm=T)*0.3048)) %>%
#   mutate(ID2=".")
# 
# names(climate_wna)[1]<-"ID1"
# 
# write.csv(climate_wna[1:4215,],"./ClimateData/FIA_climateWNA1.csv",row.names=F)
# write.csv(climate_wna[4216:8430,],"./ClimateData/FIA_climateWNA2.csv",row.names=F)
# write.csv(climate_wna[8431:12645,],"./ClimateData/FIA_climateWNA3.csv",row.names=F)
# write.csv(climate_wna[12646:16863,],"./ClimateData/FIA_climateWNA4.csv",row.names=F)
# 

# now extract climate data for those points:


#### Use PRISM climate data to create rasters of monthly climate variables in PIED study region
# from emily's code:https://github.com/emilylschultz/DemographicRangeModel/blob/master/Code/ClimateProcessing/historic.R
library(raster)

### PRISM download January 22, 2019
### January 1981 through June 2018
### (37*12) + 6 = 450 files

# Search for PRISM files
PRISM.path <-  "./PRISM_data/"
ppt.path <-  "/Users/kah/Documents/docker_pecan/pecan/PRISM_data/PRISM_ppt_stable_4kmM2_189501_198012_bil/"
ppt.path.new <-  "/Users/kah/Documents/docker_pecan/pecan/PRISM_data/PRISM_ppt_stable_4kmM3_198101_201910_bil/"


pptFiles.old <- list.files(path = ppt.path, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
pptFiles.new <- list.files(path = ppt.path.new, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
pptFiles<- c(pptFiles.old, pptFiles.new)

temp.path <-  "/Users/kah/Documents/docker_pecan/pecan/PRISM_data/PRISM_tmean_stable_4kmM3_189501_198012_bil/"
temp.path.new <-  "/Users/kah/Documents/docker_pecan/pecan/PRISM_data/PRISM_tmean_stable_4kmM3_198101_202002_bil/"

tempFiles.old <- list.files(path = temp.path, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
tempFiles.new <- list.files(path = temp.path.new, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
tempFiles<- c(tempFiles.old, tempFiles.new)

#tmpFiles <- list.files(path = PRISM.path, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)

#vpdminFiles <- list.files(path = PRISM.path, pattern = glob2rx("*vpdmin*.bil"), full.names = TRUE)
#vpdmaxFiles <- list.files(path = PRISM.path, pattern = glob2rx("*vpdmax*.bil"), full.names = TRUE)

# Stack monthly data
pptStack <- stack()
#for (i in pptFiles[1:2]) {
#  print(i)
  #pptStack <- stack(pptStack, raster(i))
#}

pptStacklist <- lapply(pptFiles, function(i){stack(raster(i))})
pptStack <- stack(pptStacklist)

tmpStack <- stack()
#for (i in tempFiles) {
#  print(i)
#  tmpStack <- stack(tmpStack, raster(i))
#}

tempStacklist <- lapply(tempFiles, function(i){stack(raster(i))})
tempStack <- stack(tempStacklist)

#vpdStack <- stack()
#for (i in 1:length(vpdmaxFiles)) {
#  print(i)
#  rast<-raster(vpdmaxFiles[i])
#  if(i == 431 | i == 432){crs(rast)<-crs(raster(vpdmaxFiles[1]))}
#  vpdStack <- stack(vpdStack, rast)
#}


# extract by plot lat and lon:
PApied <- raster("data/presenceAbsenceRasterNA.tif")
plot(PApied)

crs(pptStack)
crs(PApied)

PApied <- projectRaster(PApied, crs = crs(pptStack))

FIA_PApied <- as.data.frame(PApied,xy=TRUE)
names(FIA_PApied )<-c("lon","lat","PApied")

ppt.extracted <- data.frame(extract(pptStack, FIA_PApied[,c("lon", "lat")]))
ppt.extracted$lon <- FIA_PApied$lon
ppt.extracted$lat <- FIA_PApied$lat
ppt.extracted$name<- FIA_PApied$name


saveRDS(ppt.extracted, "data/pied_extracted.ppt.PAdata.rds")

ppt.extracted <- readRDS( "data/pied_extracted.ppt.PAdata.rds")


temp.extracted <- data.frame(raster::extract(tempStack, FIA_PApied[,c("lon", "lat")]))
temp.extracted$lon <- FIA_PApied$lon
temp.extracted$lat <- FIA_PApied$lat
temp.extracted$name<- FIA_PApied$name


saveRDS(temp.extracted, "data/pied_extracted.temp.PAdata.rds")
temp.extracted <- readRDS( "data/pied_extracted.temp.data.rds")


FIA_PApied <- as.data.frame(PApied,xy=TRUE)
names(FIA_PApied )<-c("lon","lat","PApied")

summary(FIA_PApied)

# save this as a data frame so that you could read it in without needing the raster package:
write.csv(FIA_PApied, "data/PresenceAbsernce_PIEDFIA.csv")

## Logistic models (lambda vs occurrence)
k=5
pa_monsoonp <- gam(PApied~s(Precip_JulAug,k=k),family=binomial(link = cloglog),
            data=FIA_PApied)
pa_winterp <- gam(PApied~s(Precip_NovDecJanFebMar,k=k),family=binomial(link = cloglog),
             data=FIA_PApied)
pa_springt <- gam(PApied~s(Tmean_AprMayJun,k=k),family=binomial(link = cloglog),
             data=FIA_PApied)
pa_fallt <- gam(PApied~s(Tmean_SepOct, k=k),family=binomial(link = cloglog),
            data=FIA_PApied)

## Set up data frame with binned values
ncuts=50 # number of cut to create bins

# divide into intervals based on number of cuts
chopsize_monsoonp<-cut(FIA_PApied$Precip_JulAug,ncuts)
chopsize_winterp<-cut(FIA_PApied$Precip_NovDecJanFebMar,ncuts)
chopsize_springt<-cut(FIA_PApied$Tmean_AprMayJun,ncuts)
chopsize_fallt<-cut(FIA_PApied$Tmean_SepOct,ncuts)

count_binned_lam_monsoonp<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_monsoonp),length)) # calculate number of data points in each bin
lam_monsoonp_binned<-as.vector(sapply(split(FIA_PApied$Precip_JulAug,chopsize_monsoonp),mean,na.rm=T)) # calculate mean lambda in each bin
pres_monsoonp_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_monsoonp),mean,na.rm=T)) # calculate mean prob of occurrence in each bin

count_binned_lam_winterp<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_lam_winterp),length))
lam_winterp_binned<-as.vector(sapply(split(FIA_PApied$Precip_NovDecJanFebMar,chopsize_lam_winterp),mean,na.rm=T))
pres_winterp_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_lam_winterp),mean,na.rm=T))

count_binned_lam_springt<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_lam_springt),length))
lam_springt_binned<-as.vector(sapply(split(FIA_PApied$Tmean_AprMayJun,chopsize_lam_springt),mean,na.rm=T))
pres_springt_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_lam_springt),mean,na.rm=T))

count_binned_lam_fallt<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_lam_fallt),length))
lam_fallt_binned<-as.vector(sapply(split(FIA_PApied$Tmean_SepOct,chopsize_lam_fallt),mean,na.rm=T))
pres_fallt_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_lam_fallt),mean,na.rm=T))

# combine into dataframe
pres_binned<-data.frame(count_lam=c(count_binned_lam_monsoonp,count_binned_lam_winterp,count_binned_lam_springt, count_binned_lam_fallt),
                        lam=c(lam_monsoonp_binned,lam_winterp_binned,lam_springt_binned, lam_fallt_binned),
                        pres=c(pres_monsoonp_binned,pres_winterp_binned,pres_springt_binned, pres_fallt_binned),
                        pred=c(invlogit(predict(pa_monsoonp,newdata=data.frame(lambda_monsoonp=lam_monsoonp_binned))),
                               invlogit(predict(pa_winterp,newdata=data.frame(lambda_winterp=lam_winterp_binned))),
                               invlogit(predict(pa_springt,newdata=data.frame(lambda_springt=lam_springt_binned))),
                               invlogit(predict(pa_fallt,newdata=data.frame(lambda_fallt=lam_fallt_binned)))),
                        model=c(rep("monsoonp",ncuts),rep("winterp",ncuts),rep("springt",ncuts),rep("fallt", ncuts))

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
