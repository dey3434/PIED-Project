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
#pres <- raster("./Output/tifs/PIED.clim_lambda_gam.tif")# I dont have this (KH), but is shoudl be the same size as the file emily sent
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
#ppt.extracted$name<- FIA_PApied$PApied


saveRDS(ppt.extracted, "data/pied_extracted.ppt.PAdata.rds")

ppt.extracted <- readRDS( "data/pied_extracted.ppt.PAdata.rds")


temp.extracted <- data.frame(raster::extract(tempStack, FIA_PApied[,c("lon", "lat")]))
temp.extracted$lon <- FIA_PApied$lon
temp.extracted$lat <- FIA_PApied$lat
#temp.extracted$name<- FIA_PApied$PApied


saveRDS(temp.extracted, "data/pied_extracted.temp.PAdata.rds")
temp.extracted <- readRDS( "data/pied_extracted.temp.PAdata.rds")

# now summarise the climate data to the variables we are interests in:


colnames(temp.extracted) 
years <- rep(1895:2019, each = 12)
months <- rep(1:12, 125)

colnames.clim <- c(paste0(years, "_", months), "2020_1", "2020_2")
length(temp.extracted)
length(colnames.clim)
colnames(temp.extracted)[1:1502] <- colnames.clim
#colnames(ppt.extracted)[1:1502] <- colnames.clim
temp.m <- reshape2::melt(temp.extracted, id.vars = c("lon", "lat"))
#temp.long <- temp.m %>% separate(variable, c("year", "month"))
temp.long <- temp.m %>% filter(! variable %in% c("2020_1", "2020_2"))
temp.long$variable <- as.character(temp.long$variable)
var.split <- str_split_fixed(temp.long$variable, "_", 2)#temp.long %>% separate(variable, c("year", "month"))
temp.long$year <- var.split[,1]
temp.long$month <- var.split[,2]
colnames(temp.long) <- c("lon", "lat","variable", "Tave", "year", "month")

# do the same for precipitation
colnames(ppt.extracted) 
years <- rep(1895:2018, each = 12)
months <- rep(1:12, 124)

colnames.clim <- c(paste0(years, "_", months), "2019_1", "2019_2", "2019_3", "2019_4", "2019_5", "2019_6", "2019_7", "2019_8", "2019_9", "2019_10")
length(ppt.extracted)
colnames(ppt.extracted)[1:1498] <- colnames.clim
#colnames(ppt.extracted)[1:1502] <- colnames.clim
ppt.m <- reshape2::melt(ppt.extracted, id.vars = c("lon", "lat"))
ppt.m$variable <- as.character(ppt.m$variable)
#ppt.long <- ppt.m %>% separate(variable, c("year", "month")) This kept crashing R
var.split <- str_split_fixed(ppt.m$variable, "_", 2)#temp.long %>% separate(variable, c("year", "month"))
ppt.m$year <- var.split[,1]
ppt.m$month <- var.split[,2]
ppt.long <- ppt.m
colnames(ppt.long) <- c("lon", "lat", "variable", "ppt", "year", "month")

ppt.long <- ppt.long %>% select(-variable)
temp.long <- temp.long %>% select(-variable)
temp.ppt <- left_join(temp.long, ppt.long, by =c("lon", "lat", "year", "month"))


# assign water year
wtr_yr <- function(df, start_month=9) {
  # Year offset
  offset = ifelse(as.numeric(df$month) >= start_month - 1, 1, 0)
  # Water year
  adj.year = as.numeric(df$year) + offset
  # Return the water year
  adj.year
}

temp.ppt$water_year <- wtr_yr(temp.ppt)

# May april jun - arid foresummer
foresummer <- temp.ppt %>% filter(month %in% 4:6)%>% group_by(lon, lat,  year)%>% summarise(Tmean_AprMayJun = mean(Tave),
                                                                                                  Precip_AprMayJun = mean(ppt)) 


# Jul August - monsson
monsoon <- temp.ppt %>% filter(month %in% 7:8)%>% group_by(lon, lat, year)%>% summarise(Tmean_JulAug = mean(Tave),
                                                                                              Precip_JulAug = mean(ppt)) 

# Sept Oct- Fall
fall <-temp.ppt %>% filter(month %in% 9:10)%>% group_by(lon, lat,  year)%>% summarise(Tmean_SepOct = mean(Tave),
                                                                                           Precip_SepOct = mean(ppt)) 

# december january february - winter:

winter <- temp.ppt %>% filter(month %in% c(11, 12, 1, 2, 3)) %>% group_by(lon, lat, water_year) %>% summarise(Tmean_NovDecJanFebMar = mean(Tave),
                                                                                                           Precip_NovDecJanFebMar = mean(ppt)) 
colnames(winter)[3] <- "year"

winter$year <- as.character(winter$year)

# other requested variables:
# november - march
nov_mar <- temp.ppt %>% filter(month %in% c(1,2,2,11,12)) %>% group_by(lon, lat,  water_year) %>% summarise(Tmean_NovDecJanFebMar = mean(Tave),
                                                                                                                 Precip_NovDecJanFebMar = mean(ppt)) 
colnames(nov_mar)[3] <- "year"

nov_mar$year <- as.character(nov_mar$year)

# april may june
apr_may_jun <- temp.ppt %>% filter(month %in% c(4,5,6)) %>% group_by(lon, lat,  water_year) %>% summarise(Tmean_AprMayJun  = mean(Tave),
                                                                                                               Precip_AprMayJun = mean(ppt)) 
colnames(apr_may_jun )[3] <- "year"

apr_may_jun $year <- as.character(apr_may_jun $year)

# September to October of current year

sep_oct <- temp.ppt %>% filter(month %in% c(9,10)) %>% group_by(lon, lat, year) %>% summarise(Tmean_SepOct  = mean(Tave),
                                                                                                    Precip_SepOct = mean(ppt)) 
colnames(sep_oct)[3] <- "year"

sep_oct$year <- as.character(sep_oct$year)


# water year precip
wateryr <- temp.ppt %>% filter(month %in% c(1:12)) %>% group_by(lon, lat, year) %>% summarise(Tmean_wateryr  = mean(Tave),
                                                                                              Precip_wateryr = mean(ppt)) 
colnames(wateryr)[3] <- "year"

wateryr$year <- as.character(wateryr$year)


# monthly of current year

temp.wide <- temp.ppt %>% #filter(! year %in% "T1") %>% 
  dplyr::select(lat, lon,  year, Tave, month)  %>%
  group_by(lon, lat)  %>% spread (month, Tave)                         

ppt.wide <- temp.ppt  %>% #filter(! year %in% "T1") %>% 
  dplyr::select(lat, lon, year, ppt, month)  %>%
  group_by(lon, lat)  %>% spread (month, ppt) 

colnames(ppt.wide)[4:15] <- paste0("PPT_", colnames(ppt.wide)[4:15])
colnames(temp.wide)[4:15] <- paste0("TMEAN_", colnames(ppt.wide)[4:15])

# get the previous year's temp and precipt
temp.wide.prev <- temp.wide
temp.wide.prev$nextyear <- as.numeric(temp.wide.prev$year)+1                        
temp.wide.prev$actual.year <- temp.wide.prev$year
temp.wide.prev$year <- temp.wide.prev$nextyear

colnames(temp.wide.prev)[4:15] <- paste0("PREV_",colnames(temp.wide.prev)[4:15])



temp.previous <- temp.wide.prev %>% dplyr::select(-nextyear, -actual.year)
temp.previous$year <- as.character(temp.previous$year)

ppt.wide.prev <- ppt.wide
ppt.wide.prev$nextyear <- as.numeric(ppt.wide.prev$year)+1                        
ppt.wide.prev$actual.year <- ppt.wide.prev$year
ppt.wide.prev$year <- ppt.wide.prev$nextyear

colnames(ppt.wide.prev)[4:15] <- paste0("PREV_",colnames(ppt.wide.prev)[4:15])
ppt.previous <- ppt.wide.prev %>% dplyr::select(-nextyear, -actual.year)
ppt.previous$year <- as.character(ppt.previous$year)


# combine all together:
mergeCols = c("lon", "lat", "year")
fall.wint <- left_join(winter, fall, by = mergeCols)
mons.fall.wint <- left_join(monsoon, fall.wint, by = mergeCols)
all.tmean <- left_join(foresummer, mons.fall.wint, by = mergeCols)
all.tmean.2 <- left_join(all.tmean, nov_mar, by = mergeCols)
all.tmean.3 <- left_join(all.tmean.2, sep_oct, by = mergeCols)
all.tmean.4 <- left_join(all.tmean.3, apr_may_jun, by = mergeCols)
all.tmean.5 <- left_join(all.tmean.4, ppt.wide, by = mergeCols)
all.tmean.6 <- left_join(all.tmean.5, ppt.previous, by = mergeCols)
all.tmean.7 <- left_join(all.tmean.6, temp.wide, by = mergeCols)
full.ppt.tmean <- left_join(all.tmean.7, temp.previous, by = mergeCols)
full.ppt.tmean2 <- left_join(full.ppt.tmean, wateryr, by = mergeCols)



write.csv(full.ppt.tmean2, "data/FIAPIED_ll_climate.csv", row.names = FALSE)


# calculate NORMALS FOR ALL VARIABLES
NORMALS <- full.ppt.tmean2 %>% group_by(lon, lat) %>% dplyr::summarise(MAP = mean(Precip_wateryr, na.rm = TRUE), 
                                                    MAT = mean(Tmean_wateryr, na.rm = TRUE), 
                                                    MATfall = mean(Tmean_SepOct, na.rm = TRUE), 
                                                    MATspring = mean(Tmean_AprMayJun, na.rm = TRUE), 
                                                    MAPmonsoon = mean(Precip_JulAug, na.rm = TRUE), 
                                                    MAPwinter = mean(Precip_NovDecJanFebMar, na.rm = TRUE))

# very simple maps of cliamte
ggplot(NORMALS, aes(x = lon, y = lat, fill = MAP))+geom_raster()
ggplot(NORMALS, aes(x = lon, y = lat, fill = MATfall))+geom_raster()
ggplot(NORMALS, aes(x = lon, y = lat, fill = MATspring))+geom_raster()




# now merge with the PApied data:
FIA_PApied <- read.csv( "data/PresenceAbsernce_PIEDFIA.csv")
climate.PApied <- left_join(FIA_PApied, NORMALS, by = c("lat", "lon"))
write.csv(climate.PApied, "data/Normals_PA_PIED.csv", row.names = FALSE)


FIA_PApied <- read.csv("work/PIED-Project/Normals_PA_PIED.csv")
## Logistic models (lambda vs occurrence)
k=5
pa_monsoonp <- gam(PApied~s(MAPmonsoon,k=k),family=binomial(link = cloglog),
            data=FIA_PApied)
pa_winterp <- gam(PApied~s(MAPwinter,k=k),family=binomial(link = cloglog),
             data=FIA_PApied)
pa_springt <- gam(PApied~s(MATspring,k=k),family=binomial(link = cloglog),
             data=FIA_PApied)
pa_fallt <- gam(PApied~s(MATfall, k=k),family=binomial(link = cloglog),
            data=FIA_PApied)
pa_map <- gam(PApied~s(MAP,k=k),family=binomial(link = cloglog),
            data=FIA_PApied)
pa_mat <- gam(PApied~s(MAT,k=k),family=binomial(link = cloglog),
             data=FIA_PApied)

## Set up data frame with binned values
ncuts=50 # number of cut to create bins

# divide into intervals based on number of cuts
chopsize_monsoonp<-cut(FIA_PApied$MAPmonsoon,ncuts)
chopsize_winterp<-cut(FIA_PApied$MAPwinter,ncuts)
chopsize_springt<-cut(FIA_PApied$MATspring,ncuts)
chopsize_fallt<-cut(FIA_PApied$MATfall,ncuts)
chopsize_map<-cut(FIA_PApied$MAP,ncuts)
chopsize_mat<-cut(FIA_PApied$MAT,ncuts)


count_binned_lam_monsoonp<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_monsoonp),length)) # calculate number of data points in each bin
lam_monsoonp_binned<-as.vector(sapply(split(FIA_PApied$MAPmonsoon,chopsize_monsoonp),mean,na.rm=T)) # calculate mean lambda in each bin
pres_monsoonp_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_monsoonp),mean,na.rm=T)) # calculate mean prob of occurrence in each bin

count_binned_lam_winterp<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_winterp),length))
lam_winterp_binned<-as.vector(sapply(split(FIA_PApied$MAPwinter,chopsize_winterp),mean,na.rm=T))
pres_winterp_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_winterp),mean,na.rm=T))

count_binned_lam_springt<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_springt),length))
lam_springt_binned<-as.vector(sapply(split(FIA_PApied$MATspring,chopsize_springt),mean,na.rm=T))
pres_springt_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_springt),mean,na.rm=T))

count_binned_lam_fallt<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_fallt),length))
lam_fallt_binned<-as.vector(sapply(split(FIA_PApied$MATfall,chopsize_fallt),mean,na.rm=T))
pres_fallt_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_fallt),mean,na.rm=T))

count_binned_lam_map<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_map),length))
lam_map_binned<-as.vector(sapply(split(FIA_PApied$MAP,chopsize_map),mean,na.rm=T))
pres_map_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_map),mean,na.rm=T))

count_binned_lam_mat<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_mat),length))
lam_mat_binned<-as.vector(sapply(split(FIA_PApied$MAT,chopsize_mat),mean,na.rm=T))
pres_mat_binned<-as.vector(sapply(split(FIA_PApied$PApied,chopsize_mat),mean,na.rm=T))

# combine into dataframe
pres_binned<-data.frame(count_lam=c(count_binned_lam_monsoonp,count_binned_lam_winterp,count_binned_lam_springt, count_binned_lam_fallt,
                                    count_binned_lam_map, count_binned_lam_mat),
                        lam=c(lam_monsoonp_binned,lam_winterp_binned,lam_springt_binned, lam_fallt_binned,
                              lam_map_binned, lam_mat_binned),
                        pres=c(pres_monsoonp_binned,pres_winterp_binned,pres_springt_binned, pres_fallt_binned,
                               pres_map_binned, pres_mat_binned),
                        pred=c(invlogit(predict(pa_monsoonp,newdata=data.frame(MAPmonsoon=lam_monsoonp_binned))),
                               invlogit(predict(pa_winterp,newdata=data.frame(MAPwinter=lam_winterp_binned))),
                               invlogit(predict(pa_springt,newdata=data.frame(MATspring=lam_springt_binned))),
                               invlogit(predict(pa_fallt,newdata=data.frame(MATfall=lam_fallt_binned))),
                               invlogit(predict(pa_map,newdata=data.frame(MAP=lam_map_binned))),
                               invlogit(predict(pa_mat,newdata=data.frame(MAT=lam_mat_binned)))),
                        model=c(rep("monsoonp",ncuts),rep("winterp",ncuts),rep("springt",ncuts),rep("fallt", ncuts),
                                rep("map", ncuts), rep("mat", ncuts)))

# Plot
pres_plot_monsoonp<-ggplot(data=subset(pres_binned,model=="monsoonp"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  annotate("label", x = 75, y = 0.25, 
           label = paste("Deviance=",round(pa_monsoonp$deviance,2),
                         "\nAIC =",round(pa_monsoonp$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "Monsoon Precip", y = "Probability of occurrence")


pres_plot_winterp<-ggplot(data=subset(pres_binned,model=="winterp"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  annotate("label", x = 100, y = 0.75, 
           label = paste("Deviance=",round(pa_winterp$deviance,2),
                         "\nAIC =",round(pa_winterp$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "Winter Precip", y = "Probability of occurrence")


pres_plot_springt<-ggplot(data=subset(pres_binned,model=="springt"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  annotate("label", x = 0, y = 1, 
           label = paste("Deviance=",round(pa_springt$deviance,2),
                         "\nAIC =",round(pa_springt$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "Spring Temp", y = "Probability of occurrence")


pres_plot_fallt<-ggplot(data=subset(pres_binned,model=="fallt"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  annotate("label", x = 0, y = 1, 
           label = paste("Deviance=",round(pa_fallt$deviance,2),
                         "\nAIC =",round(pa_fallt$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "Fall Temp", y = "Probability of occurrence")

pres_plot_map<-ggplot(data=subset(pres_binned,model=="map"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  annotate("label", x = 0, y = 1, 
           label = paste("Deviance=",round(pa_map$deviance,2),
                         "\nAIC =",round(pa_map$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "MAP", y = "Probability of occurrence")


pres_plot_mat<-ggplot(data=subset(pres_binned,model=="mat"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  annotate("label", x = 15, y = 1, 
           label = paste("Deviance=",round(pa_mat$deviance,2),
                         "\nAIC =",round(pa_mat$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "MAT", y = "Probability of occurrence")
