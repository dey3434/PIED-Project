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
