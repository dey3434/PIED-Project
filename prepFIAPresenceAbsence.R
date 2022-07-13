library(raster)

# Read in data
all <- read.csv("./FIAdata/TREE_COMBINED.csv", header = T, stringsAsFactors = F)

# Keep only plot and species code
all <- all[, c("PLT_CN", "SPCD")]

# Read in plots
plots <- read.csv("./FIAdata/PLOT_COMBINED.csv")
plots <- plots[, c("CN", "LAT", "LON")]

# Processing
all$LAT <- NA
all$LON <- NA
progress <- 1
no <- length(plots$CN)
for (i in plots$CN) {
  if (nrow(all[all$PLT_CN == i,]) > 0) {
    all[all$PLT_CN == i,]$LAT <- subset(plots, CN == i)$LAT
    all[all$PLT_CN == i,]$LON <- subset(plots, CN == i)$LON
  }
  progress <- progress + 1
  print(progress / no * 100)
}

# Export
write.csv(all, "./Processed/Validation/validationRecords.csv", row.names = F)

# Read in a raster for alignment
pres <- raster("./Output/tifs/PIED.clim_lambda_gam.tif")
pres <- setValues(pres, NA)

# Make records spatial
proj <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0"
allSp <- SpatialPointsDataFrame(all[, c("LON", "LAT")], as.data.frame(all[, "SPCD"]), proj4string = CRS("+proj=longlat +datum=NAD83"))

# For each tree, find in which raster cell it falls
allSp$cellNumber <- NA
progress <- 1
no <- nrow(allSp)
for (i in 1:nrow(allSp)) {
 tmp <- as.integer(extract(pres, allSp[i,], cellnumbers = T)[[1]])
 allSp[i, "cellNumber"] <- tmp
 progress <- progress + 1
 print(progress / no * 100)
} 

# For each raster cell, find whether it contains PIED (SPCD = 106)
names(allSp) <- c("SPCD", "cellNumber")
PIED <- 106
for (i in 1:ncell(pres)) {
  tmp <- subset(allSp, cellNumber == i)
  if (nrow(tmp) > 0) {
    if (PIED %in% unique(tmp$SPCD)) pres[i] <- 1
    else pres[i] <- 0
  } else pres[i] <- 0
}

# Export presence/absence raster
writeRaster(pres, "./Processed/Validation/presenceAbsenceRaster2.tif", overwrite = T)

climate_wna<-plots %>%
  group_by(PLOT) %>%
  summarise(lat=mean(LAT,na.rm=T),long=mean(LON,na.rm=T),el=(mean(ELEV,na.rm=T)*0.3048)) %>%
  mutate(ID2=".")

names(climate_wna)[1]<-"ID1"

write.csv(climate_wna[1:4215,],"./ClimateData/FIA_climateWNA1.csv",row.names=F)
write.csv(climate_wna[4216:8430,],"./ClimateData/FIA_climateWNA2.csv",row.names=F)
write.csv(climate_wna[8431:12645,],"./ClimateData/FIA_climateWNA3.csv",row.names=F)
write.csv(climate_wna[12646:16863,],"./ClimateData/FIA_climateWNA4.csv",row.names=F)

