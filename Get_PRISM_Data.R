#### Use PRISM climate data to create rasters of monthly climate variables in PIED study region
# from emily's code:https://github.com/emilylschultz/DemographicRangeModel/blob/master/Code/ClimateProcessing/historic.R
library(raster)

### PRISM download January 22, 2019
### January 1981 through June 2018
### (37*12) + 6 = 450 files

# Search for PRISM files
PRISM.path <-  "./PRISM_data/"
ppt.path <-  "./PRISM_data/PRISM_ppt_stable_4kmM2_189501_198012_bil/"
ppt.path.new <-  "./PRISM_data/PRISM_ppt_stable_4kmM3_198101_201910_bil/"


pptFiles.old <- list.files(path = ppt.path, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
pptFiles.new <- list.files(path = ppt.path.new, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
pptFiles<- c(pptFiles.old, pptFiles.new)

temp.path <-  "./PRISM_data/PRISM_tmean_stable_4kmM3_189501_198012_bil/"
temp.path.new <-  "./PRISM_data/PRISM_tmean_stable_4kmM3_198101_202002_bil/"

tempFiles.old <- list.files(path = temp.path, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
tempFiles.new <- list.files(path = temp.path.new, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
tempFiles<- c(tempFiles.old, tempFiles.new)

#tmpFiles <- list.files(path = PRISM.path, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)

#vpdminFiles <- list.files(path = PRISM.path, pattern = glob2rx("*vpdmin*.bil"), full.names = TRUE)
vpdmaxFiles <- list.files(path = PRISM.path, pattern = glob2rx("*vpdmax*.bil"), full.names = TRUE)

# Stack monthly data
pptStack <- stack()
for (i in pptFiles) {
  print(i)
  pptStack <- stack(pptStack, raster(i))
}

tmpStack <- stack()
for (i in tempFiles) {
  print(i)
  tmpStack <- stack(tmpStack, raster(i))
}

vpdStack <- stack()
for (i in 1:length(vpdmaxFiles)) {
  print(i)
  rast<-raster(vpdmaxFiles[i])
  if(i == 431 | i == 432){crs(rast)<-crs(raster(vpdmaxFiles[1]))}
  vpdStack <- stack(vpdStack, rast)
}


# extract by plot lat and lon:
PIED.data <- read.csv("data.clim copy.csv")
PIED.ll <- unique(PIED.data[,c("LAT", "LON", "name")])


ppt.extracted <- data.frame(extract(pptStack, PIED.ll [,c("LON", "LAT")]))
ppt.extracted$lon <- PIED.ll$LON
ppt.extracted$lat <- PIED.ll$LAT
ppt.extracted$name<- PIED.ll$name


saveRDS(ppt.extracted, "/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_extracted.ppt.data.rds")

ppt.extracted <- readRDS( "/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_extracted.ppt.data.rds")


temp.extracted <- data.frame(extract(tmpStack, PIED.ll [,c("LON", "LAT")]))
temp.extracted$lon <- PIED.ll$LON
temp.extracted$lat <- PIED.ll$LAT
temp.extracted$name<- PIED.ll$name


saveRDS(temp.extracted, "/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_extracted.temp.data.rds")
temp.extracted <- readRDS( "/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_extracted.temp.data.rds")


# now summarise the climate data to the variables we are interests in:


colnames(temp.extracted) 
years <- rep(1895:2019, each = 12)
months <- rep(1:12, 125)

colnames.clim <- c(paste0(years, "_", months), "2020_1", "2020_2")
length(temp.extracted)
colnames(temp.extracted)[1:1502] <- colnames.clim
#colnames(ppt.extracted)[1:1502] <- colnames.clim
temp.m <- melt(temp.extracted, id.vars = c("lon", "lat", "name"))
temp.long <- temp.m %>% separate(variable, c("year", "month"))
colnames(temp.long) <- c("lon", "lat", "name", "year", "month", "Tave")

# do the same for precipitation
colnames(ppt.extracted) 
years <- rep(1895:2018, each = 12)
months <- rep(1:12, 124)

colnames.clim <- c(paste0(years, "_", months), "2019_1", "2019_2", "2019_3", "2019_4", "2019_5", "2019_6", "2019_7", "2019_8", "2019_9", "2019_10")
length(ppt.extracted)
colnames(ppt.extracted)[1:1498] <- colnames.clim
#colnames(ppt.extracted)[1:1502] <- colnames.clim
ppt.m <- melt(ppt.extracted, id.vars = c("lon", "lat", "name"))
ppt.long <- ppt.m %>% separate(variable, c("year", "month"))
colnames(ppt.long) <- c("lon", "lat", "name", "year", "month", "ppt")

temp.ppt <- left_join(temp.long, ppt.long, by =c("lon", "lat", "name", "year", "month"))


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
foresummer <- temp.ppt %>% filter(month %in% 9:11)%>% group_by(lon, lat, name, year)%>% summarise(Tmean_MarAprMay = mean(Tave),
                                                                                                  Precip_MarAprMay = mean(ppt)) 


# Jul August - monsson
monsoon <- temp.ppt %>% filter(month %in% 7:8)%>% group_by(lon, lat, name, year)%>% summarise(Tmean_JulAug = mean(Tave),
                                                                                              Precip_JulAug = mean(ppt)) 

# Sept Oct November - Fall
fall <-temp.ppt %>% filter(month %in% 9:11)%>% group_by(lon, lat, name, year)%>% summarise(Tmean_SepOctNov = mean(Tave),
                                                                                           Precip_SepOctNov = mean(ppt)) 

# december january february - winter:

winter <- temp.ppt %>% filter(month %in% c(1,2,12)) %>% group_by(lon, lat, name, water_year) %>% summarise(Tmean_DecJanFeb = mean(Tave),
                                                                                                           Precip_DecJanFeb = mean(ppt)) 
colnames(winter)[4] <- "year"

winter$year <- as.character(winter$year)

# other requested variables:
# november - march
nov_mar <- temp.ppt %>% filter(month %in% c(1,2,2,11,12)) %>% group_by(lon, lat, name, water_year) %>% summarise(Tmean_NovDecJanFebMar = mean(Tave),
                                                                                                            Precip_NovDecJanFebMar = mean(ppt)) 
colnames(nov_mar)[4] <- "year"

nov_mar$year <- as.character(nov_mar$year)

# april may june
apr_may_jun <- temp.ppt %>% filter(month %in% c(4,5,6)) %>% group_by(lon, lat, name, water_year) %>% summarise(Tmean_AprMayJun  = mean(Tave),
                                                                                                                 Precip_AprMayJun = mean(ppt)) 
colnames(apr_may_jun )[4] <- "year"

apr_may_jun $year <- as.character(apr_may_jun $year)

# September to October of current year

sep_oct <- temp.ppt %>% filter(month %in% c(9,10)) %>% group_by(lon, lat, name, year) %>% summarise(Tmean_SepOct  = mean(Tave),
                                                                                                               Precip_SepOct = mean(ppt)) 
colnames(sep_oct)[4] <- "year"

sep_oct$year <- as.character(apr_may_jun $year)


# monthly of current year

temp.wide <- temp.ppt %>% filter(! year %in% "T1") %>% 
                           dplyr::select(lat, lon, name, year, Tave, month)  %>%
                            group_by(lon, lat, name)  %>% spread (month, Tave)                         

ppt.wide <- temp.ppt  %>% filter(! year %in% "T1") %>% 
                          dplyr::select(lat, lon, name, year, ppt, month)  %>%
                          group_by(lon, lat, name)  %>% spread (month, ppt) 
                       
colnames(ppt.wide)[5:16] <- paste0("PPT_", colnames(ppt.wide)[5:16])
colnames(temp.wide)[5:16] <- paste0("TMEAN_", colnames(ppt.wide)[5:16])

# get the previous year's temp and precipt
temp.wide.prev <- temp.wide
temp.wide.prev$nextyear <- as.numeric(temp.wide.prev$year)+1                        
temp.wide.prev$actual.year <- temp.wide.prev$year
temp.wide.prev$year <- temp.wide.prev$nextyear

colnames(temp.wide.prev)[5:16] <- paste0("PREV_",colnames(temp.wide.prev)[5:16])



temp.previous <- temp.wide.prev %>% dplyr::select(-nextyear, -actual.year)
temp.previous$year <- as.character(temp.previous$year)

ppt.wide.prev <- ppt.wide
ppt.wide.prev$nextyear <- as.numeric(ppt.wide.prev$year)+1                        
ppt.wide.prev$actual.year <- ppt.wide.prev$year
ppt.wide.prev$year <- ppt.wide.prev$nextyear

colnames(ppt.wide.prev)[5:16] <- paste0("PREV_",colnames(ppt.wide.prev)[5:16])
ppt.previous <- ppt.wide.prev %>% dplyr::select(-nextyear, -actual.year)
ppt.previous$year <- as.character(ppt.previous$year)


# combine all together:
mergeCols = c("lon", "lat", "name", "year")
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



write.csv(full.ppt.tmean, "/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_all_tmean_ppt_v2.csv", row.names = FALSE)



