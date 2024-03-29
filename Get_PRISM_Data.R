#### Use PRISM climate data to create rasters of monthly climate variables in PIED study region
# from emily's code:https://github.com/emilylschultz/DemographicRangeModel/blob/master/Code/ClimateProcessing/historic.R
library(raster)
library(reshape2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(here)
library(rFIA)
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


temp.max.path <-  "./PRISM_data/PRISM_tmax_stable_4kmM3_189501_198012_bil/"
temp.max.path.new <-  "./PRISM_data/PRISM_tmax_stable_4kmM3_198101_201910_bil/"

temp.maxFiles.old <- list.files(path = temp.max.path, pattern = glob2rx("*tmax*.bil"), full.names = TRUE)
temp.maxFiles.new <- list.files(path = temp.max.path.new, pattern = glob2rx("*tmax*.bil"), full.names = TRUE)
temp.maxFiles<- c(temp.maxFiles.old, temp.maxFiles.new)

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


tmpmaxStack <- stack()
for (i in temp.maxFiles) {
  print(i)
  tmpmaxStack <- stack(tmpmaxStack, raster(i))
}


vpdStack <- stack()
for (i in 1:length(vpdmaxFiles)) {
  print(i)
  rast<-raster(vpdmaxFiles[i])
  if(i == 431 | i == 432){crs(rast)<-crs(raster(vpdmaxFiles[1]))}
  vpdStack <- stack(vpdStack, rast)
}


# extract by plot lat and lon:


PIED.data.az <- read.csv(here("data.clim copy.csv")) # the Arizona PIED data
PIED.region.ll <- read.csv(here("locs-env-1-31-17.csv"))

PIED.region.ll %>% group_by(INV, SPCD, STATECD) %>% summarise(n())

# the locs contain lat long, aspect, slope, species code, 
PIED.region.rwl <- read.csv(here("trees-rwl-1-31-17.csv")) # note that this data is has all RWLS in columsn and a year column


# new mexico data:
NM.rwl<- read.delim(here("new-mexico-ring-width.txt"), sep = ",")
NM.meta<- read.delim(here("new-mexico-meta.txt"), sep = ",")

# connect by CN column (core cn)
NM.all <- left_join(NM.rwl, NM.meta, by = "CN")
unique(NM.meta$SPCD)

NM.meta %>% group_by(SPCD) %>% summarise(n())


NM.meta %>% filter(SPCD ==106)%>%group_by( VERIFY) %>% summarise (n())

# plot up the raw data
filterd.df <- NM.all %>% filter(SPCD ==106)%>% mutate(growth2=ifelse(RW==0,0.001,RW),loggrowth=log(RW))%>% filter(loggrowth <= -2.5)

filtered.df2<- NM.all %>% filter(SPCD ==106)%>% mutate(growth2=ifelse(RW==0,0.001,RW),loggrowth=log(RW))

# if(!exists("/Users/kah/Documents/docker_pecan/pecan/InWeUS_FIA/NM_COND.csv")){
# fiadb <- getFIA(states = "NM", dir = "InWeUS_FIA", common = FALSE, tables = c("PLOT", "TREE", "COND", "SUBPLOT"), nCores = 1)
# }else{
# fiadb <- readFIA(dir = "InWeUS_FIA")
# }


NM.PLOT <- read.csv("/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/NM_PLOT.csv")


# join with lat long
NM.PLOT.small <- unique(NM.PLOT[,c("PLOT", "LAT", "LON", "ELEV", "STATECD", "COUNTYCD")])

NM.meta.ll<- left_join( NM.all,NM.PLOT.small, by = c("COUNTYCD", "STATECD", "PLOT"))

summary(NM.meta.ll %>% filter(!is.na(LON)))
# make a map of all of these:
all_states <- map_data("state")
states <- subset(all_states, region %in% c( "arizona", "utah", "new mexico", "colorado", 
                                            "idaho", "wyoming", "montana", "nevada", 
                                            "california", "oregon", "washington", "texas", "kansas", 
                                            "nebraska", "north dakota", "south dakota") )
coordinates(states)<-~long+lat
class(states)
proj4string(states) <-CRS("+proj=longlat +datum=NAD83")
mapdata <- states
mapdata<-data.frame(mapdata)

PIED.region.ll$SPCD <- as.character(PIED.region.ll$SPCD)

png(height = 6, width = 7, units = "in", res = 200, "new_regional_cores_by_SPCD.png")
ggplot(data = PIED.region.ll, aes(x = LON, y = LAT, color = SPCD))+geom_point()+
  geom_polygon(data=data.frame(mapdata), aes(x=long, y=lat, group=group),
               colour = "darkgrey", fill = NA)+theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), legend.position = "bottom")+ylab("Latitude")+
  xlab("Longitude")+coord_cartesian(xlim = c(-118, -103), ylim = c(32, 49))+facet_wrap(~SPCD)
dev.off()


PIED.rwl.m <- melt(PIED.region.rwl %>% dplyr::select(-X), id.vars = "Year")
colnames(PIED.rwl.m) <- c("year", "series", "growth")
PIED.region.ll$series <- as.character(PIED.region.ll$series)


all.outside.AZ <- left_join(PIED.region.ll,PIED.rwl.m, by = c( "series"))

# get only the PIED
PIED.NM <- NM.meta.ll %>% filter(SPCD %in% "106")
PIED.outside.AZ <- all.outside.AZ %>% filter(SPCD %in% "106")
PIPO.NM <- NM.meta.ll %>% filter(SPCD %in% "122")
PIPO.outside.AZ <- all.outside.AZ %>% filter(SPCD %in% "122")


# just get all the Lat long datasets for the 
# join with the current PIED data
PIED.data.az$dataset <- "AZ"
PIED.outside.AZ$dataset <- "outside AZ"

PIED.az.ll <- unique(PIED.data.az[,c("LAT", "LON", "name", "dataset")])
PIED.outside.ll <- unique(PIED.outside.AZ[,c("LAT", "LON", "series", "dataset")])
colnames(PIED.outside.ll)[3] <- "name"
PIED.nm.ll <- unique(PIED.NM[,c("LAT", "LON", "CN")])
PIED.nm.ll$dataset <- "NM"
colnames(PIED.nm.ll)[3] <- "name"
PIED.ll <- rbind(PIED.az.ll, PIED.outside.ll, PIED.nm.ll)

PIED.outside.AZ %>% filter(year >1999 & !is.na(growth))
max(unique(PIED.outside.AZ$year)) # 2014
max(unique(PIED.NM$Year))

length(unique(PIED.ll$name))
#1818
#length(unique(PIED.ll$))


ppt.extracted <- data.frame(extract(pptStack, PIED.ll [,c("LON", "LAT")]))
ppt.extracted$lon <- PIED.ll$LON
ppt.extracted$lat <- PIED.ll$LAT
ppt.extracted$name<- PIED.ll$name


saveRDS(ppt.extracted, "pied_extracted.ppt.data_v4.rds")

ppt.extracted <- readRDS(  "/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_extracted.ppt.data_v4.rds")


temp.extracted <- data.frame(extract(tmpStack, PIED.ll [,c("LON", "LAT")]))
temp.extracted$lon <- PIED.ll$LON
temp.extracted$lat <- PIED.ll$LAT
temp.extracted$name<- PIED.ll$name


saveRDS(temp.extracted, here("FIA_inc_data","pied_extracted.temp.data_v4.rds"))
temp.extracted <- readRDS( "/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_extracted.temp.data_v4.rds")

tempmax.extracted <- data.frame(raster::extract(tmpmaxStack, PIED.ll [,c("LON", "LAT")]))
tempmax.extracted$lon <- PIED.ll$LON
tempmax.extracted$lat <- PIED.ll$LAT
tempmax.extracted$name<- PIED.ll$name


saveRDS(tempmax.extracted, here("FIA_inc_data","pied_extracted.tmax.data_v4.rds"))
tmax.extracted <- readRDS( "/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_extracted.tmax.data_v4.rds")

# now summarise the climate data to the variables we are interests in:

colnames(tmax.extracted) 
years <- rep(1895:2019, each = 12)
months <- rep(1:12, 125)

colnames.clim <- c(paste0(years, "_", months))
colnames.clim <- colnames.clim[1:1498]
length(tmax.extracted)
colnames(tmax.extracted)[1:1498] <- colnames.clim
#colnames(ppt.extracted)[1:1502] <- colnames.clim
tmax.m <- melt(tmax.extracted, id.vars = c("lon", "lat", "name"))

year.df <- data.frame(variable = c(unique(tmax.m$variable), "2019_11", "2019_12"), 
           year = rep(1895:2019, each = 12), 
           month = rep(1:12, 125))
tmax.long <- left_join(year.df, tmax.m, by = "variable")
tmax.long <- tmax.long[,c("lon", "lat", "name", "year","month", "value")]
colnames(tmax.long) <- c("lon", "lat", "name", "year", "month", "Tave")



colnames(temp.extracted) 
years <- rep(1895:2019, each = 12)
months <- rep(1:12, 125)

colnames.clim <- c(paste0(years, "_", months), "2020_1", "2020_2")
length(temp.extracted)
colnames(temp.extracted)[1:1502] <- colnames.clim
#colnames(ppt.extracted)[1:1502] <- colnames.clim
temp.m <- melt(temp.extracted, id.vars = c("lon", "lat", "name"))

year.df <- data.frame(variable = unique(temp.m$variable)[1:1500], 
                      year = rep(1895:2019, each = 12), 
                      month = rep(1:12, 125))
temp.long <- left_join(year.df, temp.m, by = "variable")
temp.long <- temp.long[,c("lon", "lat", "name", "year","month", "value")]

#temp.long <- temp.m %>% separate(variable, c("year", "month"))
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
#ppt.long <- ppt.m %>% separate(variable, c("year", "month"))
year.df <- data.frame(variable = unique(temp.m$variable)[1:1498], 
                      year = rep(1895:2019, each = 12)[1:1498], 
                      month = rep(1:12, 125)[1:1498])
ppt.long <- left_join(year.df, ppt.m, by = "variable")
ppt.long <- ppt.long[,c("lon", "lat", "name", "year","month", "value")]


colnames(ppt.long) <- c("lon", "lat", "name", "year", "month", "ppt")

temp.ppt <- left_join(temp.long, ppt.long, by =c("lon", "lat", "name", "year", "month"))

write.csv(temp.ppt, "temp.ppt.pied.csv", row.names = FALSE)

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
foresummer <- temp.ppt %>% filter(month %in% 3:5)%>% group_by(lon, lat, name, year)%>% summarise(Tmean_MarAprMay = mean(Tave),
                                                                                                  Precip_MarAprMay = mean(ppt)) 

foresummer$year <- as.character(foresummer$year)
# Jul August - monsson
monsoon <- temp.ppt %>% filter(month %in% 7:8)%>% group_by(lon, lat, name, year)%>% summarise(Tmean_JulAug = mean(Tave),
                                                                                              Precip_JulAug = mean(ppt)) 

monsoon$year <- as.character(monsoon$year)
# Sept Oct November - Fall
fall <-temp.ppt %>% filter(month %in% 9:11)%>% group_by(lon, lat, name, year)%>% summarise(Tmean_SepOctNov = mean(Tave),
                                                                                           Precip_SepOctNov = mean(ppt)) 

fall$year <- as.character(fall$year)
# december january february - winter:

winter <- temp.ppt %>% filter(month %in% c(1,2,12)) %>% group_by(lon, lat, name, water_year) %>% summarise(Tmean_DecJanFeb = mean(Tave),
                                                                                                           Precip_DecJanFeb = mean(ppt)) 
colnames(winter)[4] <- "year"

winter$year <- as.character(winter$year)

# other requested variables:
# november - march
nov_mar <- temp.ppt %>% filter(month %in% c(1,2,3,11,12)) %>% group_by(lon, lat, name, water_year) %>% summarise(Tmean_NovDecJanFebMar = mean(Tave),
                                                                                                                 Precip_NovDecJanFebMar = mean(ppt)) 
colnames(nov_mar)[4] <- "year"

nov_mar$year <- as.character(nov_mar$year)

# april may june
apr_may_jun <- temp.ppt %>% filter(month %in% c(4,5,6)) %>% group_by(lon, lat, name, water_year) %>% summarise(Tmean_AprMayJun  = mean(Tave),
                                                                                                               Precip_AprMayJun = mean(ppt)) 
colnames(apr_may_jun )[4] <- "year"

apr_may_jun $year <- as.character(apr_may_jun $year)

# September to October of current year

sep_oct <- temp.ppt %>% filter(month %in% c(9,10)) %>% group_by(lon, lat, name, water_year) %>% summarise(Tmean_SepOct  = mean(Tave),
                                                                                                    Precip_SepOct = mean(ppt, na.rm=TRUE)) 
colnames(sep_oct)[4] <- "year"

sep_oct$year <- as.character(sep_oct$year)

total_avg_wateryr <- temp.ppt  %>% group_by(lon, lat, name, water_year) %>% summarise(Tmean  = mean(Tave),
                                                                                      wateryr_precip = sum(ppt, na.rm=TRUE)) 
colnames(total_avg_wateryr)[4] <- "year"

total_avg_wateryr$year <- as.character(total_avg_wateryr$year)

# monthly of current year

temp.wide <- temp.ppt %>% ungroup()%>% 
  dplyr::select(lat, lon, name, year, Tave, month)  %>%
  group_by(lon, lat, name)  %>% spread (month, Tave)                         

ppt.wide <- temp.ppt  %>% 
  dplyr::select(lat, lon, name, year, ppt, month)  %>%
  group_by(lon, lat, name)  %>% spread (month, ppt) 

colnames(ppt.wide)[5:16] <- paste0("PPT_", colnames(ppt.wide)[5:16])
colnames(temp.wide)[5:16] <- paste0("TMEAN_", colnames(temp.wide)[5:16])

# get the previous year's temp and precipt
temp.wide.prev <- temp.wide
temp.wide.prev$nextyear <- as.numeric(temp.wide.prev$year)+1                        
temp.wide.prev$actual.year <- temp.wide.prev$year
temp.wide.prev$year <- temp.wide.prev$nextyear

colnames(temp.wide.prev)[5:16] <- paste0("PREV_",colnames(temp.wide.prev)[5:16])



temp.previous <- temp.wide.prev %>% dplyr::select(-nextyear, -actual.year)
temp.previous$year <- as.character(temp.previous$year)
temp.wide$year <- as.character(temp.wide$year)


ppt.wide.prev <- ppt.wide
ppt.wide.prev$nextyear <- as.numeric(ppt.wide.prev$year)+1                        
ppt.wide.prev$actual.year <- ppt.wide.prev$year
ppt.wide.prev$year <- ppt.wide.prev$nextyear

colnames(ppt.wide.prev)[5:16] <- paste0("PREV_",colnames(ppt.wide.prev)[5:16])
ppt.previous <- ppt.wide.prev %>% dplyr::select(-nextyear, -actual.year)
ppt.previous$year <- as.character(ppt.previous$year)
ppt.wide$year <- as.character(ppt.wide$year)

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
all.tmean.8 <- left_join(all.tmean.7, total_avg_wateryr)
full.ppt.tmean <- left_join(all.tmean.8, temp.previous, by = mergeCols)



#ggplot(full.ppt.tmean, aes(lon, lat, color = PREV_TMEAN_PPT_6))+geom_point()


full.ppt.tmean$tmp_yr <- rowMeans(full.ppt.tmean[,c("TMEAN_1","TMEAN_10","TMEAN_11",            
                                                    "TMEAN_12","TMEAN_2","TMEAN_3",               
                                                    "TMEAN_4","TMEAN_5","TMEAN_6",               
                                                    "TMEAN_7","TMEAN_8", "TMEAN_9")])


full.ppt.tmean$ppt_yr <- rowSums(full.ppt.tmean[,c("PPT_1","PPT_10","PPT_11",            
                                                    "PPT_12","PPT_2","PPT_3",               
                                                    "PPT_4","PPT_5","PPT_6",               
                                                    "PPT_7","PPT_8", "PPT_9")])

full.ppt.tmean$ppt_yr_mean <- rowMeans(full.ppt.tmean[,c("PPT_1","PPT_10","PPT_11",            
                                                   "PPT_12","PPT_2","PPT_3",               
                                                   "PPT_4","PPT_5","PPT_6",               
                                                   "PPT_7","PPT_8", "PPT_9")])


climatenormals <- full.ppt.tmean %>% group_by(lat, lon, name) %>% summarise(ppt_norm = mean(ppt_yr, na.rm = TRUE), 
                                                                            tmp_norm = mean(tmp_yr, na.rm = TRUE),
                                                                            ppt_norm_mean = mean(ppt_yr_mean, na.rm =TRUE),
                                                                            ppt_norm_wateryear = mean(wateryr_precip, na.rm =TRUE),
                                                                            temp_norm_wateryear = mean(Tmean, na.rm =TRUE)) 



ggplot(climatenormals, aes(ppt_norm, ppt_norm_wateryear))+geom_point()+geom_abline(aes(slope = 1, intercept = 0))
ggplot(climatenormals, aes(ppt_norm, ppt_norm_mean))+geom_point()+geom_abline(aes(slope = 1, intercept = 0))

full.ppt.tmean.norms <- right_join(full.ppt.tmean, climatenormals, by = c("lon", "lat", "name"))



write.csv(full.ppt.tmean.norms, "pied_all_tmean_ppt_v7.csv", row.names = FALSE)



# connect with the rest of the dataset:

full.ppt.tmean.norms <- read.csv(here("FIA_inc_data","pied_all_tmean_ppt_v7.csv"))

head(PIED.data.az)
head(PIED.nm.ll)

head(PIED.outside.AZ)
unique(PIED.outside.AZ$TREE)

PIED.outside.AZ$name <- PIED.outside.AZ$series
head(PIED.outside.AZ)
#PIED.outside.AZ$name <- paste0(PIED.outside.AZ$COUNTYCD, "0", PIED.outside.AZ$PLOT, "0", PIED.outside.AZ$SUBP, "0", PIED.outside.AZ$STATECD)
# arizona name is the county, a 0, plot, a 0 and su


colnames(PIED.outside.AZ) 
colnames(PIED.data.az)
PIED.outside.AZ$DIA_measured <- as.numeric(PIED.outside.AZ$DIA) *2.54 # get diameter in cm

unique(PIED.outside.AZ$name)
PIED.outside.AZ.nona <- PIED.outside.AZ[!is.na(PIED.outside.AZ$growth),]
# get the previous year's diameter from tree rings:
df.list <- list()

# this isnt working because for some names we have multiple measurements of grwoth for that year:

for(i in unique(PIED.outside.AZ.nona$name)){
  
  #get_prev<- function(i){
  df <- PIED.outside.AZ.nona[PIED.outside.AZ.nona$name %in% i, ]
  df$DIA_cm <- NA
  df$DIA_prev <- NA
  
  if(length(df$X > 1)){
    for(t in max(df$year, na.rm = TRUE):min(df$year, na.rm = TRUE)){
      if(df[df$year == t,]$year  == max(df$year, na.rm = TRUE)){
        df[df$year == t,]$DIA_cm <- df[df$year == t,]$DIA_measured 
      }else{
        df[df$year == t,]$DIA_cm <- df[df$year == (t+1),]$DIA_cm - (((df[df$year == t,]$growth)/10)*2) #prev diameter = dia in cm - (growth (mm)/(10mm/1cm) *2)
        df[df$year == t,]$DIA_prev <- df[df$year == (t+1),]$DIA_cm 
      }
    }
  }
  df.list[[i]] <- df
}

#df.list <- lapply(unique(PIED.outside.AZ$name), get_prev)

PIED.outside.prev_DBH <- do.call(rbind, df.list)

summary(PIED.outside.prev_DBH$DIA_cm)


##################################################################################
# Make sure we are getting the correct DBH values for the PIED AZ data
##################################################################################
PIED.data.az$DIA_measured <- as.numeric(PIED.data.az$DIA) *2.54 # get diameter in cm

unique(PIED.data.az$name)

# Back calculating tree diameter tree rings:
df.list <- list()

# this isnt working because for some names we have multiple measurements of grwoth for that year for some pied trees:
# in AZ there are 62 trees PIED trees that have multiple values of growth per year (out of 520 trees total)
# lets remove these here, but leaving a note to look at these values in tellervo
dups2remove <- unique(PIED.data.az[duplicated(PIED.data.az[,c("name", "year")]),]$name)
PIED.data.az.nodup <- filter(PIED.data.az, !name  %in% dups2remove )

# make sure we remove all the extra years with NA values (most of these are from rwl formatting)
PIED.data.az.nodup.nona <- PIED.data.az.nodup[!is.na(PIED.data.az.nodup$growth),]

# messy for loop to calculate the previous years diameter and add to df
for(i in 1:length(unique(PIED.data.az.nodup.nona$name))){
  print(i)
  # select a df with the tree of interest i
  df <- PIED.data.az.nodup.nona[PIED.data.az.nodup.nona$name %in% unique(PIED.data.az.nodup.nona$name)[i], ]
  
  df$DIA_cm <- NA # make a DIA_cm column
  df$DIA_prev <- NA # make a DIA_cm column
  
  if(length(df$X > 1)){ # if the df has values in it..
    for(t in max(df$year, na.rm = TRUE):min(df$year, na.rm = TRUE)){ 
      if(df[df$year == t,]$year  == max(df$year, na.rm = TRUE)){# if t = the last year of rw records, set the DIA_cm = measured DIA (which should already be in cm)
        df[df$year == t,]$DIA_cm <- df[df$year == t,]$DIA_measured 
      }else{ # else set the DIA_cm as the previous years diameter minus 2*tree ring growth. Since tr growth is in mm, need to also convert to cm
        df[df$year == t,]$DIA_cm <- df[df$year == (t+1),]$DIA_cm - (((df[df$year == t,]$growth)/10)*2) # diameter = dia in cm - (growth (mm)/(10mm/1cm) *2)
        df[df$year == t,]$DIA_prev <- df[df$year == (t+1),]$DIA_cm 
      }
    }
  }
  df.list[[i]] <- df
}

#df.list <- lapply(unique(PIED.data.az.nodup$name), get_prev)

PIED.arizona.prev_DBH <- do.call(rbind, df.list)

summary(PIED.arizona.prev_DBH $DIA_cm)
hist(PIED.arizona.prev_DBH $DIA_cm)
hist(PIED.arizona.prev_DBH$DIA_prev)

hist(PIED.outside.prev_DBH$DIA_prev)

PIED.az <- PIED.arizona.prev_DBH  %>% dplyr::select(name, year, LAT, LON, X, growth, PLOT.x, SUBP, ELEV, DIA, SLOPE, ASPECT,DIA, DIA_prev,dataset)
PIED.az$STATECD <- 4# add STATECD for AZ

PIED.other <- PIED.outside.prev_DBH %>% dplyr::select(name, year, LAT, LON, X, growth,  PLOT, SUBP, ELEV,  SLOPE, ASPECT, DIA_measured, DIA_prev, dataset, STATECD)
# fix some column naming between PIED and PIPO
colnames(PIED.az)[7]<- "PLOT"
colnames(PIED.other)[12:13] <- c("DIA", "DIA_prev")



##################################################################################
# Estimate diameter for PIED NM data
##################################################################################

PIED.NM$DIA_measured <- as.numeric(PIED.NM$DIA) *2.54 # get diameter in cm
PIED.NM$name <- PIED.NM$CN
PIED.NM$year <- PIED.NM$Year
PIED.NM$growth <- PIED.NM$RW
unique(PIED.NM$name)

# Back calculating tree diameter tree rings:
df.list <- list()

# this isnt working because for some names we have multiple measurements of grwoth for that year for some pied trees:
# in AZ there are 62 trees PIED trees that have multiple values of growth per year (out of 520 trees total)
# lets remove these here, but leaving a note to look at these values in tellervo
dups2remove <- unique(PIED.NM[duplicated(PIED.NM[,c("name", "year")]),]$name)
PIED.NM.nodup <- filter(PIED.NM, !name  %in% dups2remove )

# make sure we remove all the extra years with NA values (most of these are from rwl formatting)
PIED.NM.nodup.nona <- PIED.NM.nodup[!is.na(PIED.NM.nodup$growth),]


# messy for loop to calculate the previous years diameter and add to df
for(i in 1:length(unique(PIED.NM.nodup.nona$name))){
  print(i)
  #get_prev<- function(i){
  df <- PIED.NM.nodup.nona[PIED.NM.nodup.nona$name %in% unique(PIED.NM.nodup.nona$name)[i], ]
  df$DIA_cm <- NA
  df$DIA_prev <- NA
  
  if(length(df$year > 1)){
    for(t in max(df$year, na.rm = TRUE):min(df$year, na.rm = TRUE)){
      if(df[df$year == t,]$year  == max(df$year, na.rm = TRUE)){
        df[df$year == t,]$DIA_cm <- df[df$year == t,]$DIA_measured 
      }else{
        df[df$year == t,]$DIA_cm <- df[df$year == (t+1),]$DIA_cm - (((df[df$year == t,]$growth)/10)*2) #prev diameter = dia in cm - (growth (mm)/(10mm/1cm) *2)
        df[df$year == t,]$DIA_prev <- df[df$year == (t+1),]$DIA_cm 
      }
    }
  }
  df.list[[i]] <- df
}

#df.list <- lapply(unique(PIED.NM.nodup$name), get_prev)

PIED.NM.prev_DBH <- do.call(rbind, df.list)

summary(PIED.NM.prev_DBH $DIA_cm)
hist(PIED.NM.prev_DBH  $DIA_cm)
hist(PIED.NM.prev_DBH $DIA_prev)

hist(PIED.outside.prev_DBH$DIA_prev)
PIED.NM.prev_DBH $dataset <- "New Mexico"

PIED.nm <- PIED.NM.prev_DBH  %>% dplyr::select(name, year, LAT, LON, growth, PLOT, SUBP, ELEV, DIA,  DIA, DIA_prev,dataset, STATECD)

PIED.az <- PIED.arizona.prev_DBH  %>% dplyr::select(name, year, LAT, LON,  growth, PLOT.x, SUBP, ELEV, DIA,DIA, DIA_prev,dataset)
PIED.az$STATECD <- 4# add STATECD for AZ

PIED.other <- PIED.outside.prev_DBH %>% dplyr::select(name, year, LAT, LON, growth,  PLOT, SUBP, ELEV,   DIA_measured, DIA_prev, dataset, STATECD)
# fix some column naming between PIED and PIPO
colnames(PIED.az)[6]<- "PLOT"
colnames(PIED.other)[9:10] <- c("DIA", "DIA_prev")

unique(PIED.nm$name)
unique(PIED.az$name)
unique(PIED.other$name)

colnames(PIED.az)
colnames(PIED.other)
colnames(PIED.nm)
### combine all

PIED.all <- rbind(PIED.az, PIED.other, PIED.nm)

summary(unique(PIED.all$name) %in% unique(PIED.ll$name))
class(PIED.all$name)

# seven trees have negative diameter for at least 1 year--leading to 147 years in the dataset with negetive diameters
neg.dbh <- filter(PIED.all, DIA_prev <= 0)
unique(neg.dbh[,c("name", "year")])

# lets remove these
PIED.all.pos <- filter(PIED.all, DIA_prev > 0)
unique(PIED.all.pos$name)
PIED.all.pos$ST_PLT <- paste0(PIED.all.pos$STATECD, "_", PIED.all.pos$PLOT)


write.csv(PIED.all.pos, here("FIA_inc_data","pied_all_growth_v6.csv"))

hist(PIED.all$DIA_prev)

PIED.all <- read.csv("/Users/kah/Documents/docker_pecan/pecan/FIA_inc_data/pied_all_growth_v6.csv")


ggplot(PIED.all, aes(growth))+geom_histogram()+facet_wrap(~dataset)
full.ppt.tmean.norms <- read.csv(here("pied_all_tmean_ppt_v7.csv"))

new.merged.growth <- left_join(PIED.all, full.ppt.tmean.norms, by.x = c("name", "year", "LON", "LAT"), by.y = c("name", "year", "lon", "lat"))

# we lose 75 trees
# missing Lat Lons for thosw
lost <- PIED.all %>% filter(is.na(LON) | is.na(LAT) | is.na(growth) | is.na(DIA_prev))
unique(lost$name)

length(unique(PIED.all$name))
length(unique(new.merged.growth$ST_PLT))

grow.na.omit <- na.omit(new.merged.growth)
length(unique(grow.na.omit$ST_PLT))
length(unique(grow.na.omit$PLOT))
length(unique(grow.na.omit$name))

write.csv(grow.na.omit, "PIED_grow_na_omit.csv", row.names = FALSE)

df<- PIED.all %>% group_by(name) %>% summarise(minyear = min(year), 
                                               maxyear = max(year))

min(df$minyear)
max(df$minyear)
max(df$maxyear)

grow.monsoon <- na.omit(new.merged.growth) %>% 
  mutate_at(scale, .vars = vars(ppt_norm, tmp_norm, DIA_prev))%>%
  arrange(STATECD, PLOT,SUBP,name) %>%
  
  mutate(PlotCD=as.numeric(factor(ST_PLT, levels = unique(ST_PLT))),treeCD=as.numeric(factor(name,levels=unique(name))),
         growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2))%>%
  group_by(PlotCD)%>%
  mutate_at(scale, .vars = vars(Precip_JulAug, Precip_NovDecJanFebMar, Tmean_AprMayJun, Tmean_SepOct)) %>%
  ungroup()  

ggplot(grow.monsoon, aes(loggrowth))+geom_histogram()+facet_wrap(~dataset)

# look at the number of trees per plot
# look 

summary(grow.monsoon)
grow.monsoon %>% group_by(STATECD) %>% summarise(n_distinct(name))
#STATECD `n_distinct(name)`
# <int>              <int>
# Arizona     4                457
# Colorado       8                 43
# New Mexico      35                754
# Utah      49                304



length(unique(grow.monsoon$name))
length(unique(grow.monsoon$PlotCD))
length(unique(PIED.all$name))

summary(grow.monsoon$DIA_prev)
summary(grow.monsoon$LON)
summary(PIED.all$LON)
summary(full.ppt.tmean.norms$lon)

nona.clim <- full.ppt.tmean.norms[!is.na(full.ppt.tmean.norms$lon),]
nona.PIED <- PIED.all[!is.na(PIED.all$LON),]

ggplot(PIED.all, aes(x = LON, y = LAT, color = STATECD))+geom_point()

