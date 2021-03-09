# Read_PIED_data_format.R
# code to read in the Pinus edulis data, and format for our stan model.


#----------------------------------------------------------------------------------------------
#  Read in the data 
#----------------------------------------------------------------------------------------------

PIED.all <- read.csv("data/pied_all_growth_v4.csv") #read in the pied growth data
full.ppt.tmean.norms <- read.csv(url("https://de.cyverse.org/dl/d/329EBCF7-817F-497C-BEA7-A2D906493392/pied_all_tmean_ppt_v3.csv"))
grow.new <- merge(PIED.all, full.ppt.tmean.norms, by.x = c("name", "year", "LON", "LAT"), by.y = c("name", "year", "lon", "lat"))


#----------------------------------------------------------------------------------------------
#  Scale covariates & split into testing and training
#----------------------------------------------------------------------------------------------
# here we are doing global scaling...need to update to local scaling

grow.monsoon<-na.omit(grow.new) %>% 
  mutate_at(scale, .vars = vars(Precip_JulAug, Precip_NovDecJanFebMar, Tmean_AprMayJun, Tmean_SepOct, tmp_norm, ppt_norm)) %>%
  arrange(PLOT,SUBP,name) %>%
  mutate(PlotCD=as.numeric(factor(PLOT, levels = unique(PLOT))),treeCD=as.numeric(factor(name,levels=unique(name))),
         growth2=ifelse(growth==0,0.001,growth),loggrowth=log(growth2))


# spilt into testing and training data
split=0.20
trainIndex <- createDataPartition(grow.monsoon$name, p=split, list=FALSE)
grow_test <- grow.monsoon[trainIndex,]
grow_train <- grow.monsoon[-trainIndex,]

#----------------------------------------------------------------------------------------------
#  Set up data inputs for stan model
#----------------------------------------------------------------------------------------------

# xG and xGtest is big matrices with all the covariates and interactions 
xG<-as.matrix(cbind(grow_train$ppt_norm, grow_train$tmp_norm, grow_train$Precip_JulAug, grow_train$Precip_NovDecJanFebMar, 
                    grow_train$Tmean_AprMayJun, grow_train$Tmean_SepOct, grow_train$DIA_prev,
                    grow_train$ppt_norm*grow_train$tmp_norm, grow_train$ppt_norm*grow_train$Precip_JulAug, 
                    grow_train$ppt_norm*grow_train$DIA_prev, grow_train$ppt_norm*grow_train$Precip_NovDecJanFebMar,
                    grow_train$ppt_norm*grow_train$Tmean_AprMayJun, grow_train$ppt_norm*grow_train$Tmean_SepOct,
                    grow_train$tmp_norm*grow_train$DIA_prev, grow_train$tmp_norm*grow_train$Precip_JulAug,
                    grow_train$tmp_norm*grow_train$Precip_NovDecJanFebMar, grow_train$tmp_norm*grow_train$Tmean_AprMayJun,
                    grow_train$tmp_norm*grow_train$Tmean_SepOct, grow_train$DIA_prev*grow_train$Precip_JulAug,
                    grow_train$DIA_prev*grow_train$Precip_NovDecJanFebMar, grow_train$DIA_prev*grow_train$Tmean_AprMayJun,
                    grow_train$DIA_prev*grow_train$Tmean_SepOct, grow_train$Precip_JulAug*grow_train$Precip_NovDecJanFebMar,
                    grow_train$Precip_JulAug*grow_train$Tmean_AprMayJun, grow_train$Precip_JulAug*grow_train$Tmean_SepOct,
                    grow_train$Precip_NovDecJanFebMar*grow_train$Tmean_AprMayJun, grow_train$Precip_NovDecJanFebMar*grow_train$Tmean_SepOct,
                    grow_train$Tmean_AprMayJun*grow_train$Tmean_SepOct))
xGtest<-as.matrix(cbind(grow_test$ppt_norm, grow_test$tmp_norm, grow_test$Precip_JulAug, grow_test$Precip_NovDecJanFebMar, 
                        grow_test$Tmean_AprMayJun, grow_test$Tmean_SepOct, grow_test$DIA_prev,
                        grow_test$ppt_norm*grow_test$tmp_norm, grow_test$ppt_norm*grow_test$Precip_JulAug, 
                        grow_test$ppt_norm*grow_test$DIA_prev, grow_test$ppt_norm*grow_test$Precip_NovDecJanFebMar,
                        grow_test$ppt_norm*grow_test$Tmean_AprMayJun, grow_test$ppt_norm*grow_test$Tmean_SepOct,
                        grow_test$tmp_norm*grow_test$DIA_prev, grow_test$tmp_norm*grow_test$Precip_JulAug,
                        grow_test$tmp_norm*grow_test$Precip_NovDecJanFebMar, grow_test$tmp_norm*grow_test$Tmean_AprMayJun,
                        grow_test$tmp_norm*grow_test$Tmean_SepOct, grow_test$DIA_prev*grow_test$Precip_JulAug,
                        grow_test$DIA_prev*grow_test$Precip_NovDecJanFebMar, grow_test$DIA_prev*grow_test$Tmean_AprMayJun,
                        grow_test$DIA_prev*grow_test$Tmean_SepOct, grow_test$Precip_JulAug*grow_test$Precip_NovDecJanFebMar,
                        grow_test$Precip_JulAug*grow_test$Tmean_AprMayJun, grow_test$Precip_JulAug*grow_test$Tmean_SepOct,
                        grow_test$Precip_NovDecJanFebMar*grow_test$Tmean_AprMayJun, grow_test$Precip_NovDecJanFebMar*grow_test$Tmean_SepOct,
                        grow_test$Tmean_AprMayJun*grow_test$Tmean_SepOct))

# yG and yGtest are the log transformed growth estimates for each tree & year
yG <- as.vector(grow_train$loggrowth)
yGtest <- as.vector(grow_test$loggrowth)

# setting up indices
nG <- nrow(grow_train)
nGtest<-nrow(grow_test)
plot<-grow_train$PlotCD
nplot<-length(unique(grow_train$PlotCD))
K<-ncol(xG)
tree<-grow_train$treeCD
treetest<-grow_test$treeCD
ntree<-length(unique(grow_train$treeCD))
plotfortree<-grow_train %>%
  group_by(treeCD) %>%
  summarize(Plot=mean(PlotCD))
plotfortree<-plotfortree$Plot


# data list for stan
pied_dat <- list(K = K, nG = nG, nGtest = nGtest, yG = yG, xG = xG, xGtest = xGtest, plot = plot, 
                 nplot = nplot, tree = tree, treetest = treetest, ntree = ntree, 
                 plotfortree = plotfortree)

