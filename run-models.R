# Restore the renv environment
# renv::snapshot()    # Used to create the renv environment. Not to be run by the user
renv::restore()       # Used to restore the renv environment to improve reproducibility.

# Setup the data and results folder
library(here) ## make sure the RStudio project is active
if (!dir.exists(here::here("data"))) {
  dir.create(here::here("data"))
}
if (!dir.exists(here::here("results"))) {
  dir.create(here::here("results"))
}
if (!dir.exists(here::here("images"))) {
  dir.create(here::here("images"))
}
for (i in 1:10) {
  if (!dir.exists(here::here("images", paste("model", i-1, sep="_")))) {
    dir.create(here::here("images", paste("model", i-1, sep="_")))
  }
}
# Download the data and save in folder
# if (file.exists(here::here("data", "PIED_data.csv"))) {
#   PIED.all <- read.csv(here::here("data", "PIED_data.csv"))
# } else {
#   PIED.all <- read.csv(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_all_growth_v7.csv"))  
#   write.csv(PIED.all, file = here::here("data", "PIED_data.csv"), row.names = FALSE)
# }
# 
# if (file.exists(here::here("data", "climate_data.csv"))) {
#   full.ppt.tmean.norms <- read.csv(here::here("data", "climate_data.csv"))
# } else {
#   full.ppt.tmean.norms <- read.csv(url("https://data.cyverse.org/dav-anon/iplant/home/smdey/data/pied_all_tmean_ppt_v6.csv"))
#   write.csv(full.ppt.tmean.norms, file = here::here("data", "climate_data.csv"), row.names = FALSE)
# }


# Run the models with RStudioAPI
library(rstudioapi)
rstudioapi::jobRunScript("R/model_0.R")
rstudioapi::jobRunScript("R/model_1.R")
rstudioapi::jobRunScript("R/model_2.R")
rstudioapi::jobRunScript("R/model_3.R")
rstudioapi::jobRunScript("R/model_4.R")
rstudioapi::jobRunScript("R/model_5.R")
rstudioapi::jobRunScript("R/model_6.R")
rstudioapi::jobRunScript("R/model_7.R")
rstudioapi::jobRunScript("R/model_8.R")
rstudioapi::jobRunScript("R/model_9.R")




