library(stringr)
library(tidyr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#_____ Script setup _____#
scripts.loc <- "ARIM-analysis/"
mod.nam <- "mod_int_paths"
development <- FALSE # Set to TRUE for running test model with only develop.spp, and FALSE to run the full model.
mod <- R.utils::loadObject("mod_int_paths")
source(str_c(scripts.loc, "Data_processing_path.R"))
#________________________#

### Parameter summaries ###
## Intermediate path regressions ##
