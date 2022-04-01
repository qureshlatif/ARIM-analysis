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
source(str_c(scripts.loc, "Data_processing_path_flatten_marg.R"))
#________________________#

### Save parameter summaries ###
write.csv(mod$summary, "Mech_pathway_parameters.csv")

## Plot mean mech values over time by development stratum ##
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

# Compile index #
Dev_level <- rep("High", length(grid.list))
Dev_level[which(X.PSI.raw[, 1, "Dev_lo"] == 1)] <- "Low"
Dev_level[which(X.PSI.raw[, 1, "Dev_bg"] == 1)] <- "Background"

# Decommissioned well pads 3km #
sum_tab <- X.PSI.raw[,,"WellD_3km"]
dimnames(sum_tab)[[2]] <- years

sum_grid <- cov_grid %>%
  dplyr::group_by(Development, Year) %>%
  summarise(Well_1km_mean = mean(WellA_1km),
            Well_3x3km_mean = mean(WellA_3x3km),
            Roads_mean = mean(Road_1km)) %>%
  mutate(x = ifelse(Development == "HI", Year + 0.1,
                    ifelse(Development == "LO", Year, Year - 0.1)))

# Well pads #
ggplot(sum_grid, aes(x = x, y = Well_1km_mean, color = Development)) +
  geom_point() +
  scale_x_continuous(breaks = 2010:2019, labels = 2010:2019) +
  xlab("Year")

ggplot(sum_grid, aes(x = x, y = Well_3x3km_mean, color = Development)) +
  geom_point() +
  scale_x_continuous(breaks = 2010:2019, labels = 2010:2019) +
  xlab("Year")

# Roads #
ggplot(sum_grid, aes(x = Year, y = Roads_mean, color = Development)) +
  geom_point() +
  scale_x_continuous(breaks = 2010:2019, labels = 2010:2019) +
  xlab("Year")
