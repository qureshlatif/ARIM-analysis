library(stringr)
library(tidyr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#_____ Script setup _____#
scripts.loc <- "ARIM-analysis/"
mod.nam <- "mod_path"
development <- reduce.data.aug <- FALSE
mod <- R.utils::loadObject(mod.nam)
source(str_c(scripts.loc, "Data_processing.R"))
#________________________#

### Save parameter summaries ###
#write.csv(mod$summary, "Mech_pathway_parameters.csv")

## Plot mean mech values over time by development stratum ##
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

# Compile index #
Dev_level <- rep("High", n.grdyr)
Dev_level[which(X.PSI.raw[, "Dev_lo"] == 1)] <- "Low"
Dev_level[which(X.PSI.raw[, "Dev_bg"] == 1)] <- "Background"

Year <- as.character(yearID.grdyr + 2009)

# Well pads #
sum_tab <- data.frame(X.PSI.raw[, c("Well_1km", "Well_3km", "Road_1km")]) %>%
  mutate(Development = Dev_level,
         Year = Year) %>%
  dplyr::group_by(Development, Year) %>%
  summarise(Well_1km = mean(Well_1km),
            Well_3km = mean(Well_3km),
            Road_1km = mean(Road_1km)) %>%
  mutate(Year = as.integer(Year))

# Well pads #
ggplot(sum_tab, aes(x = Year, y = Well_1km, color = Development)) +
  geom_point() +
  scale_x_continuous(breaks = 2010:2019, labels = 2010:2019) +
  xlab("Year")

ggplot(sum_tab, aes(x = Year, y = Well_3km, color = Development)) +
  geom_point() +
  scale_x_continuous(breaks = 2010:2019, labels = 2010:2019) +
  xlab("Year")

# Roads #
ggplot(sum_tab, aes(x = Year, y = Road_1km, color = Development)) +
  geom_point() +
  scale_x_continuous(breaks = 2010:2019, labels = 2010:2019) +
  xlab("Year")
