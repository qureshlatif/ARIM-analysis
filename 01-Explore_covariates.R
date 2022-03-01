library(stringr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

#   # Attach Mutter et al. road and well pad densities #
# cov_grid_compMutt <- cov_grid %>%
#   left_join(
#     read.csv("covariates/ARIM_2010-2012_MutterEtAl_roads.csv", header = T) %>%
#       rename(RD_Mutt = Road_density,
#              Well_Mutt = Well_density),
#     by = c("Grid" = "Transect", "Year" = "Year")
#   )
# cor(cov_grid_compMutt %>% select(WellA_1km:Well_Mutt), use = "complete")
# Correlation between WellA_3x3km and Mutter et al.'s well pad density is r = 0.6.
# Correlation between Road_1km and Mutter et al.'s road density is r = 0.4.

cor(cov_grid %>% select(PJ_area:NDVI), use = "complete") # Just for exploring
cor(cov_grid %>% select(PJ_area, WellA_3x3km, WellD_3x3km, Road_1km, NDVI),
    use = "complete") %>% # All of these can be included in analysis
  write.csv("Cor_grid.csv", row.names = T)
QSLpersonal::VIF(cov_grid %>% select(PJ_area, WellA_3x3km, WellD_3x3km, Road_1km, NDVI)) # Max VIF = 1.340998 for road density.

# plot(cov_grid %>% filter(Year == 2010) %>% pull(Road_length_km), # Old. Need to load raw data for this.
#      cov_grid %>% filter(Year == 2010) %>% pull(Road_length_2009_km),
#      xlab = "Tiger Roads", ylab = "ODonnel Roads")

cor(cov_point %>% select(Sage, Herb:NDVI), use = "complete") # Just exploring.
# Drop NDVI - correlated with Herb at r = 0.81
cor(cov_point %>% select(Sage, Herb:Road_125m, vrm_125m, TPI_min:Road_1km),
    use = "complete") %>% # All of these can be included in analysis.
  write.csv("Cor_point.csv", row.names = T)
QSLpersonal::VIF(cov_point %>% select(Sage, Herb:Road_125m, vrm_125m, WellA_1km:Road_1km)) # Max VIF = 1.722444

# plot(cov_point %>% filter(Year == 2010) %>% pull(Road_length_125m), # Old. Need to load raw data for this.
#      cov_point %>% filter(Year == 2010) %>% pull(Road_length_2009_125m),
#      xlab = "Tiger Roads", ylab = "ODonnel Roads")

## Tabulate summary values by development category ##
cols <- c("HI", "LO", "BG")
rows <- c(names(cov_grid)[c(6:11, 14)],
          names(cov_point)[c(9:14, 17:20, 25, 26)])
out <- matrix("", nrow = length(rows),
              ncol = length(cols),
              dimnames = list(rows, cols))
out[1:7, "HI"] <- (cov_grid %>%
                     filter(Development == "HI") %>%
                     QSLpersonal::SumStats_df(vars = rows[1:7]))[,1]
out[1:7, "LO"] <- (cov_grid %>%
                     filter(Development == "LO") %>%
                     QSLpersonal::SumStats_df(vars = rows[1:7]))[,1]
out[1:7, "BG"] <- (cov_grid %>%
                     filter(Development == "BG") %>%
                     QSLpersonal::SumStats_df(vars = rows[1:7]))[,1]

out[8:19, "HI"] <- (cov_point %>%
                     filter(Development == "HI") %>%
                      QSLpersonal::SumStats_df(vars = rows[8:19]))[,1]
out[8:19, "LO"] <- (cov_point %>%
                     filter(Development == "LO") %>%
                      QSLpersonal::SumStats_df(vars = rows[8:19]))[,1]
out[8:19, "BG"] <- (cov_point %>%
                     filter(Development == "BG") %>%
                      QSLpersonal::SumStats_df(vars = rows[8:19]))[,1]

write.csv(out, "SumStats.csv", row.names = T)

## Plot well pad counts and road density over time ##
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

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
