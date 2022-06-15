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

dat <- Cov_grid[,1,]
for(t in years[-1]) dat <- rbind(dat, Cov_grid[,as.character(t),])
#cor(cov_grid %>% select(PJ_area:NDVI), use = "complete") # Just for exploring
cor(dat) %>% # All of these can be included in analysis
  write.csv("Cor_grid.csv", row.names = T)
QSLpersonal::VIF(dat %>% as.data.frame() %>%
                   select(Dev_lo:NDVI, Well_3km, Road_1km)) # Max VIF = 1.411699 for Dev_lo.

# plot(cov_grid %>% filter(Year == 2010) %>% pull(Road_length_km), # Old. Need to load raw data for this.
#      cov_grid %>% filter(Year == 2010) %>% pull(Road_length_2009_km),
#      xlab = "Tiger Roads", ylab = "ODonnel Roads")

gridID <- cov_point[, "gridIndex"]
dat <- Cov_point %>%
  abind::abind(Cov_grid[gridID, , c("Well_1km", "Road_1km")], along = 3)
dat2 <- dat[,1,]
for(t in years[-1]) dat2 <- rbind(dat2, dat[,as.character(t),])
#cor(cov_point %>% select(Sage, Herb:NDVI), use = "complete") # Just exploring.
# Drop NDVI - correlated with Herb at r = 0.81
cor(dat2 %>% as.data.frame() %>% select(Well_1km, Road_125m:TPI_min, Road_1km),
    use = "complete") %>% # All of these can be included in analysis.
  write.csv("Cor_point.csv", row.names = T)
QSLpersonal::VIF(dat2 %>% as.data.frame() %>%
                   select(Road_125m:TPI_min, Road_1km)) # Max VIF = 1.722444
rm(dat, dat2)

# plot(cov_point %>% filter(Year == 2010) %>% pull(Road_length_125m), # Old. Need to load raw data for this.
#      cov_point %>% filter(Year == 2010) %>% pull(Road_length_2009_125m),
#      xlab = "Tiger Roads", ylab = "ODonnel Roads")

## Tabulate summary values by development category ##
cols <- c("HI", "LO", "BG")
rows <- c(dimnames(Cov_grid)[[3]][-c(1:2)],
          dimnames(Cov_point)[[3]])
out <- matrix("", nrow = length(rows),
              ncol = length(cols),
              dimnames = list(rows, cols))

# Elongate covariate arrays into data frames
dat.grid <- Cov_grid[,1,]
for(t in years[-1]) dat.grid <- rbind(dat.grid, Cov_grid[,as.character(t),])
dat.grid <- dat.grid %>% data.frame() %>%
  mutate(ID = str_c(rep(1:length(grid.list), length(years)),
                    rep(years, each = length(grid.list)), sep = "_")) %>%
  select(ID, Dev_bg:Road_1km)
dat.point <- Cov_point[,1,]
for(t in years[-1]) dat.point <- rbind(dat.point, Cov_point[,as.character(t),])
gridID <- as.integer(as.factor(cov_point$Grid))
dat.point <- dat.point %>% data.frame() %>%
  mutate(ID = str_c(rep(gridID, length(years)),
                    rep(years, each = length(gridID)), sep = "_"))
dat.point <- dat.grid %>%
  select(ID, Dev_lo, Dev_bg) %>%
  left_join(dat.point, by = "ID")

# Calculate summary stats
sum.fn <- function(x) str_c(mean(x, na.rm = T) %>% round(digits = 2),
                            " (",
                            sd(x, na.rm = T) %>% round(digits = 2),
                            ", ",
                            min(x, na.rm = T) %>% round(digits = 2),
                            "-",
                            max(x, na.rm = T) %>% round(digits = 2),
                            ")")
dat.grid.sum <- dat.grid %>%
  dplyr::group_by(Dev_lo, Dev_bg) %>%
  summarise(across(PJ_area:Road_1km, sum.fn))
out[names(dat.grid.sum)[-c(1:2)], c("HI", "BG", "LO")] <- dat.grid.sum %>%
  ungroup() %>%
  select(PJ_area:Road_1km) %>%
  as.matrix() %>% t()

dat.point.sum <- dat.point %>%
  dplyr::group_by(Dev_lo, Dev_bg) %>%
  summarise(across(Well_125m:TPI_point, sum.fn))
out[names(dat.point.sum)[-c(1:2)], c("HI", "BG", "LO")] <- dat.point.sum %>%
  ungroup() %>%
  select(Well_125m:TPI_point) %>%
  as.matrix() %>% t()

write.csv(out, "SumStats.csv", row.names = T)

## Survey dates by elevation ##
dat_elev <- read.csv("covariates/ARIM_elevations.csv", header = T, stringsAsFactors = F)
tab_sum <- tab.datetime %>%
  mutate(TrnsctN = str_sub(Point_year, 1, -9)) %>%
  left_join(dat_elev, by = "TrnsctN") %>%
  mutate(Elev_class = ifelse(Elev_m < 2000, "lt2000",
                             ifelse(Elev_m >= 2000 & Elev_m < 2300, "2000-2300",
                                    ifelse(Elev_m >= 2300 & Elev_m < 2600, "2300-2600", "gt2600")))) %>%
  dplyr::group_by(Elev_class) %>%
  summarise(minDOY = round(quantile(DOY, prob = 0.01, type = 8)),
            maxDOY = round(quantile(DOY, prob = 0.99, type = 8)),
            meanDOY = round(mean(DOY, na.rm = T)))
tab_ref <- tab.datetime %>%
  mutate(Date = str_sub(Date, 6, -1)) %>%
  arrange(DOY) %>%
  select(DOY, Date) %>%
  distinct() %>%
  dplyr::group_by(DOY) %>%
  summarise(first = first(Date), last = last(Date))
tab_sum <- tab_sum %>%
  left_join(
    tab_ref %>%
      select(DOY, first),
    by = c("minDOY" = "DOY")
  ) %>%
  left_join(
    tab_ref %>%
      select(DOY, last),
    by = c("maxDOY" = "DOY")
  ) %>%
  left_join(
    tab_ref %>%
      rename(mean = first) %>%
      select(DOY, mean),
    by = c("meanDOY" = "DOY")
  )
View(tab_sum)
