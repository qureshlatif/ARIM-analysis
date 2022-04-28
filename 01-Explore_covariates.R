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
  abind::abind(Cov_grid[gridID, , c("WellA_1km", "Road_1km")], along = 3)
dat2 <- dat[,1,]
for(t in years[-1]) dat2 <- rbind(dat2, dat[,as.character(t),])
#cor(cov_point %>% select(Sage, Herb:NDVI), use = "complete") # Just exploring.
# Drop NDVI - correlated with Herb at r = 0.81
cor(dat2 %>% as.data.frame() %>% select(WellA_125m:TPI_min, WellA_1km:Road_1km),
    use = "complete") %>% # All of these can be included in analysis.
  write.csv("Cor_point.csv", row.names = T)
QSLpersonal::VIF(dat2 %>% as.data.frame() %>%
                   select(WellA_125m:vrm_125m, WellA_1km:Road_1km)) # Max VIF = 1.722444
rm(dat, dat2)

# plot(cov_point %>% filter(Year == 2010) %>% pull(Road_length_125m), # Old. Need to load raw data for this.
#      cov_point %>% filter(Year == 2010) %>% pull(Road_length_2009_125m),
#      xlab = "Tiger Roads", ylab = "ODonnel Roads")

## Tabulate summary values by development category ##
cols <- c("HI", "LO", "BG", "HILO.p", "HIBG.p", "LOBG.p")
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
  select(ID, Dev_lo:Road_1km)
dat.point <- Cov_point[,1,]
for(t in years[-1]) dat.point <- rbind(dat.point, Cov_point[,as.character(t),])
dat.point <- dat.point %>% data.frame() %>%
  mutate(ID = str_c(rep(gridID, length(years)),
                    rep(years, each = length(gridID)), sep = "_"))
dat.point <- dat.grid %>%
  select(ID, Dev_lo, Dev_hi) %>%
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
  dplyr::group_by(Dev_lo, Dev_hi) %>%
  summarise(across(PJ_area:Road_1km, sum.fn))
out[1:7, c("BG", "HI", "LO")] <- dat.grid.sum %>%
  ungroup() %>%
  select(PJ_area:Road_1km) %>%
  as.matrix() %>% t()

dat.point.sum <- dat.point %>%
  dplyr::group_by(Dev_lo, Dev_hi) %>%
  summarise(across(WellA_125m:TPI_point, sum.fn))
out[8:nrow(out), c("BG", "HI", "LO")] <- dat.point.sum %>%
  ungroup() %>%
  select(WellA_125m:TPI_point) %>%
  as.matrix() %>% t()

# Calculate p-values
dlo <- as.numeric(Cov_grid[,,"Dev_lo"])
dhi <- as.numeric(Cov_grid[,,"Dev_hi"])
for(i in 1:7) {
  x <- as.numeric(Cov_grid[,,rows[i]])[which(dhi == 1)]
  y <- as.numeric(Cov_grid[,,rows[i]])[which(dlo == 1)]
  out[i, "HILO.p"] <-
    round(t.test(x, y, alternative = "two.sided")$p.value, digits = 5)
  x <- as.numeric(Cov_grid[,,rows[i]])[which(dhi == 1)]
  y <- as.numeric(Cov_grid[,,rows[i]])[which(dlo == 0 & dhi == 0)]
  out[i, "HIBG.p"] <-
    round(t.test(x, y, alternative = "two.sided")$p.value, digits = 5)
  x <- as.numeric(Cov_grid[,,rows[i]])[which(dlo == 1)]
  y <- as.numeric(Cov_grid[,,rows[i]])[which(dlo == 0 & dhi == 0)]
  out[i, "LOBG.p"] <-
    round(t.test(x, y, alternative = "two.sided")$p.value, digits = 5)
}
dlo <- as.numeric(Cov_grid[gridID,,"Dev_lo"])
dhi <- as.numeric(Cov_grid[gridID,,"Dev_hi"])
for(i in 8:nrow(out)) {
  x <- as.numeric(Cov_point[,,rows[i]])[which(dhi == 1)]
  y <- as.numeric(Cov_point[,,rows[i]])[which(dlo == 1)]
  out[i, "HILO.p"] <-
    round(t.test(x, y, alternative = "two.sided")$p.value, digits = 5)
  x <- as.numeric(Cov_point[,,rows[i]])[which(dhi == 1)]
  y <- as.numeric(Cov_point[,,rows[i]])[which(dlo == 0 & dhi == 0)]
  out[i, "HIBG.p"] <-
    round(t.test(x, y, alternative = "two.sided")$p.value, digits = 5)
  x <- as.numeric(Cov_point[,,rows[i]])[which(dlo == 1)]
  y <- as.numeric(Cov_point[,,rows[i]])[which(dlo == 0 & dhi == 0)]
  out[i, "LOBG.p"] <-
    round(t.test(x, y, alternative = "two.sided")$p.value, digits = 5)
}

write.csv(out, "SumStats.csv", row.names = T)

