library(stringr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")

cov_grid <- read.csv("covariates/Grid_Covariates_ARIM.csv", stringsAsFactors = F, header = T) %>%
  rename(TransectNum = TrnsctN) %>%
  tibble::as_tibble() %>%
  mutate(PJ_area = str_sub(PJ_area, 1, -2) %>% as.numeric,
         Road_length_2009_km = suppressWarnings(Road_length_2009_km %>% as.numeric()))

cov_point <- read.csv("covariates/Point_Covariates_ARIM.csv", stringsAsFactors = F, header = T) %>%
  tibble::as_tibble() %>%
  rename(TransectNum = TrnsctN) %>%
  mutate(Sage = suppressWarnings(as.numeric(Sage)), # 123 values are missing for each of these.
         Litter = suppressWarnings(as.numeric(Litter)),
         Herb = suppressWarnings(as.numeric(Herb)),
         AHerb = suppressWarnings(as.numeric(AHerb)))

#   # Attach Mutter et al. road and well pad densities #
# cov_grid_compMutt <- cov_grid %>%
#   left_join(
#     read.csv("covariates/ARIM_2010-2012_MutterEtAl_roads.csv", header = T) %>%
#       rename(RD_Mutt = Road_density,
#              Well_Mutt = Well_density),
#     by = c("TransectNum" = "Transect", "Year" = "Year")
#   )
# cor(cov_grid_compMutt %>% select(Well_count_1km, Well_count_3x3km, Road_length_km, Road_length_2009_km:Well_Mutt), use = "complete")

cor(cov_grid %>% select(PJ_area, Well_count_1km:Road_length_2009_km), use = "complete")
cor(cov_grid %>% filter(Year == 2010) %>%
      select(PJ_area, Well_count_1km:Road_length_2009_km), use = "complete")
plot(cov_grid %>% filter(Year == 2010) %>% pull(Road_length_km),
     cov_grid %>% filter(Year == 2010) %>% pull(Road_length_2009_km),
     xlab = "Tiger Roads", ylab = "ODonnel Roads")

cor(cov_point %>% select(Sage:min_elevation_difference), use = "complete")
cor(cov_point %>% filter(Year == 2010) %>%
      select(Sage:min_elevation_difference), use = "complete")
plot(cov_point %>% filter(Year == 2010) %>% pull(Road_length_125m),
     cov_point %>% filter(Year == 2010) %>% pull(Road_length_2009_125m),
     xlab = "Tiger Roads", ylab = "ODonnel Roads")

plot(cov_point %>% filter(Year == 2010) %>% pull(elevation_difference),
     cov_point %>% filter(Year == 2010) %>% pull(min_elevation_difference),
     xlab = "Point elevation difference", ylab = "Min nbrhd elev difference")

## Tabulate summary values by development category ##
cov_grid <- cov_grid %>%
  mutate(Development = ifelse(str_detect(TransectNum, "ARIM-HI"), "HI",
                              ifelse(str_detect(TransectNum, "ARIM-LO"), "LO", "BG")))
cov_point <- cov_point %>%
  mutate(Development = ifelse(str_detect(TransectNum, "ARIM-HI"), "HI",
                              ifelse(str_detect(TransectNum, "ARIM-LO"), "LO", "BG")))

cols <- c("HI", "LO", "BG")
rows <- c(names(cov_grid)[c(7, 9:11, 13:14)],
          names(cov_point)[c(9:18)])
out <- matrix("", nrow = length(rows),
              ncol = length(cols),
              dimnames = list(rows, cols))
out[1:6, "HI"] <- (cov_grid %>%
                     filter(Development == "HI") %>%
                     select(matches(rows[1:6])) %>%
                     QSLpersonal::SumStats_df())[,1]
out[1:6, "LO"] <- (cov_grid %>%
                     filter(Development == "LO") %>%
                     select(matches(rows[1:6])) %>%
                     QSLpersonal::SumStats_df())[,1]
out[1:6, "BG"] <- (cov_grid %>%
                     filter(Development == "BG") %>%
                     select(matches(rows[1:6])) %>%
                     QSLpersonal::SumStats_df())[,1]

out[7:16, "HI"] <- (cov_point %>%
                     filter(Development == "HI") %>%
                     select(matches(rows[7:16])) %>%
                     QSLpersonal::SumStats_df())[,1]
out[7:16, "LO"] <- (cov_point %>%
                     filter(Development == "LO") %>%
                     select(matches(rows[7:16])) %>%
                     QSLpersonal::SumStats_df())[,1]
out[7:16, "BG"] <- (cov_point %>%
                     filter(Development == "BG") %>%
                     select(matches(rows[7:16])) %>%
                     QSLpersonal::SumStats_df())[,1]

## Plot well pad counts and road density over time ##
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

sum_grid <- cov_grid %>%
  dplyr::group_by(Development, Year) %>%
  summarise(Well_count_mn = mean(Well_count_1km),
            #Well_count_lo = mean(Well_count_1km) -
            #  1.96 * (sd(Well_count_1km) / sqrt(sum(!is.na(Well_count_1km)))),
            #Well_count_hi = mean(Well_count_1km) +
            #  1.96 * (sd(Well_count_1km) / sqrt(sum(!is.na(Well_count_1km)))),
            #nWell_count = sum(!is.na(Well_count_1km)),
            Roads_mn = mean(Road_length_km)#,
            #Roads_lo = quantile(Road_length_km, prob = 0.05, type = 8),
            #Roads_hi = quantile(Road_length_km, prob = 0.95, type = 8)
            ) %>%
  mutate(x = ifelse(Development == "HI", Year + 0.1,
                    ifelse(Development == "LO", Year, Year - 0.1)))

# Well pads #
ggplot(sum_grid, aes(x = x, y = Well_count_mn, color = Development)) +
  geom_point() +
  #geom_errorbar(aes(ymin = Well_count_lo,
  #                  ymax = Well_count_hi)) +
  scale_x_continuous(breaks = 2010:2019, labels = 2010:2019) +
  xlab("Year")

ggplot(sum_grid, aes(x = Year, y = Roads_mn, color = Development)) +
  geom_point() +
  scale_x_continuous(breaks = 2010:2019, labels = 2010:2019) +
  xlab("Year")
