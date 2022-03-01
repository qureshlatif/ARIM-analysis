library(stringr)
library(BCRDataAPI)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")

################## Inputs ####################
# Get latest AOU checklist with tax names and order #
aou.checklist <- read.csv("C:/Users/Quresh.Latif/files/data/NACC_list_bird_species_downloaded_20220204.csv",
                          header = T, stringsAsFactors = F) %>% as_tibble() %>%
  mutate(tax_ord = row_number())
# Recommended citation: Chesser, R. T., S. M. Billerman, K. J. Burns,
  #C. Cicero, J. L. Dunn, B. E. Hernández-Baños, A. W. Kratter, I. J. Lovette,
  #N. A. Mason, P. C. Rasmussen, J. V. Remsen, Jr., D. F. Stotz, and K. Winker. 2021.
  #Check-list of North American Birds (online). American Ornithological Society.
  #http://checklist.aou.org/taxa

spp.exclude <- aou.checklist %>%
  filter(!order %in% c("Apodiformes", "Caprimulgiformes", "Charadriiformes",
                       "Caprimulgiformes", "Columbiformes", "Coraciiformes",
                       "Cuculiformes", "Nyctibiiformes", "Passeriformes",
                       "Piciformes")) %>%
  pull(common_name)
spp.exclude <- c(spp.exclude,
                 "Squirrel, Red", "Squirrel, Abert's", "Pika",
                 "American Avocet", "Forster's Tern", "Franklin's Gull",
                 "Long-billed Curlew", "Marbled Godwit", "Ring-billed Gull",
                 "Willet", "Black-necked Stilt", "Wilson's Phalarope",
                 "California Gull", "Caspian Tern")
strata <- c("WY-ARIM-HI", "WY-ARIM-LO", "WY-BCR10-BU", "WY-BCR10-CA", "WY-BCR10-CO", "WY-BCR10-KE",
            "WY-BCR10-LA", "WY-BCR10-PI", "WY-BCR10-RA", "WY-BCR10-RO", "WY-BCR10-WO")
years <- 2010:2019
SampDesign <- c("IMBCR", "GRTS")
##############################################

#### Transect list ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('TransectNum|str'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ","))))
BCRDataAPI::group_by(c("TransectNum"))
grab <- BCRDataAPI::get_data()
write.csv(grab, "Transects.csv", row.names = F)

#### Get primary habitat types for sample transects ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'primaryHabitat',
                          'HabitatCommonName'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ","))))
BCRDataAPI::group_by(c('TransectNum', 'Point', 'Year', 'primaryHabitat',
                       'HabitatCommonName'))
grab <- BCRDataAPI::get_data()
#grab %>% dplyr::group_by(primaryHabitat, HabitatCommonName) %>%
#  summarise(n = n()) %>% arrange(desc(n)) %>%
#  View()
PH.list <- grab %>% dplyr::group_by(primaryHabitat, HabitatCommonName) %>%
  summarise(n = n()) %>% arrange(desc(n)) %>%
  pull(primaryHabitat)
PH.list <- PH.list[1:7]

#### Compile species list ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('BirdCode|str',
                          'Species|str')
)
BCRDataAPI::filter_on('SelectionMethod in IMBCR,GRTS')
BCRDataAPI::filter_on('Migrant = 0')
BCRDataAPI::filter_on('BirdCode <> NOBI')
BCRDataAPI::filter_on('BCR = 10')
BCRDataAPI::filter_on(str_c('primaryHabitat in ', str_c(PH.list, collapse = ",")))
BCRDataAPI::group_by(c('BirdCode', 'Species'))
grab <- BCRDataAPI::get_data()

spp.out <- grab %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(!Species %in% spp.exclude)

# Collapsing sub-species and renamed species #
ss <- BCRDataAPI::subspecies()
spp.out <- spp.out %>%
  mutate(BirdCode = ss[BirdCode] %>% as.character)%>%
  dplyr::group_by(BirdCode) %>%
  mutate(min_length = min(nchar(Species))) %>%
  mutate(Species = str_sub(Species, 1, min_length)) %>%
  select(BirdCode, Species) %>%
  # Additional adjustment:
  mutate(Species = ifelse(BirdCode == "WEFL",
                              "Cordilleran Flycatcher", Species)) %>%
  mutate(BirdCode = ifelse(BirdCode == "WEFL", "COFL", BirdCode)) %>%
  #_____________________#
  ungroup %>%
  unique

#sum(!spp.out$Species %in% aou.checklist$common_name) # check - should be zero
spp.out <- spp.out %>%
  rename(common_name = Species) %>%
  left_join(aou.checklist, by = "common_name") %>%
  select(tax_ord, BirdCode, common_name, species) %>%
  arrange(tax_ord)

# Join guilds classified for report #
spp_guilds <- read.csv("Spp_guilds_prior.csv", header = T, stringsAsFactors = F) %>%
  mutate(Common_Name = ifelse(Common_Name == "Gray Jay", "Canada Jay", Common_Name)) %>%
  mutate(Common_Name = ifelse(Common_Name == "Western Scrub-Jay", "Woodhouse's Scrub-Jay", Common_Name)) %>%
  mutate(Common_Name = ifelse(Common_Name == "Sage Sparrow", "Sagebrush Sparrow", Common_Name)) %>%
  mutate(Common_Name = ifelse(Common_Name == "McCown's Longspur", "Thick-billed Longspur", Common_Name))
#spp.out %>% filter(!common_name %in% spp_guilds$Common_Name) %>% pull(common_name) # Review this....
#spp_guilds %>% filter(!Common_Name %in% spp.out$common_name) %>% pull(Common_Name) %>% sort() # ...and this.
spp.out <- spp.out %>%
  left_join(
    spp_guilds %>%
      select(Common_Name, Guilds) %>%
      rename(common_name = Common_Name,
             Guild = Guilds),
    by = "common_name"
  )
rm(spp_guilds)
spp.out <- spp.out %>%
  mutate(Guild = ifelse(common_name == "Cassin's Sparrow", "Range", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Orchard Oriole", "Generalist", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Yellow Grosbeak", "Generalist", Guild))

# Remove additional implausible members of the metacommunity (based on review of BNA range maps and habitat accounts) #
#spp.out <- spp.out %>%
#  filter(!BirdCode %in% c("RUHU", "PAWR", "OLWA", "AMPI", "WWCR", "SABS"))

spp.excluded <- grab %>%
  select(BirdCode, Species) %>%
  unique %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(Species %in% spp.exclude) %>%
  select(BirdCode, Species) %>%
  unique %>%
  rename(common_name = Species) %>%
  left_join(aou.checklist, by = "common_name") %>%
  select(tax_ord, BirdCode, common_name, species) %>%
  arrange(tax_ord)

#### Detection data ####
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'PointLatitude|num',
                          'PointLongitude|num',
                          'easting|int',
                          'northing|int',
                          'zone|int',
                          'Year|int',
                          'Date|str',
                          'PointVisitStartTime|str',
                          'Stratum|str',
                          'radialDistance|int',
                          'CL_count|int',
                          'BirdCode|str',
                          'Species|str',
                          'How|str',
                          'Sex|str',
                          'TimePeriod|int'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ",")),
                        str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                        str_c('Year in ', str_c(years, collapse = ",")),
                        str_c('BirdCode in ', str_c(spp.out$BirdCode, collapse = ",")),
                        'ninetynine = 0',
                        'eightyeight = 0',
                        'How <> F',
                        'Sex <> J',
                        'Migrant = 0',
                        'TimePeriod > -1',
                        'radialDistance < 125'))
grab <- BCRDataAPI::get_data(interpolate_effort = T) %>%
  mutate(BirdCode = ss[BirdCode] %>% as.character)

point.coords <- grab %>%
  select(TransectNum, Point, easting, northing, zone) %>%
  distinct()
point.list <- unique(str_c(point.coords$TransectNum,
                           str_pad(point.coords$Point, width = 2, side = "left", pad = "0"),
                           sep = "-")) %>%
  sort
grid.list <- unique(point.coords$TransectNum) %>% sort

## Point X years surveyed ##
pointXyears.list <- unique(str_c(grab$TransectNum,
                                 str_pad(grab$Point, width = 2,
                                         side = "left", pad = "0"),
                                 grab$Year, sep = "-")) %>% sort

## Add number of detections and count summaries to spp.out by stratum ##
smry <- grab %>% select(BirdCode, TransectNum, Point, Year) %>%
  unique %>% dplyr::group_by(BirdCode) %>% count() %>%
  rename(Detections = n)
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

smry <- grab %>% select(BirdCode, CL_count) %>%
  dplyr::group_by(BirdCode) %>%
  summarise(sumCount = sum(CL_count))
spp.out <- spp.out %>% left_join(smry, by = "BirdCode")

spp.out <- spp.out %>% # replace NAs with zeros
  mutate_at(vars(Detections, sumCount), (function(x) replace(x, is.na(x), 0)))

maxDetPossible <- length(pointXyears.list) # max possible by stratum
names(spp.out)[which(names(spp.out) == "Detections")] <-
  str_c("Detections (max = ", maxDetPossible, ")")

write.csv(spp.out, "Spp_list.csv", row.names = F)
rm(smry)

## Add number of detections and count summaries to excluded species ##
smry <- grab %>% select(BirdCode, TransectNum, Point, Year) %>% unique %>%
  dplyr::group_by(BirdCode) %>% count() %>%
  rename(Detections = n)
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

smry <- grab %>% select(BirdCode, TransectNum, Point, CL_count) %>%
  dplyr::group_by(BirdCode) %>%
  summarise(sumCount = sum(CL_count))
spp.excluded <- spp.excluded %>% left_join(smry, by = "BirdCode")

spp.excluded <- spp.excluded %>% # replace NAs with zeros
  mutate_at(vars(Detections, sumCount),
            (function(x) replace(x, is.na(x), 0)))

write.csv(spp.excluded, "Spp_excluded.csv", row.names = F)
rm(smry)

bird_data <- grab %>%  # Store bird survey data for later use.
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year))

## Covariates ##
# # Apply slope filter on units represented in data from GIS shop
# PY.GIS <- read.csv("covariates/Point_Covariates_ARIM.csv", stringsAsFactors = F, header = T) %>%
#   mutate(PY = str_c(TrnsctN, str_pad(Point, width = 2, side = "left", pad = "0"), Year, sep = "-")) %>%
#   pull(PY)
# pointXyears.list <- pointXyears.list[which(pointXyears.list %in% PY.GIS)]
# rm(PY.GIS)
# point.list <- str_sub(pointXyears.list, 1, -6) %>% unique() %>% sort
# grid.list <- str_sub(pointXyears.list, 1, -9) %>% unique() %>% sort

# GIS #
cov_grid <- data.frame(Grid = pointXyears.list %>% str_sub(1, -9),
                             Year = pointXyears.list %>%
                         str_sub(-4, -1) %>% as.integer(),
                             stringsAsFactors = F) %>%
  mutate(gridIndex = Grid %>% as.factor %>% as.integer,
         YearInd = Year %>% as.factor %>% as.integer,
         Development = ifelse(str_detect(Grid, "ARIM-HI"), "HI",
                             ifelse(str_detect(Grid, "ARIM-LO"), "LO", "BG"))) %>%
  distinct() %>%
  left_join(read.csv("covariates/Grid_Covariates_ARIM_2.csv",
                     stringsAsFactors = F, header = T) %>%
              tibble::as_tibble() %>%
              mutate(PJ_area = str_sub(PJ_area, 1, -2) %>% as.numeric) %>%
              rename(TransectNum = TrnsctN,
                     Road_1km_2009 = Road_length_2009_Odonnell_km,
                     Road_1km_2019 = Road_length_2019_km,
                     WellA_1km = Well_count_1km,
                     WellA_3x3km = Well_count_3x3km,
                     WellD_1km = PA_Well_count_1km,
                     WellD_3x3km = PA_Well_count_3x3km) %>%
              mutate(Road_1km = ifelse(Year == 2010,
                                       Road_1km_2009,
                                       Road_1km_2019)) %>%
              select(TransectNum, Year, PJ_area,
                     WellA_1km, WellA_3x3km, WellD_1km, WellD_3x3km,
                     Road_1km, Road_1km_2009, Road_1km_2019,
                     NDVI),
    by = c("Grid" = "TransectNum", "Year" = "Year"))

  # Check number of years sampled for each grid cell #
# out <- matrix(0, nrow = sum(str_sub(grid.list, 4, 7) == "ARIM"),
#               ncol = length(years),
#               dimnames = list(grid.list[which(str_sub(grid.list, 4, 7) == "ARIM")],
#                               as.character(years)))
# ind <- which(str_sub(cov_grid$Grid, 4, 7) == "ARIM")
# for(i in ind) out[cov_grid$Grid[i], as.character(cov_grid$Year[i])] <- 1

cov_point <- data.frame(Point_year = pointXyears.list,
                        Grid = pointXyears.list %>% str_sub(1, -9),
                        Point = pointXyears.list %>% str_sub(1, -6),
                        Year = pointXyears.list %>%
                          str_sub(-4, -1) %>% as.integer(),
                        stringsAsFactors = F) %>%
  mutate(gridIndex = Grid %>% as.factor %>% as.integer,
         pointIndex = Point %>% as.factor %>% as.integer,
         YearInd = Year %>% as.factor %>% as.integer,
         Development = ifelse(str_detect(Point, "ARIM-HI"), "HI",
                              ifelse(str_detect(Point, "ARIM-LO"), "LO", "BG"))) %>%
  distinct() %>%
  dplyr::left_join(
    read.csv("covariates/Point_Covariates_ARIM_2.csv", stringsAsFactors = F, header = T) %>%
      tibble::as_tibble() %>%
      rename(TransectNum = TrnsctN) %>%
      mutate(Point = str_c(TransectNum,
                           str_pad(Point, width = 2, side = "left", pad = "0"),
                           sep = "-"),
             Sage = suppressWarnings(as.numeric(Sage)), # 123 values are missing for each of these.
             Litter = suppressWarnings(as.numeric(Litter)),
             Herb = suppressWarnings(as.numeric(Herb)),
             AHerb = suppressWarnings(as.numeric(AHerb)),
             vrm_125m = suppressWarnings(as.numeric(vrm_125m))) %>%
      rename(Road_125m_2009 = Roads_2009_km_Odonnell,
             Road_125m_2019 = Road_length_2019_km,
             TPI_min = Min_elev_diff,
             TPI_mean = Mean_elevation_difference,
             TPI_point = Point_elev_diff,
             WellA_125m = Well_count_125m
             ) %>%
      mutate(Road_125m = ifelse(Year == 2010, Road_125m_2009, Road_125m_2019),
             WellA_125m = ifelse(WellA_125m > 0, 1, 0)) %>% # Convert to presence / absence of active well pad (Only 17 points with decom well pad within 125 m).
      select(Point, Year, Sage:AHerb, WellA_125m, Road_125m,
             Road_125m_2009:vrm_125m, TPI_mean, TPI_point, TPI_min),
    by = c("Point", "Year")
  )

cov_point <- cov_point %>%
  left_join(
    cov_grid %>%
      select(Grid:Year, WellA_1km, WellD_1km, Road_1km, NDVI),
    by = c("Grid", "Year")
  )

# IMBCR woody vegetation cover #
BCRDataAPI::reset_api()
BCRDataAPI::set_api_server('analysis.api.bcr.eco')
BCRDataAPI::add_columns(c('TransectNum|str',
                          'Point|int',
                          'Year|int',
                          'o_canopy_percent|int',
                          'shrub_cover|int'
))

BCRDataAPI::filter_on(c(str_c('Stratum in ', str_c(strata, collapse = ",")),
                        str_c('SelectionMethod in ', str_c(SampDesign, collapse = ",")),
                        'TransectExcludeAnalysis = FALSE'
))
BCRDataAPI::group_by(c('TransectNum',
                       'Point',
                       'Year',
                       'o_canopy_percent',
                       'shrub_cover'))
grab <- BCRDataAPI::get_data() %>%
  mutate(Point_year = str_c(TransectNum, "-",
                            str_pad(Point, width = 2, pad = "0",
                                    side = "left"), "-", Year)) %>%
  rename(CanCov = o_canopy_percent,
         ShrubCov = shrub_cover) %>%
  mutate(CanCov = ifelse(CanCov == -1, NA, CanCov),
         ShrubCov = ifelse(ShrubCov == -1, NA, ShrubCov))

cov_point <- cov_point %>%
  left_join(
    grab %>% select(Point_year, CanCov, ShrubCov),
    by = "Point_year"
  )

## Fill 123 missing ground cover values ##
library(randomForest)
cov_point <- cov_point %>%
  left_join(
    cov_point %>%
      group_by(Grid, Year) %>%
      summarise(Sage_grid = mean(Sage, na.rm = T),
                Litter_grid = mean(Litter, na.rm = T),
                Herb_grid = mean(Herb, na.rm = T),
                AHerb_grid = mean(AHerb, na.rm = T)),
    by = c("Grid", "Year")
  )

v.fill <- c("Sage", "Herb", "AHerb")
v.inform <- c("vrm_125m", "TPI_mean", "TPI_point", "TPI_min", "Sage_grid", "Litter_grid", "Herb_grid", "AHerb_grid")
for(v in 1:length(v.fill)) {
  ind.missing <- which(is.na(cov_point[,v.fill[v]]))
  ind.known <- which(!is.na(cov_point[,v.fill[v]]))
  rf <- randomForest(as.formula(str_c(v.fill[v], "~", str_c(v.inform, collapse = "+"))), data = (cov_point %>% slice(ind.known)))
  cov_point[ind.missing, v.fill[v]] <- predict(rf, newdata = (cov_point %>% slice(ind.missing)))
}
cov_point <- cov_point %>% select(Point_year:ShrubCov)
rm(v.fill, v.inform, ind.missing, ind.known, rf, v)

## Trim dates, compile day of year & start time in minutes ##
library(lubridate)
tab.datetime <- bird_data %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year)) %>%
  filter(Point_year %in% pointXyears.list) %>%
  arrange(Point_year) %>%
  select(Point_year, PointLatitude, PointLongitude, Date, PointVisitStartTime) %>%
  distinct() %>%
  mutate(Date = str_sub(Date, 6, -14) %>% dmy) %>%
  mutate(DOY = yday(Date)) %>%
  mutate(PointVisitStartTime = PointVisitStartTime %>%
           replace(which(PointVisitStartTime == "0"), NA)) %>%
  mutate(HR = PointVisitStartTime %>% str_sub(1, -3) %>% as.integer()) %>%
  mutate(HR = HR %>% str_pad(width = 2, side = "left", pad = "0")) %>%
  mutate(MIN = PointVisitStartTime %>% str_sub(-2, -1) %>% str_pad(width = 2, side = "left", pad = "0")) %>%
  mutate(Time = str_c(HR, MIN, "00", sep = ":")) %>%
  mutate(dateTime = str_c(Date, Time, sep = " ")) %>%
  mutate(Time_ssr = QSLpersonal::tssr(PointLatitude, PointLongitude, dateTime)) %>%
  select(Point_year, PointLatitude, PointLongitude, DOY, Time_ssr)

## Compile multidimensional detection data array ##
spp.list <- spp.out$BirdCode

bird_data <- bird_data %>%
  mutate(Point_year = str_c(TransectNum, "-", str_pad(Point, width = 2, pad = "0", side = "left"), "-", Year))
Y.mat <- matrix(NA, nrow = length(pointXyears.list), ncol = length(spp.list),
               dimnames = list(pointXyears.list, spp.list))
TR.mat <- matrix(6, nrow = length(pointXyears.list), ncol = length(spp.list),
               dimnames = list(pointXyears.list, spp.list))
for(sp in 1:length(spp.list)) {
  obs <- bird_data %>% filter(BirdCode == spp.list[sp] & Point_year %in% pointXyears.list)
  if(nrow(obs) > 0) {
    Y.mat[, sp] <- (pointXyears.list %in% obs$Point_year) %>% as.integer
    tvec <- tapply(obs$TimePeriod, obs$Point_year, min)
    tvec <- tvec[order(names(tvec))]
    TR.mat[which(pointXyears.list %in% obs$Point_year), sp] <- tvec
  } else {
    Y.mat[, sp] <- 0
  }
}

## Compile covariate arrays ##
cov.names <- c("gridIndex", "YearInd", "pointInd", "DayOfYear",
               "Time_ssr", "Develop_low", "Develop_high",
               names(cov_point)[-c(1:8)])
Cov_point <- matrix(NA, nrow = length(pointXyears.list),
                    ncol = length(cov.names),
                    dimnames = list(pointXyears.list, cov.names))
Cov_point[, "gridIndex"] <- pointXyears.list %>% str_sub(1, -9) %>% as.factor %>% as.integer
Cov_point[, "YearInd"] <- pointXyears.list %>% str_sub(-4, -1) %>% as.factor %>% as.integer
Cov_point[, "pointInd"] <- pointXyears.list %>% str_sub(1, -6) %>% as.factor %>% as.integer
Cov_point[, "DayOfYear"] <- tab.datetime %>% arrange(Point_year) %>% pull(DOY)
Cov_point[, "Time_ssr"] <- tab.datetime %>% arrange(Point_year) %>% pull(Time_ssr)
Cov_point[, "Develop_low"] <- cov_point %>% arrange(Point_year) %>% pull(Development) %>%
  (function(x) (x == "LO") * 1)
Cov_point[, "Develop_high"] <- cov_point %>% arrange(Point_year) %>% pull(Development) %>%
  (function(x) (x == "HI") * 1)
Cov_point[, cov.names[-c(1:7)]] <- (cov_point %>%
  arrange(Point_year) %>%
  select(matches(cov.names[-c(1:7)])) %>%
  data.matrix())[,cov.names[-c(1:7)]]

cov.names <- c("Develop_low", "Develop_high", names(cov_grid)[-c(1:5)])
Cov_grid <- array(NA, dim = c(length(grid.list), length(years), length(cov.names)),
              dimnames = list(grid.list, years, cov.names))
Cov_grid[, , "Develop_low"] <- str_detect(grid.list, "ARIM-LO") * 1
Cov_grid[, , "Develop_high"] <- str_detect(grid.list, "ARIM-HI") * 1
for(i in 1:nrow(cov_grid)) Cov_grid[cov_grid$gridIndex[i], cov_grid$YearInd[i], -c(1:2)] <-
  (cov_grid %>% slice(i) %>%
     select(matches(cov.names[-c(1:2)])) %>%
     data.matrix())[,cov.names[-c(1:2)]] %>% as.numeric()

rm(obs, maxDetPossible, sp, ss, tvec, grab, cov.names, i)
save.image("Data_compiled.RData")
