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

#############
# Bird data #
#############

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
  mutate(Guild = ifelse(common_name == "Cassin's Sparrow", "Grassland", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Orchard Oriole",   "Generalist", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Vaux's Swift",     "Montane", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Cassin's Vireo",   "Woodland", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Pacific Wren",     "Generalist", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Purple Martin",    "Generalist", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Bushtit",          "Generalist", Guild)) %>%
  mutate(Guild = ifelse(common_name == "Western Bluebird", "Generalist", Guild)) %>%
  filter(!common_name %in% c("Gray Vireo", "Varied Thrush",
                             "Sprague's Pipit", "Baird's Sparrow",
                             "Northern Waterthrush", "Nashville Warbler",
                             "Yellow Grosbeak", "Greater Yellowlegs",
                             "Pileated Woodpecker", "Chestnut-backed Chickadee",
                             "Boreal Chickadee")) # Excluded: breeding range nowhere near study area.
  
spp.excluded <- grab %>%
  select(BirdCode, Species) %>%
  unique %>%
  filter(!str_sub(BirdCode, 1, 2) == "UN") %>%
  filter(Species %in% spp.exclude |
           Species %in% c("Gray Vireo", "Varied Thrush",
                          "Sprague's Pipit", "Baird's Sparrow",
                          "Northern Waterthrush", "Nashville Warbler",
                          "Yellow Grosbeak", "Greater Yellowlegs",
                          "Pileated Woodpecker", "Chestnut-backed Chickadee",
                          "Boreal Chickadee")) %>% # Excluded: breeding range nowhere near study area.
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

##############
# Covariates #
##############

### GIS ###
# Grid cells #
cov_grid <- data.frame(Grid = pointXyears.list %>% str_sub(1, -9),
                       stringsAsFactors = F) %>%
  mutate(gridIndex = Grid %>% as.factor %>% as.integer,
         Development = ifelse(str_detect(Grid, "ARIM-HI"), "HI",
                              ifelse(str_detect(Grid, "ARIM-LO"), "LO", "BG"))) %>%
  distinct() %>%
  arrange(gridIndex) %>%
  left_join(
    read.csv("covariates/ARIM_grid_covariates_all_years.csv",
             stringsAsFactors = F, header = T) %>%
      tibble::as_tibble() %>%
      mutate(across(PJ_area_2010:PJ_area_2019, function(x) str_sub(x, 1, -2) %>% as.numeric)) %>%
      rename(Grid = TrnsctN,
             Road_1km_2009 = Road_length_Odonnell_km,
             Road_1km_2019 = Road_length_2019_km) %>%
      select(Grid, PJ_area_2010:Road_1km_2019),
    by = "Grid"
  )

Cov_grid <- abind::abind(
  Dev_bg = (cov_grid$Development == "BG") %>% as.integer %>%
    array(dim = c(length(grid.list), length(years))),
  Dev_lo = (cov_grid$Development == "LO") %>% as.integer %>%
    array(dim = c(length(grid.list), length(years))),
  PJ_area = cov_grid %>% select(PJ_area_2010:PJ_area_2019) %>%
    data.matrix,
  NDVI = cov_grid %>% select(NDVI_2010:NDVI_2019) %>%
    data.matrix,
  WellA_1km = cov_grid %>% select(starts_with("Well_count_1km")) %>%
    data.matrix,
  WellA_3km = cov_grid %>% select(starts_with("Well_count_3km")) %>%
    data.matrix,
  WellD_1km = cov_grid %>% select(starts_with("PA_Well_count_1km")) %>%
    data.matrix,
  WellD_3km = cov_grid %>% select(starts_with("PA_Well_count_3km")) %>%
    data.matrix,
  Road_1km = cov_grid$Road_1km_2009 %>%
    cbind(cov_grid$Road_1km_2019 %>%
            matrix(nrow = length(grid.list), ncol = (length(years) - 1))),
  along = 3
)
dimnames(Cov_grid)[[2]] <- years

# Points #
cov_point <- data.frame(Grid = pointXyears.list %>% str_sub(1, -9),
                        Point = pointXyears.list %>% str_sub(1, -6),
                        stringsAsFactors = F) %>%
  mutate(gridIndex = Grid %>% as.factor %>% as.integer,
         pointIndex = Point %>% as.factor %>% as.integer,
         Development = ifelse(str_detect(Point, "ARIM-HI"), "HI",
                              ifelse(str_detect(Point, "ARIM-LO"), "LO", "BG"))) %>%
  distinct() %>%
  arrange(pointIndex) %>%
  dplyr::left_join(
    read.csv("covariates/ARIM_point_covariates_all_years.csv", stringsAsFactors = F, header = T) %>%
      tibble::as_tibble() %>%
      rename(Grid = TrnsctN) %>%
      mutate(Point = str_c(Grid,
                           str_pad(str_split(Waypont, "-", simplify = T)[,3], width = 2, side = "left", pad = "0"),
                           sep = "-"),
             across(Sage_2010:Herb_2019, function(x) suppressWarnings(as.numeric(x)))) %>%
      rename(Road_125m_2009 = Roads_Odonnell_2009_km,
             Road_125m_2019 = Road_length_2019_km,
             TPI_min = min_elev_diff,
             TPI_mean = mean_elev_diff,
             TPI_point = elev_diff) %>%
      select(Point, Well_count_2010:PA_Well_count_2019, Road_125m_2009, Road_125m_2019,
             Sage_2010:Herb_2019, vrm_125m:TPI_mean),
    by = "Point"
  )

Cov_point <- abind::abind(
  WellA_125m = cov_point %>% select(starts_with("Well_count")) %>%
    data.matrix,
  Road_125m = cov_point$Road_125m_2009 %>%
    cbind(cov_point$Road_125m_2009 %>%
            matrix(nrow = length(point.list), ncol = (length(years) - 1))),
  Sage = cov_point %>% select(starts_with("Sage")) %>%
    data.matrix(),
  Herb = cov_point %>% select(starts_with("Herb")) %>%
    data.matrix(),
  AHerb = cov_point %>% select(starts_with("AHerb")) %>%
    data.matrix(),
  vrm_125m = cov_point$vrm_125m %>%
    matrix(nrow = length(point.list), ncol = length(years)),
  TPI_min = cov_point$TPI_min %>%
    matrix(nrow = length(point.list), ncol = length(years)),
  TPI_mean = cov_point$TPI_mean %>%
    matrix(nrow = length(point.list), ncol = length(years)),
  TPI_point = cov_point$TPI_point %>%
    matrix(nrow = length(point.list), ncol = length(years)),
  along = 3
)
dimnames(Cov_point)[[2]] <- years

 # Missing value imputation for ground cover variables
library(randomForest)
library(QSLpersonal)
Cov_point_long <- Cov_point[,1,]
for(t in 2:length(years)) Cov_point_long <- rbind(Cov_point_long, Cov_point[,2,])
#sum(which(is.na(Cov_point)) != which(is.na(Cov_point_long))) # Should be zero.
grpID <- rep(cov_point$gridIndex, length(years)) * rep(years, each = nrow(cov_point))
vars <- c("Sage", "Herb", "AHerb")
for(v in vars) {
  xsum <- tapply(Cov_point_long[,v], grpID, function(x) mean(x, na.rm = T))[as.character(grpID)]
  dat <- cbind(Cov_point_long, xsum)
  dat <- Impute_missing_covs_rf(dat, v, dimnames(dat)[[2]][-which(dimnames(dat)[[2]] %in% vars)])
  Cov_point_long[,v] <- dat[,v]
}
Cov_point[which(is.na(Cov_point))] <- Cov_point_long[which(is.na(Cov_point))]
rm(Cov_point_long, grpID, vars, v, dat, xsum, t)

# Detection #
cov_pntyr <- data.frame(Point_year = pointXyears.list,
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
  distinct()

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

cov_pntyr <- cov_pntyr %>%
  left_join(
    grab %>% select(Point_year, CanCov, ShrubCov),
    by = "Point_year"
  ) %>%
  left_join(
    tab.datetime %>%
      select(Point_year, DOY, Time_ssr),
    by = "Point_year"
  )

rm(obs, maxDetPossible, sp, ss, tvec, grab)
save.image("Data_compiled.RData")
