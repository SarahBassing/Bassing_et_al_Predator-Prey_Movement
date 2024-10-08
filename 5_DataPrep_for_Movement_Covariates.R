  #'  ======================================================
  #'  Code from Bassing et al. 2024 "Predator-prey space-use  
  #'  and landscape features influence animal movement  
  #'  behaviors in a large-mammal community". Ecology.
  #'  
  #'  Covariate extraction for movement analysis
  #'  Washington Predator-Prey Project
  #'   =====================================================
  #'  Script to extract covariate data at each GPS collar location and interpolated
  #'  locations generated by the crawlWrap function. 
  #'  Relevant covariates include:
  #'    -Terrain Ruggedness Index (TRI; 30m res)
  #'    -Habitat openness (30m res)
  #'    -Snow cover (derived from GEE snow on and snow disappearance dates; 500m, daily res)
  #'    -Distance to nearest road (30m res)
  #'    -Prey habitat selection (30m res)
  #'    -Predator habitat selection (30m res)
  #'  
  #'  Notes about covariates:
  #'  Snow Disappearance Date & Snow On Date rasters were generated using Google 
  #'  Earth Engine (GEE) and scripts provided by T. Ganz & K. Breen. Snow 
  #'  Disappearance Date script adapted from code written & edited by R.Crumley 
  #'  (2019) & E.Mar (2018), and based on Nolin et al 2021. 
  #'  
  #'  Snow On Date script adapted from GEE code based on https://developers.google.com/earth-engine/tutorials/community/identifying-first-day-no-snow
  #'  
  #'  Distance to nearest road raster was generated by T. Ganz in ArcGIS using 
  #'  the Cascadia Biodiversity Watch roads layer. 
  #'  
  #'  Habitat openness raster was generated using the Cascadia Biodiveristy Watch 
  #'  landcover type raster. 
  #'  
  #'  Predator and prey habitat selection rasters were generated by estimating 
  #'  species- and season-specific resource selection functions and predicting 
  #'  relative probability of selection across study areas per year with the 
  #'  Resource_Selection_Function_Models.R script. Covariates informing RSFs are 
  #'  different than those included here.
  #'  
  #'  NOTE: This script is provided for transparency and reproducibility but it 
  #'  will not run with data provided on Dryad (animal relocation data are sensitive; 
  #'  contact Director of the Science Division with the Washington Dept. of Fish 
  #'  and Wildlife for raw data).
  #'  ======================================================
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(sf)
  library(stars)
  library(rgeos)
  library(raster)
  library(terra)
  library(parallel)
  library(doParallel)
  library(future.apply)
  library(tidyverse)
  
  #'  Load crwOut animal location data for each species (takes a hot minute)
  load("./Outputs/Telemetry_crwOut/crwOut_ALL.RData")
  
  #'  Read in spatial data
  TRI <- rast("./Shapefiles/WA DEM rasters/WPPP_TRI.tif")
  dist2Road_major <- rast("./Shapefiles/Cascadia_layers/roadsForTaylor/WPPP_dist_to_major_rds_2855.tiff")
  dist2Road_minor <- rast("./Shapefiles/Cascadia_layers/roadsForTaylor/WPPP_dist_to_minor_rds_2855.tiff")
  percopen18 <- rast("./Shapefiles/Cascadia_layers/PercOpen18.tif")
  percopen19 <- rast("./Shapefiles/Cascadia_layers/PercOpen19.tif")
  percopen20 <- rast("./Shapefiles/Cascadia_layers/PercOpen20.tif")
  snow_on18 <- rast("./Shapefiles/Snow/SOD_WPPP_2018_500m.tif")
  snow_on19 <- rast("./Shapefiles/Snow/SOD_WPPP_2019_500m.tif")
  snow_on20 <- rast("./Shapefiles/Snow/SOD_WPPP_2020_500m.tif")
  snow_off1819 <- rast("./Shapefiles/Snow/SDD2019_WPPP_500m.tif")
  snow_off1920 <- rast("./Shapefiles/Snow/SDD2020_WPPP_500m.tif")
  snow_off2021 <- rast("./Shapefiles/Snow/SDD2021_WPPP_500m.tif")
  #'  RSFs: stack of 3 rasters, 1 per year
  md_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/md_smr_RSFstack_global.tif") 
  md_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/md_wtr_RSFstack_global.tif")
  elk_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/elk_smr_RSFstack_global.tif")
  elk_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/elk_wtr_RSFstack_global.tif")
  wtd_smr_rsf <- rast("./Shapefiles/Predicted_RSFs/wtd_smr_RSFstack_global.tif")
  wtd_wtr_rsf <- rast("./Shapefiles/Predicted_RSFs/wtd_wtr_RSFstack_global.tif")
  coug_smr_OK_rsf <- rast("./Shapefiles/Predicted_RSFs/coug_smr_OK_RSFstack.tif")
  coug_smr_NE_rsf <- rast("./Shapefiles/Predicted_RSFs/coug_smr_NE_RSFstack.tif")
  coug_wtr_OK_rsf <- rast("./Shapefiles/Predicted_RSFs/coug_wtr_OK_RSFstack.tif")
  coug_wtr_NE_rsf <- rast("./Shapefiles/Predicted_RSFs/coug_wtr_NE_RSFstack.tif")
  wolf_smr_OK_rsf <- rast("./Shapefiles/Predicted_RSFs/wolf_smr_OK_RSFstack.tif")
  wolf_smr_NE_rsf <- rast("./Shapefiles/Predicted_RSFs/wolf_smr_NE_RSFstack.tif")
  wolf_wtr_OK_rsf <- rast("./Shapefiles/Predicted_RSFs/wolf_wtr_OK_RSFstack.tif")
  wolf_wtr_NE_rsf <- rast("./Shapefiles/Predicted_RSFs/wolf_wtr_NE_RSFstack.tif")
  
  #'  Identify projections & resolutions of spatial layers
  crs(TRI, describe=TRUE, proj=TRUE) 
  crs(dist2Road_major, describe=TRUE, proj=TRUE) 
  crs(percopen18, describe=TRUE, proj=TRUE) 
  crs(snow_on18, describe=TRUE, proj=TRUE) 
  crs(snow_off1819, describe=TRUE, proj=TRUE) 
  crs(md_smr_rsf, describe=TRUE, proj=TRUE) 
  
  res(TRI)
  res(dist2Road_major)
  res(percopen18)
  res(snow_on18)
  res(snow_off1819)
  res(md_smr_rsf)
  
  #'  Stack similar covariate rasters (use c() to stack with the terra package)
  road_stack <- c(dist2Road_major, dist2Road_minor)
  open_stack <- c(percopen18, percopen19, percopen20)
  snowOn_stack <- c(snow_on18, snow_on19, snow_on20) 
  snowOff_stack <- c(snow_off1819, snow_off1920, snow_off2021)
  
  #'  Define projections based on raster projections
  sa_proj <- projection("EPSG:2855")
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Make crwOut data spatial sf objects
  spatial_locs <- function(locs) {
    move <- locs[[2]]
    sf_locs <- st_as_sf(move, coords = c("mu.x", "mu.y"), crs = sa_proj) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      dplyr::select(c("ID", "FullID", "Sex", "Season", "StudyArea", "AnimalID", "time", "obs")) %>%
      #'  Fill in missing values for AnimalID, etc. where locations interpolated
      fill(FullID, .direction = "down") %>%
      fill(Sex, .direction = "down") %>%
      fill(Season, .direction = "down") %>%
      fill(StudyArea, .direction = "down") %>%
      fill(AnimalID, .direction = "down")
    return(sf_locs)
  }
  sf_locs <- lapply(crwOut_ALL, spatial_locs)
  

  #'  -------------------------------------------
  ####  COVARIATE EXTRACTION & CALCULATIONS  ####
  #'  -------------------------------------------
  #'  Takes quite awhile. Go get a cup of coffee.
  
  #'  Monitor time
  start.time <- Sys.time()

  #'  Master function to extract and manipulate covariate data for each species
  cov_extract <- function(locs) { 
    
    #'  Reproject sf objects to match WGS84 rasters to match certain rasters
    locs_wgs84 <- st_transform(locs, crs = wgs84)
    
    #'  --------------------------------------------------
    ####  1. Extract data from habitat and road rasters  #### 
    #'  --------------------------------------------------
    #'  Note: most of these are stacked rasters, 1 raster per year
    #'  Be sure to match data and raster projections
    #'  Extract from terra package requires location data be vectorized (vect())
    tri <- terra::extract(TRI, vect(locs_wgs84)) %>% 
      transmute(obs = ID, TRI = WPPP_TRI)
    dist2Road <- terra::extract(road_stack, vect(locs))
    #'  Find distance to nearest road between major and minor roadways- takes a minute
    #'  apply min function row-wise (1 = across each row), excluding 1st column of df
    dist2Road$NearestRd <- apply(dist2Road[-1], 1, min)
    dist2Road <- transmute(dist2Road, obs = ID, Dist2Road = NearestRd)
    percopen <- terra::extract(open_stack, vect(locs_wgs84)) %>%
      transmute(obs = ID, PercOpen18 = PercOpen18, PercOpen19 = PercOpen19, PercOpen20 = PercOpen20) %>%
      full_join(locs, by = "obs") %>%
      mutate(
        PercOpen = ifelse(Season == "Summer18", PercOpen18, PercOpen20),
        PercOpen = ifelse(Season == "Summer19", PercOpen19, PercOpen),
        PercOpen = ifelse(Season == "Winter1819", PercOpen18, PercOpen),
        PercOpen = ifelse(Season == "Winter1920", PercOpen19, PercOpen)
      ) %>%
      dplyr::select(c(obs, ID, AnimalID, Season, StudyArea, PercOpen)) 

    #'  ----------------------------------------------------------------
    ####  2. Extract data from snow-on and snow-disappearance rasters  ####
    #'  ----------------------------------------------------------------
    #'  Need to extract snow-on and snow-off dates for each location, then compare
    #'  this period of snow cover to timing of animal location and determine if
    #'  there was > 0% snow cover at the location when the animal was present.
    
    #'  Snow-On-Date (based on first day looking backwards from Dec. 31 that 
    #'  pixel has 0% snow cover)
    snowon <- terra::extract(snowOn_stack, vect(locs_wgs84)) 
    names(snowon) <- c("obs", "calDoy18", "calDoy19", "calDoy20") 
    snowon <- snowon %>%
      mutate(
        #'  Some snow-on observations are "zero" due to MODIS masking cloud cover, etc. 
        #'  Apply a location's mean calDoy using rowMean across the three years 
        #'  of sampling when observation is missing (0) for a single year 
        calDoy18 = ifelse(calDoy18 == 0, round(rowMeans(snowon[,c(3:4)], na.rm = TRUE), 0), calDoy18),
        calDoy19 = ifelse(calDoy19 == 0, round(rowMeans(snowon[,c(2,4)], na.rm = TRUE), 0), calDoy19),
        calDoy20 = ifelse(calDoy20 == 0, round(rowMeans(snowon[,c(2:3)], na.rm = TRUE), 0), calDoy20),
        #'  Change 0s to NAs so missing observations don't bias averaging below
        calDoy18 = ifelse(calDoy18 == 0, NA, calDoy18),
        calDoy19 = ifelse(calDoy19 == 0, NA, calDoy19),
        calDoy20 = ifelse(calDoy20 == 0, NA, calDoy20))
    #'  Apply annual mean calDoy using colMeans (& round to nearest whole number) 
    #'  when calDoy is masked across all years (all 0s) 
    #'  FYI, colMean includes missing observations in averaging if this step is
    #'  done in the same mutate command above. So annoying!
    snowon <- snowon %>%
      mutate(
        SnowOn18 = ifelse(is.na(calDoy18), round(colMeans(snowon[c("calDoy18")], na.rm = TRUE), 0), calDoy18),
        SnowOn19 = ifelse(is.na(calDoy19), round(colMeans(snowon[c("calDoy19")], na.rm = TRUE), 0), calDoy19),
        SnowOn20 = ifelse(is.na(calDoy20), round(colMeans(snowon[c("calDoy20")], na.rm = TRUE), 0), calDoy20),
        #'  calDoy represents the last julian day with NO snow (origin day is Jan. 1)
        #'  Need to add 1 to calDoy to represent first julian day WITH > 0% snow
        SOD18 = SnowOn18 + 1,
        SOD19 = SnowOn19 + 1,
        SOD20 = SnowOn20 + 1, 
        #'  Convert SOD (days since Jan. 1) to actual dates
        #'  Make into character b/c ifelse below does odd formatting otherwise
        SOD18_Date = as.character(as.Date(SOD18, origin = "2018-01-01")),
        SOD19_Date = as.character(as.Date(SOD19, origin = "2019-01-01")),
        SOD20_Date = as.character(as.Date(SOD20, origin = "2020-01-01"))) 

    #'  Snow-Disappearance-Date (first julian day moving forward where pixel had
    #'  0% snow cover, origin day is Oct. 1 based on "Water Year")
    snowoff <- terra::extract(snowOff_stack, vect(locs_wgs84))
    names(snowoff) <- c("obs", "SDD1819", "SDD1920", "SDD2021")
    snowoff <- snowoff %>%
      mutate(
        #'  Some snow disappearance dates are missing due to MODIS masking cloud cover, etc.
        #'  Change 0s to NAs so missing observations don't bias averaging below
        SDD1819 = ifelse(SDD1819 == 0, NA, SDD1819),
        SDD1920 = ifelse(SDD1920 == 0, NA, SDD1920),
        SDD2021 = ifelse(SDD2021 == 0, NA, SDD2021))
    #'  Apply annual mean calDoy using colMeans (& round to nearest whole number) 
    #'  to observations where calDoy was masked. Not using row-wise averages for 
    #'  a given location like above b/c many observations were missing for 2 years 
    #'  which biases row-wise average low using rowMeans.
    #'  FYI, colMean includes missing observations in averaging if this step is
    #'  done in the same mutate command above. So annoying!
    snowoff <- snowoff %>%
      mutate(
        SnowOff1819 = ifelse(is.na(SDD1819), round(colMeans(snowoff[c("SDD1819")], na.rm = TRUE), 0), SDD1819),
        SnowOff1920 = ifelse(is.na(SDD1920), round(colMeans(snowoff[c("SDD1920")], na.rm = TRUE), 0), SDD1920),
        SnowOff2021 = ifelse(is.na(SDD2021), round(colMeans(snowoff[c("SDD2021")], na.rm = TRUE), 0), SDD2021),
        #'  Convert SDD (days since Oct. 1) to actual dates
        #'  Note origin date is start of "Water Year", not January 1
        #'  Make into character b/c ifelse below does odd formatting otherwise
        SDD1819_Date = as.character(as.Date(SnowOff1819, origin = "2018-10-01")),
        SDD1920_Date = as.character(as.Date(SnowOff1920, origin = "2019-10-01")),
        SDD2021_Date = as.character(as.Date(SnowOff2021, origin = "2020-10-01"))) 
      
    #'  Identify whether there was > 0% snow cover at each animal location
    snowcover <- snowon %>%
      full_join(snowoff, by = "obs") %>%
      #'  Save only dates and arrange by snow on & disappearance dates per year
      dplyr::select(c(obs, SOD18_Date, SDD1819_Date, SOD19_Date, SDD1920_Date,
                      SOD20_Date, SDD2021_Date)) %>%
      #'  Join with crwOut data
      full_join(locs, by = "obs") %>% 
      mutate(
        Date = as.Date(time, tz = "America/Los_Angeles"),
        #'  Create single column for snow-on dates across years
        SOD = ifelse(Season == "Summer18", SOD18_Date, SOD20_Date), # don't actually need Snow for summer locs but helpful for column book-keeping
        SOD = ifelse(Season == "Summer19", SOD19_Date, SOD),        
        SOD = ifelse(Season == "Winter1819", SOD18_Date, SOD),
        SOD = ifelse(Season == "Winter1920", SOD19_Date, SOD),
        SOD = as.Date(SOD, format = "%Y-%m-%d"),
        #'  Create single column for snow-disappearance dates across years
        SDD = ifelse(Season == "Summer18", SDD1819_Date, SDD2021_Date), # don't actually need Snow for summer locs but helpful for column book-keeping
        SDD = ifelse(Season == "Summer19", SDD1920_Date, SDD),          
        SDD = ifelse(Season == "Winter1819", SDD1819_Date, SDD),
        SDD = ifelse(Season == "Winter1920", SDD1920_Date, SDD),
        SDD = as.Date(SDD, format = "%Y-%m-%d"),
        #'  Indicate if location Date falls between snow-on & snow-disappearance date
        SnowCover = ifelse(Date >= SOD & Date < SDD, 1, 0)
      ) %>%
      dplyr::select(c(obs, ID, AnimalID, Season, StudyArea, SOD, Date, SDD, SnowCover))
    
    #'  ----------------------------------------------------------
    ####  3. Species and season specific RSFs (stacked by year)  ####
    #'  ----------------------------------------------------------
    #'  Note: NAs arise for locations if RSF is specific to only one study area
    
    ####  MULE DEER  ####
    md_smr <- terra::extract(md_smr_rsf, vect(locs))
    names(md_smr) <- c("obs", "MD_smr18", "MD_smr19", "MD_smr20")
    md_wtr <- terra::extract(md_wtr_rsf, vect(locs))
    names(md_wtr) <- c("obs", "MD_wtr1819", "MD_wtr1920", "MD_wtr2021")
    md_rsf <- full_join(md_smr, md_wtr, by = "obs") %>%
      full_join(locs, by = "obs") %>%
      mutate(
        MD_smr = ifelse(Season == "Summer18", MD_smr18, MD_smr20),
        MD_smr = ifelse(Season == "Summer19", MD_smr19, MD_smr),
        MD_wtr = ifelse(Season == "Winter1819", MD_wtr1819, MD_wtr2021),
        MD_wtr = ifelse(Season == "Winter1920", MD_wtr1920, MD_wtr)
      ) %>%
      dplyr::select(c(obs, ID, AnimalID, Season, StudyArea, MD_smr, MD_wtr)) 
    ####  ELK  ####
    elk_smr <- terra::extract(elk_smr_rsf, vect(locs))
    names(elk_smr) <- c("obs", "ELK_smr18", "ELK_smr19", "ELK_smr20")
    elk_wtr <- terra::extract(elk_wtr_rsf, vect(locs))
    names(elk_wtr) <- c("obs", "ELK_wtr1819", "ELK_wtr1920", "ELK_wtr2021")
    elk_rsf <- full_join(elk_smr, elk_wtr, by = "obs") %>%
      full_join(locs, by = "obs") %>%
      mutate(
        ELK_smr = ifelse(Season == "Summer18", ELK_smr18, ELK_smr20),
        ELK_smr = ifelse(Season == "Summer19", ELK_smr19, ELK_smr),
        ELK_wtr = ifelse(Season == "Winter1819", ELK_wtr1819, ELK_wtr2021),
        ELK_wtr = ifelse(Season == "Winter1920", ELK_wtr1920, ELK_wtr)
      ) %>%
      dplyr::select(c(obs, ID, AnimalID, Season, StudyArea, ELK_smr, ELK_wtr)) 
    ####  WHITE-TAILED DEER  ####
    wtd_smr <- terra::extract(wtd_smr_rsf, vect(locs))
    names(wtd_smr) <- c("obs", "WTD_smr18", "WTD_smr19", "WTD_smr20")
    wtd_wtr <- terra::extract(wtd_wtr_rsf, vect(locs))
    names(wtd_wtr) <- c("obs", "WTD_wtr1819", "WTD_wtr1920", "WTD_wtr2021")
    wtd_rsf <- full_join(wtd_smr, wtd_wtr, by = "obs") %>%
      full_join(locs, by = "obs") %>%
      mutate(
        WTD_smr = ifelse(Season == "Summer18", WTD_smr18, WTD_smr20),
        WTD_smr = ifelse(Season == "Summer19", WTD_smr19, WTD_smr),
        WTD_wtr = ifelse(Season == "Winter1819", WTD_wtr1819, WTD_wtr2021),
        WTD_wtr = ifelse(Season == "Winter1920", WTD_wtr1920, WTD_wtr)
      ) %>%
      dplyr::select(c(obs, ID, AnimalID, Season, StudyArea, WTD_smr, WTD_wtr)) 
    ####  COUGAR  ####
    coug_smr_OK <- terra::extract(coug_smr_OK_rsf, vect(locs))
    names(coug_smr_OK) <- c("obs", "COUG_smr18_OK", "COUG_smr19_OK", "COUG_smr20_OK")
    coug_smr_NE <- terra::extract(coug_smr_NE_rsf, vect(locs))
    names(coug_smr_NE) <- c("obs", "COUG_smr18_NE", "COUG_smr19_NE", "COUG_smr20_NE")
    coug_wtr_OK <- terra::extract(coug_wtr_OK_rsf, vect(locs))
    names(coug_wtr_OK) <- c("obs", "COUG_wtr1819_OK", "COUG_wtr1920_OK", "COUG_wtr2021_OK")
    coug_wtr_NE <- terra::extract(coug_wtr_NE_rsf, vect(locs))
    names(coug_wtr_NE) <- c("obs", "COUG_wtr1819_NE", "COUG_wtr1920_NE", "COUG_wtr2021_NE")
    coug_rsf <- full_join(coug_smr_OK, coug_smr_NE, by = "obs") %>%
      full_join(coug_wtr_OK, by = "obs") %>%
      full_join(coug_wtr_NE, by = "obs") %>%
      full_join(locs, by = "obs") %>%
      mutate(
        #'  Keep in mind this approach provides a meaningless value for summer/winter
        #'  RSF when the Season is winter/summer - be sure to drop this column later
        COUG_smr_OK = ifelse(Season == "Summer18", COUG_smr18_OK, COUG_smr20_OK),
        COUG_smr_OK = ifelse(Season == "Summer19", COUG_smr19_OK, COUG_smr_OK),
        COUG_smr_NE = ifelse(Season == "Summer18", COUG_smr18_NE, COUG_smr20_NE),
        COUG_smr_NE = ifelse(Season == "Summer19", COUG_smr19_NE, COUG_smr_NE),
        COUG_wtr_OK = ifelse(Season == "Winter1819", COUG_wtr1819_OK, COUG_wtr2021_OK),
        COUG_wtr_OK = ifelse(Season == "Winter1920", COUG_wtr1920_OK, COUG_wtr_OK),
        COUG_wtr_NE = ifelse(Season == "Winter1819", COUG_wtr1819_NE, COUG_wtr2021_NE),
        COUG_wtr_NE = ifelse(Season == "Winter1920", COUG_wtr1920_NE, COUG_wtr_NE)
      ) %>%
      dplyr::select(c(obs, ID, AnimalID, Season, StudyArea, COUG_smr_OK, COUG_smr_NE, COUG_wtr_OK, COUG_wtr_NE)) 
    ####  WOLF  ####
    wolf_smr_OK <- terra::extract(wolf_smr_OK_rsf, vect(locs))
    names(wolf_smr_OK) <- c("obs", "WOLF_smr18_OK", "WOLF_smr19_OK", "WOLF_smr20_OK")
    wolf_smr_NE <- terra::extract(wolf_smr_NE_rsf, vect(locs))
    names(wolf_smr_NE) <- c("obs", "WOLF_smr18_NE", "WOLF_smr19_NE", "WOLF_smr20_NE")
    wolf_wtr_OK <- terra::extract(wolf_wtr_OK_rsf, vect(locs))
    names(wolf_wtr_OK) <- c("obs", "WOLF_wtr1819_OK", "WOLF_wtr1920_OK", "WOLF_wtr2021_OK")
    wolf_wtr_NE <- terra::extract(wolf_wtr_NE_rsf, vect(locs))
    names(wolf_wtr_NE) <- c("obs", "WOLF_wtr1819_NE", "WOLF_wtr1920_NE", "WOLF_wtr2021_NE")
    wolf_rsf <- full_join(wolf_smr_OK, wolf_smr_NE, by = "obs") %>%
      full_join(wolf_wtr_OK, by = "obs") %>%
      full_join(wolf_wtr_NE, by = "obs") %>%
      full_join(locs, by = "obs") %>%
      mutate(
        WOLF_smr_OK = ifelse(Season == "Summer18", WOLF_smr18_OK, WOLF_smr20_OK),
        WOLF_smr_OK = ifelse(Season == "Summer19", WOLF_smr19_OK, WOLF_smr_OK),
        WOLF_smr_NE = ifelse(Season == "Summer18", WOLF_smr18_NE, WOLF_smr20_NE),
        WOLF_smr_NE = ifelse(Season == "Summer19", WOLF_smr19_NE, WOLF_smr_NE),
        WOLF_wtr_OK = ifelse(Season == "Winter1819", WOLF_wtr1819_OK, WOLF_wtr2021_OK),
        WOLF_wtr_OK = ifelse(Season == "Winter1920", WOLF_wtr1920_OK, WOLF_wtr_OK),
        WOLF_wtr_NE = ifelse(Season == "Winter1819", WOLF_wtr1819_NE, WOLF_wtr2021_NE),
        WOLF_wtr_NE = ifelse(Season == "Winter1920", WOLF_wtr1920_NE, WOLF_wtr_NE)
      ) %>%
      dplyr::select(c(obs, ID, AnimalID, Season, StudyArea, WOLF_smr_OK, WOLF_smr_NE, WOLF_wtr_OK, WOLF_wtr_NE)) 
    
    #'  Join all extracted covariates together
    crwOut_covs <- full_join(tri, dist2Road, by = "obs") %>%
      full_join(percopen, by = "obs") %>%
      full_join(snowcover, by = c("obs", "ID", "AnimalID", "Season", "StudyArea")) %>% 
      full_join(md_rsf, by = c("obs", "ID", "AnimalID", "Season", "StudyArea")) %>%
      full_join(elk_rsf, by = c("obs", "ID", "AnimalID", "Season", "StudyArea")) %>%
      full_join(wtd_rsf, by = c("obs", "ID", "AnimalID", "Season", "StudyArea")) %>%
      full_join(coug_rsf, by = c("obs", "ID", "AnimalID", "Season", "StudyArea")) %>%
      full_join(wolf_rsf, by = c("obs", "ID", "AnimalID", "Season", "StudyArea")) %>%
      full_join(locs, by = c("obs", "ID", "AnimalID", "Season", "StudyArea")) %>%
      dplyr::select(c("obs", "ID", "AnimalID", "Season", "StudyArea", "time", 
                      "Date", "Dist2Road", "PercOpen", "SnowCover", "TRI", "MD_smr",       
                      "MD_wtr", "ELK_smr", "ELK_wtr", "WTD_smr", "WTD_wtr",
                      "COUG_smr_OK", "COUG_smr_NE", "COUG_wtr_OK", "COUG_wtr_NE", 
                      "WOLF_smr_OK", "WOLF_smr_NE", "WOLF_wtr_OK", "WOLF_wtr_NE"))
                      
    return(crwOut_covs)
  }
  #'  Run list of species location data through function 
  spp_telem_covs <- lapply(sf_locs, cov_extract) 
  
  #'  End time keeping
  end.time <- Sys.time()
  #'  How long did this take?
  difftime(end.time, start.time, units = "hours")
  
  
  #'  List covariate data sets by season and study area (for predators)
  #'  Note the order: 1) MD smr, 2) MD wtr, 3) ELK smr, 4) ELK wtr, 5) WTD smr, 
  #'  6) WTD wtr, 7) COUG smr OK, 8) COUG wtr OK, 9) COUG smr NE, 10) COUG wtr NE, 
  #'  11) WOLF smr OK, 12) WOLF wtr OK, 13) WOLF smr NE, 14) WOLF wtr NE
  smr_covs <- list(spp_telem_covs[[1]], spp_telem_covs[[3]], spp_telem_covs[[5]], 
                   spp_telem_covs[[7]], spp_telem_covs[[9]], spp_telem_covs[[11]], 
                   spp_telem_covs[[13]])
  wtr_covs <- list(spp_telem_covs[[2]], spp_telem_covs[[4]], spp_telem_covs[[6]], 
                   spp_telem_covs[[8]], spp_telem_covs[[10]], spp_telem_covs[[12]], 
                   spp_telem_covs[[14]])
  
  #'  Remove summer/winter RSF values from winter/summer data sets, respectively.
  #'  Be aware there are some locations where landscape covariates were extracted
  #'  but missing RSF values b/c those locations were masked in the RSF analyses,
  #'  producing NAs for all RSF values at these locations
  remove_wtr_covs <- function(covs) {
    smr_data <- dplyr::select(covs, -c("MD_wtr", "ELK_wtr", "WTD_wtr", "COUG_wtr_OK", 
                                       "COUG_wtr_NE", "WOLF_wtr_OK", "WOLF_wtr_NE")) %>% 
      #'  Add columns for location hour, re-factored hour, & whether it's day or 
      #'  night based on sunrise/sunset times over course of each month
      mutate(hour = as.integer(strftime(time, format = "%H", tz="Etc/GMT+8")),
             hour_fix = ifelse(hour == 0 | hour == 2, 2, 22),
             hour_fix = ifelse(hour == 4 | hour == 6, 6, hour_fix),
             hour_fix = ifelse(hour == 8 | hour == 10, 10, hour_fix),
             hour_fix = ifelse(hour == 12 | hour == 14, 14, hour_fix),
             hour_fix = ifelse(hour == 16 | hour  == 18, 18, hour_fix),
             # hour2 = as.numeric(as.factor(hour)),
             hour3 = ifelse(hour <= 2, 1, 6),
             hour3 = ifelse(hour >= 4 & hour <= 6, 2, hour3),
             hour3 = ifelse(hour >= 8 & hour <= 10, 3, hour3),
             hour3 = ifelse(hour >= 12 & hour <= 14, 4, hour3),
             hour3 = ifelse(hour >= 16 & hour <= 18, 5, hour3),
             month = as.integer(strftime(time, format = "%m", tz="Etc/GMT+8")),
             #' For July - Aug, hours between 5am and 9pm are daytime (1), else nightime (0)
             daytime = ifelse(month < 9 & hour < 5 | hour > 21, 1, 0),
             #'  For Sept, hours between 7am and 7pm are daytime (1), else nightime (0)
             daytime = ifelse(month == 9 & hour < 7 | hour > 19, 1, daytime))
    names(smr_data) <- c("obs", "ID", "AnimalID", "Season", "StudyArea", "time", "Date", 
                     "Dist2Road", "PercOpen", "SnowCover", "TRI", "MD_RSF", "ELK_RSF", 
                     "WTD_RSF", "COUG_RSF_OK","COUG_RSF_NE", "WOLF_RSF_OK", "WOLF_RSF_NE",    
                     "hour", "hour_fix", "hour3", "month", "daytime")
    return(smr_data)
  }
  smr_telem_data <- lapply(smr_covs, remove_wtr_covs)
  
  
  remove_smr_covs <- function(covs) {
    wtr_data <- dplyr::select(covs, -c("MD_smr", "ELK_smr", "WTD_smr", "COUG_smr_OK", 
                                       "COUG_smr_NE", "WOLF_smr_OK", "WOLF_smr_NE")) %>% 
      #'  Add columns for location hour, re-factored hour, & whether it's day or 
      #'  night based on sunrise/sunset times over course of each month
      mutate(hour = as.integer(strftime(time, format = "%H", tz="Etc/GMT+8")),
             hour_fix = ifelse(hour == 0 | hour == 2, 2, 22),
             hour_fix = ifelse(hour == 4 | hour == 6, 6, hour_fix),
             hour_fix = ifelse(hour == 8 | hour == 10, 10, hour_fix),
             hour_fix = ifelse(hour == 12 | hour == 14, 14, hour_fix),
             hour_fix = ifelse(hour == 16 | hour  == 18, 18, hour_fix),
             # hour2 = as.numeric(as.factor(hour)),
             hour3 = ifelse(hour <= 2, 1, 6),
             hour3 = ifelse(hour >= 4 & hour <= 6, 2, hour3),
             hour3 = ifelse(hour >= 8 & hour <= 10, 3, hour3),
             hour3 = ifelse(hour >= 12 & hour <= 14, 4, hour3),
             hour3 = ifelse(hour >= 16 & hour <= 18, 5, hour3),
             month = as.integer(strftime(time, format = "%m", tz="Etc/GMT+8")),
             #' For Dec, hours between 7:30am and 4:20pm are daytime (1), else nightime (0)
             daytime = ifelse(month == 12 & hour < 7 | hour > 16, 1, 0),
             #' For Jan, hours between 8am and 4:20pm are daytime (1), else nightime (0)
             daytime = ifelse(month == 1 & hour < 8 | hour > 16, 1, daytime),
             #'  For Feb, hours between 7:30am and 5pm are daytime (1), else nightime (0)
             daytime = ifelse(month == 2 & hour < 7 | hour > 17, 1, daytime))
    names(wtr_data) <- c("obs", "ID", "AnimalID", "Season", "StudyArea", "time", "Date", 
                     "Dist2Road", "PercOpen", "SnowCover", "TRI", "MD_RSF", "ELK_RSF", 
                     "WTD_RSF", "COUG_RSF_OK", "COUG_RSF_NE", "WOLF_RSF_OK", "WOLF_RSF_NE", 
                     "hour", "hour_fix", "hour3", "month", "daytime")
    return(wtr_data)
  }
  wtr_telem_data <- lapply(wtr_covs, remove_smr_covs)
  
  #'  List covariates by study area
  #'  Include only species collared in respective study area
  #'  Okanogan index: 1) MD, 4) COUG, 6) WOLF
  OK_covs <- list(smr_telem_data[[1]], wtr_telem_data[[1]], smr_telem_data[[4]], wtr_telem_data[[4]],
                  smr_telem_data[[6]], wtr_telem_data[[6]])
  #'  Northeast index: 2) ELK, 3) WTD, 5) COUG, 7) WOLF
  NE_covs <- list(smr_telem_data[[2]], wtr_telem_data[[2]], smr_telem_data[[3]], wtr_telem_data[[3]],
                  smr_telem_data[[5]], wtr_telem_data[[5]], smr_telem_data[[7]], wtr_telem_data[[7]])
  
  #'  Remove MD vs ELK/WTD and COUG/WOLF OK vs NE data from NE & OK covariates, respectively
  #'  These columns of pure NAs will cause problems when merging covariate data
  #'  with crwOut data in HMM script
  #'  Remove RSF data for species collared only in NE from OK data sets
  nix_NE_RSFs <- function(OK_covs) {
    no_NE_RSFs <- dplyr::select(OK_covs, -c(ELK_RSF, WTD_RSF, COUG_RSF_NE, WOLF_RSF_NE)) %>%
      #'  Rename study areas-specific cougar and wolf columns
      rename(COUG_RSF = COUG_RSF_OK,
             WOLF_RSF = WOLF_RSF_OK)
    return(no_NE_RSFs)
  }
  skinny_OK <- lapply(OK_covs, nix_NE_RSFs)
  
  ats_covs <- lapply(ats_telem_data, nix_NE_RSFs)
  
  #'  Remove RSF data for species collared only in OK from NE data sets
  nix_OK_RSFs <- function(NE_covs) {
    no_OK_RSFs <- dplyr::select(NE_covs, -c(MD_RSF, COUG_RSF_OK, WOLF_RSF_OK)) %>%
      #'  Rename study areas-specific cougar and wolf columns
      rename(COUG_RSF = COUG_RSF_NE,
             WOLF_RSF = WOLF_RSF_NE)
    return(no_OK_RSFs)
  }
  skinny_NE <- lapply(NE_covs, nix_OK_RSFs)
  
  
  #'  Merge back into one giant list keeping track of list order!
  #'  1) MD smr OK,    2) MD wtr OK,    3) ELK smr NE,   4) ELK wtr NE,   5) WTD smr NE, 
  #'  6) WTD wtr NE,   7) COUG smr OK,  8) COUG wtr OK,  9) COUG smr NE,  10) COUG wtr NE, 
  #'  11) WOLF smr OK, 12) WOLF wtr OK, 13) WOLF smr NE, 14) WOLF wtr NE
  spp_telem_covs <- list(skinny_OK[[1]], skinny_OK[[2]], skinny_NE[[1]], skinny_NE[[2]],
                         skinny_NE[[3]], skinny_NE[[4]], skinny_OK[[3]], skinny_OK[[4]], 
                         skinny_NE[[5]], skinny_NE[[6]], skinny_OK[[5]], skinny_OK[[6]],
                         skinny_NE[[7]], skinny_NE[[8]])

  
  #'  Save!
  save(spp_telem_covs, file = "./Outputs/Telemetry_covs/spp_telem_covs.RData")
  
  
  #'  Next stop: Hidden_Markov_Movement_Models.R to run HMMs

  