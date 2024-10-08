  #'  =====================================================
  #'  Code from Bassing et al. 2024 "Predator-prey space-use  
  #'  and landscape features influence animal movement  
  #'  behaviors in a large-mammal community". Ecology.
  #'  
  #'  Data Prep for Resource Selection Functions (RSFs)
  #'  Washington Predator-Prey Project
  #'  =====================================================
  #'  Script to randomly sample "available" points from the home range of each 
  #'  collared individual and extract covariate data at each "used" and "available"
  #'  animal location. Final output to be used as input data for RSFs.  
  #'  
  #'  NOTE: This script is provided for transparency and reproducibility but it 
  #'  will not run with data provided on Dryad (animal relocation data are sensitive; 
  #'  contact Director of the Science Division with the Washington Dept. of Fish 
  #'  and Wildlife for raw data). 
  #'  ======================================================
  
  #'  Clear memory
  rm(list = ls())  

  #'  Load packages for selecting available points
  library(tidyverse)
  library(sp)
  library(sf)
  library(rgeos)
  library(raster)
  library(lme4)
  library(adehabitatHR)
  library(rlist)
  
  #'  Load telemetry data
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDis_noMig.RData")
  #'  Relocation data excludes obvious dispersal events and relocation during
  #'  mule deer migrations because these long-distance movements are very different
  #'  from the typical within-home range movements of interest to this studying.
  
  #'  Functions to filter species-specific data by study area
  NE_collars <- function(locs) {
    tracks <- filter(locs, StudyArea == "NE")
    return(tracks)
  }
  OK_collars <- function(locs) {
    tracks <- filter(locs, StudyArea == "OK")
    return(tracks)
  }
  
  #'  Drop mule deer lists (no collared mule deer in NE study area)
  #'  Drop elk and white-tailed deer lists (no collared elk or wtd in OK study area)
  noMD <- spp_all_tracks[-c(1:2)]
  noELKWTD <- spp_all_tracks[-c(3:6)]
  
  #'  Run lists through respective functions
  NE_tracks <- lapply(noMD, NE_collars)
  OK_tracks <- lapply(noELKWTD, OK_collars)

  
  #'  Prep lists to generate available locations per individual
  #'  =========================================================
  #'  Function to pull out unique animal IDs
  unq_id <- function(locs) {
    animal <- as.data.frame(locs$AnimalID) %>%
      unique()
    colnames(animal) <- "AnimalID"
    return(animal)
  }
  animal_ID <- lapply(spp_all_tracks, unq_id)
  NE_ID <- lapply(NE_tracks, unq_id)
  OK_ID <- lapply(OK_tracks, unq_id)
  
  #'  Function to split data into list of dataframes by unique animal ID and year
  #'  (Use AnimalID if ignoring year aspect of data)
  split_animal <- function(locs) {
    ind_animal <- group_split(locs, locs$FullID) 
    return(ind_animal)
  }
  animal_split <- lapply(spp_all_tracks, split_animal)
  NE_split <- lapply(NE_tracks, split_animal)
  OK_split <- lapply(OK_tracks, split_animal)
  
  
  #'  Generate random "available" locations for each individual
  #'  =========================================================
  #'  Set projection for spatial analyses
  sa_proj <- projection("EPSG:2855")  # NAD83(HARN) / Washington North
  
  #'  Load water bodies shapefile (needed to mask unavailable habitat)
  waterbody <- sf::st_read("./Shapefiles/WA_DeptEcology_HydroWA", layer = "WPPP_waterbody") %>%
    st_transform(crs = sa_proj)
  #'  Identify large bodies of water (anything larger than 1 sq-km in size)
  bigwater <- waterbody[waterbody$AreSqKm > 1,]
  #'  Convert to an sp object instead of sf
  bigwater_sp <- as(st_geometry(bigwater), "Spatial")

  
  #'  3rd ORDER SELECTION Average 1:20 used:available
  #'  -----------------------------------------------
  #'  Function to randomly select "available" points within each animal's seasonal 
  #'  home range (utilization distributions)
  avail_pts_3rd <- function(locs, ex, plotit = F) {
    #'  1. Make each animal's locations spatial
    #'  ---------------------------------------------------------
    locs_sp <- SpatialPoints(locs[,c("x", "y")], proj4string = CRS(sa_proj)) 
    #'  Plot to make sure step 1 worked
    if(plotit) {
      plot(locs_sp, col = "blue", pch = 19, main = locs$FullID)
    }
    
    #'  2. Create UDs for each animal following Bivariate normal utilization distributions
    #'  ----------------------------------------------------------
    #'  Estimate KDEs for each home range, extend the spatial extent by 1.5 (or 2) 
    #'  when estimating UDs (throws an error about grid being too small to 
    #'  estimate home range otherwise)
    UD <- kernelUD(locs_sp, kern = "bivnorm", extent = ex) 
    UD50 <- getverticeshr(UD, 50)
    UD75 <- getverticeshr(UD, 75)
    UD90 <- getverticeshr(UD, 90)
    UD95 <- getverticeshr(UD, 95)
    #'  Intersect and clip water body polygons from MCP polygons so large bodies
    #'  of water are not available to collared animals
    UD95_clip <- rgeos::gDifference(UD95, bigwater_sp)
    #'  Plot to make sure step 2 worked
    if(plotit) {
      plot(UD95, border = "darkgreen", col = NA)
      plot(locs_sp, col = "blue", pch = 19, cex = 0.75, add = T)
      plot(UD50, border = "green", col = NA, add = T)
      plot(UD75, border = "green", col = NA, add = T)
      plot(UD90, border = "green", col = NA, add = T)
      plot(UD95_clip, border = "red", col = NA, add = T)
    }
    
    #'  3. Randomly select points within each home range
    #'  ------------------------------------------------
    #'  Sampling 20 available points to every 1 used point
    #'  Identify number of used points per individual
    nused <- nrow(locs) 
    #'  Multiply by 20 
    navailable <- nused*20
    rndpts <- spsample(UD95_clip, navailable, type = "random")
    #'  Turn them into spatial points
    rndpts_sp <- SpatialPoints(rndpts, proj4string = CRS(sa_proj))
    #' Plot to make sure step 3 worked
    if(plotit) {
      plot(rndpts_sp, col = "red", pch = 19, cex = 0.75)
    }
    
    #'  4. Make list of locations non-spatial
    #'  -------------------------------------
    rndpts_df <- as.data.frame(rndpts_sp) 
    ID <- unique(droplevels(locs$AnimalID))
    Season <- unique(locs$Season)
    rndpts_df$ID <- ID
    rndpts_df$Season <- Season
    
    return(rndpts_df)

  }
  #'  Run lists by species and season through function
  md_smr_df <- lapply(animal_split[[1]], avail_pts_3rd, ex = 1.5, T) 
  md_wtr_df <- lapply(animal_split[[2]], avail_pts_3rd, ex = 1.5, T) 
  elk_smr_df <- lapply(animal_split[[3]], avail_pts_3rd, ex = 1.5, T) 
  elk_wtr_df <- lapply(animal_split[[4]], avail_pts_3rd, ex = 1.5, T) 
  wtd_smr_df <- lapply(animal_split[[5]], avail_pts_3rd, ex = 1.5, T) 
  wtd_wtr_df <- lapply(animal_split[[6]], avail_pts_3rd, ex = 2, T) # Note the different extent- fails otherwise
  coug_smr_df <- lapply(animal_split[[7]], avail_pts_3rd, ex = 1.5, T) 
  coug_wtr_df <- lapply(animal_split[[8]], avail_pts_3rd, ex = 1.5, T) 
  wolf_smr_df <- lapply(animal_split[[9]], avail_pts_3rd, ex = 1.5, T) 
  wolf_wtr_df <- lapply(animal_split[[10]], avail_pts_3rd, ex = 1.5, T) 
  
  #'  Convert to dataframes instead of lists
  md_smr_df <- do.call(rbind.data.frame, md_smr_df)
  md_wtr_df <- do.call(rbind.data.frame, md_wtr_df)
  elk_smr_df <- do.call(rbind.data.frame, elk_smr_df)
  elk_wtr_df <- do.call(rbind.data.frame, elk_wtr_df)
  wtd_smr_df <- do.call(rbind.data.frame, wtd_smr_df)
  wtd_wtr_df <- do.call(rbind.data.frame, wtd_wtr_df)
  coug_smr_df <- do.call(rbind.data.frame, coug_smr_df)
  coug_wtr_df <- do.call(rbind.data.frame, coug_wtr_df)
  wolf_smr_df <- do.call(rbind.data.frame, wolf_smr_df)
  wolf_wtr_df <- do.call(rbind.data.frame, wolf_wtr_df)
  
  #'  Gather into one big list
  md_available <- list(md_smr_df, md_wtr_df)
  elk_available <- list(elk_smr_df, elk_wtr_df)
  wtd_available <- list(wtd_smr_df, wtd_wtr_df)
  coug_available <- list(coug_smr_df, coug_wtr_df)
  wolf_available <- list(wolf_smr_df, wolf_wtr_df)
  
  #'  Extract covariate data for each available location
  #'  ==================================================
  #'  This will take awhile!

  #'  Load packages for covariate extraction
  library(sf)
  library(sp)
  library(rgeos)
  library(raster)
  library(parallel)
  library(future.apply)
  library(tidyverse)
  
  #'  Read in spatial data
  #'  Terrain rasters
  dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m.tif")
  slope <- raster("./Shapefiles/WA DEM rasters/WPPP_slope_aspect.tif", band = 1)
  #'  Canopy cover change
  canopy18 <- raster("./Shapefiles/Global_Forest_Change/treecov_2018.tif")
  canopy19 <- raster("./Shapefiles/Global_Forest_Change/treecov_2019.tif")
  canopy20 <- raster("./Shapefiles/Global_Forest_Change/treecov_2020.tif")
  #'  Cascadia Biodiveristy Watch rasters & shapefiles
  landcov18 <- raster("./Shapefiles/Cascadia_layers/landcover_2018.tif")
  landcov19 <- raster("./Shapefiles/Cascadia_layers/landcover_2019.tif")
  landcov20 <- raster("./Shapefiles/Cascadia_layers/landcover_2020.tif")
  dist_open18 <- raster("./Shapefiles/Cascadia_layers/Dist2OpenEdge18.tif")
  dist_open19 <- raster("./Shapefiles/Cascadia_layers/Dist2OpenEdge19.tif")
  dist_open20 <- raster("./Shapefiles/Cascadia_layers/Dist2OpenEdge20.tif")
  dist_close18 <- raster("./Shapefiles/Cascadia_layers/Dist2ForestEdge18.tif")
  dist_close19 <- raster("./Shapefiles/Cascadia_layers/Dist2ForestEdge19.tif")
  dist_close20 <- raster("./Shapefiles/Cascadia_layers/Dist2ForestEdge20.tif")
  rdden <- raster("./Shapefiles/Cascadia_layers/roadsForTaylor/RoadDensity_1km.tif")
  #'  Distance to water
  h2o <- raster("./Shapefiles/WA_DeptEcology_HydroWA/WPPP_dist_to_water_2855.tif")
  #'  Human modified landscapes
  HM <- raster("./Shapefiles/Additional_WPPP_Layers/WPPP_gHM.tif")
  
  #'  Create raster stacks of data with same dimensions and resolution
  terra_stack <- stack(dem, slope)
  dist2edge_stack18 <- stack(dist_close18, dist_open18) # includes open & closed habitats
  dist2edge_stack19 <- stack(dist_close19, dist_open19) 
  dist2edge_stack20 <- stack(dist_close20, dist_open20) 
  
  #'  Note the different projections
  projection(dem)
  projection(canopy20)
  projection(landcov18)
  projection(rdden)
  projection(h2o)
  projection(dist_close18) # slightly different from sa_proj
  projection(dist_close20)
  
  #'  Identify projections & resolutions of relevant features
  sa_proj <- projection("EPSG:2855")  # NAD83(HARN) / Washington North
  dist_proj <- projection("+proj=lcc +lat_0=47 +lon_0=-120.833333333333 +lat_1=47.5 +lat_2=48.7333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs")
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  

  #'  Function to make available points data a spatial sf object
  #'  Actually make it a SpatialPointsDF with sp if running in parallel
  #'  No clue why it won't work as an sf object but whatever
  spatial_locs <- function(locs) {
    sf_locs <- SpatialPointsDataFrame(data = locs, coords = locs[,c("x", "y")], proj4string = CRS(sa_proj))
    # sf_locs <- st_as_sf(locs, coords = c("x", "y"), crs = sa_proj)
    return(sf_locs)
  }
  #'  Run used & available locations through function
  used_locs <- lapply(spp_all_tracks, spatial_locs)
  #'  3rd-order Selection Covariates
  md_locs <- lapply(md_available, spatial_locs)
  elk_locs <- lapply(elk_available, spatial_locs)
  wtd_locs <- lapply(wtd_available, spatial_locs)
  coug_locs <- lapply(coug_available, spatial_locs)
  wolf_locs <- lapply(wolf_available, spatial_locs)
  

  #'  COVARIATE EXTRACTION & CALCULATIONS  
  #'  ===========================================
  #'  Takes forever but running in parallel helps 
  
  #'  Monitor time
  start.time <- Sys.time()
  
  #'  Identify how many cores I want to use
  detectCores(logical = FALSE)
  cl <- parallel::makeCluster(4) 
  #'  Run in parallel on local computer with specified number of cores
  plan(cluster, workers = cl)
  
  #'  Master function to extract and manipulate covaraite data for each species
  cov_extract <- function(locs) {
    #'  Reproject to WGS84 to match unprojected rasters
    locs_wgs84 <- spTransform(locs, wgs84)
    locs_dist <- spTransform(locs, dist_proj)
    
    #'  1. Extract data from terrain & anthropogenic rasters
    #'  ----------------------------------------------------
    terrain <- raster::extract(terra_stack, locs_wgs84, df = TRUE) # WGS84
    modified <- raster::extract(HM, locs_wgs84, df = TRUE) # WGS84
    road_den <- raster::extract(rdden, locs, df = TRUE) # NAD83(HARN)
    water <- raster::extract(h2o, locs, df = TRUE) # NAD83(HARN)
    #'  Merge into a single data frame of covariates
    join_covs <- terrain %>%
      full_join(modified, by = "ID") %>%
      full_join(road_den, by = "ID") %>%
      full_join(water, by = "ID") %>%
      transmute(
        obs = ID,
        Elev = round(WPPP_DEM_30m, digits = 2),
        Slope = round(WPPP_slope_aspect, digits = 2),
        TPI = round(WPPP_TPI, digits = 2),
        HumanMod = round(WPPP_gHM, digits = 2),
        RoadDen = round(RoadDensity_1km, digits = 2),
        Dist2Water = round(WPPP_dist_to_water_2855, digits = 2)
      ) %>%
      #'  Need to change NAs to 0 for road density (if NA it means there are no
      #'  roads within that 1km pixel and raster pixel was empty)
      mutate(
        RoadDen = ifelse(is.na(RoadDen), 0, RoadDen),
      )
    #'  Make location and animal ID information non-spatial
    #'  Be sure to remove geometry if this is an sf object
    animal <- as.data.frame(locs) 
    #'  Merge animal/time information with covariates
    covs <- as.data.frame(cbind(animal, join_covs))

    #'  2. Extract annual canopy cover data
    #'  -----------------------------------
    animal$obs <- as.numeric(seq(1:nrow(animal)))
    #'  Extract annual canopy cover
    cover18 <- raster::extract(canopy18, locs_wgs84, df = TRUE) %>% # WGS84
      transmute(
        obs = ID,
        Canopy = round(treecov_2018, digits = 2)
      ) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer18" | Season == "Winter1819")
    cover19 <- raster::extract(canopy19, locs_wgs84, df = TRUE) %>%
      transmute(
        obs = ID,
        Canopy = round(treecov_2019, digits = 2)
      ) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer19" | Season == "Winter1920")
    cover20 <- raster::extract(canopy20, locs_wgs84, df = TRUE) %>%
      transmute(
        obs = ID,
        Canopy = round(treecov_2020, digits = 2)
      ) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer20" | Season == "Winter2021")
    canopy1820 <- rbind(cover18, cover19, cover20)
    
    #'  3. Extract Distance to Edge data, based on Cascadia landcover rasters
    #'  reclassified to Forest (forest + woody wetland) and non-Forest. Extracts
    #'  distance to edge of forest and edge of open habitat separately, then 
    #'  takes the minimum distance between the two for an overall Distance to 
    #'  Edge value, regardless of landcover type of focal pixel.
    #'  ------------------------------------------------------------------------
    animal$obs <- as.numeric(seq(1:nrow(animal)))
    #'  Extract distance to edge data at each location
    Dist2Edge18 <- raster::extract(dist2edge_stack18, locs_dist, df = TRUE) %>% # Dist2Edge projection
      transmute(
        obs = ID,
        Dist2Forest = round(Dist2ForestEdge18, 2),
        Dist2Open = round(Dist2OpenEdge18, 2)
      ) %>%
      #'  Combine distance to forest & open edges so one distance to edge value
      #'  for each pixel (essentially drops the NAs and merges into single column)
      rowwise() %>%
      mutate(Dist2Edge = min(Dist2Forest, Dist2Open, na.rm=TRUE)) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer18" | Season == "Winter1819")
    Dist2Edge19 <- raster::extract(dist2edge_stack19, locs_dist, df = TRUE) %>%
      transmute(
        obs = ID,
        Dist2Forest = round(Dist2ForestEdge19, 2),
        Dist2Open = round(Dist2OpenEdge19, 2)
      ) %>%
      rowwise() %>%
      mutate(Dist2Edge = min(Dist2Forest, Dist2Open, na.rm=TRUE)) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer19" | Season == "Winter1920")
    Dist2Edge20 <- raster::extract(dist2edge_stack20, locs_dist, df = TRUE) %>%
      transmute(
        obs = ID,
        Dist2Forest = round(Dist2ForestEdge20, 2),
        Dist2Open = round(Dist2OpenEdge20, 2)
      ) %>%
      rowwise() %>%
      mutate(Dist2Edge = min(Dist2Forest, Dist2Open, na.rm=TRUE)) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer20" | Season == "Winter2021")
    Dist2Edge <- rbind(Dist2Edge18, Dist2Edge19, Dist2Edge20)
    
    #'  4. Extract landcover data from Cascadia landcover rasters. 
    #'  ------------------------------------------------------------------------
    animal$obs <- as.numeric(seq(1:nrow(animal)))
    #'  Extract landcover type at each location (categorical)
    landcover18 <- raster::extract(landcov18, locs_wgs84, df = TRUE) %>% # WGS84
      transmute(
        obs = ID,
        landcov = landcover_2018
      ) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer18" | Season == "Winter1819")
    landcover19 <- raster::extract(landcov19, locs_wgs84, df = TRUE) %>%
      transmute(
        obs = ID,
        landcov = landcover_2019
      ) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer19" | Season == "Winter1920")
    landcover20 <- raster::extract(landcov20, locs_wgs84, df = TRUE) %>%
      transmute(
        obs = ID,
        landcov = landcover_2020
      ) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer20" | Season == "Winter2021")
    landcover <- rbind(landcover18, landcover19, landcover20) %>%
      mutate(
        landcover_type = ifelse(landcov == 101, "Water", landcov),
        landcover_type = ifelse(landcov == 111, "Glacier", landcover_type),
        landcover_type = ifelse(landcov == 121, "Barren", landcover_type),
        landcover_type = ifelse(landcov == 201, "Wetland", landcover_type),
        landcover_type = ifelse(landcov == 201, "Wetland", landcover_type),
        landcover_type = ifelse(landcov == 202, "Woody Wetland", landcover_type),
        landcover_type = ifelse(landcov == 211, "Mesic Grass", landcover_type),
        landcover_type = ifelse(landcov == 212, "Xeric Grass", landcover_type),
        landcover_type = ifelse(landcov == 221, "Mesic Shrub", landcover_type),
        landcover_type = ifelse(landcov == 222, "Xeric Shrub", landcover_type),
        landcover_type = ifelse(landcov == 230, "Forest", landcover_type),
        landcover_type = ifelse(landcov == 301, "Agriculture", landcover_type),
        landcover_type = ifelse(landcov == 331, "Commercial", landcover_type), 
        landcover_type = ifelse(landcov == 332, "Developed", landcover_type)  
      )
    
    #'  5. Join all covatiates together & clean up for inclusion in HMM
    telem_covs <- covs %>%
      full_join(canopy1820, by = c("obs", "ID", "Season")) %>%
      full_join(Dist2Edge, by = c("obs", "ID", "Season")) %>%
      full_join(landcover, by = c("obs", "ID", "Season")) %>%
      transmute(
        ID = ID,
        Season = Season,
        #'  Create a year column based on the season the data were collected
        Year = ifelse(Season == "Summer18" | Season == "Winter1819", "Year1", Season),
        Year = ifelse(Season == "Summer19" | Season == "Winter1920", "Year2", Year),
        Year = ifelse(Season == "Summer20" | Season == "Winter2021", "Year3", Year),
        Elev = Elev,
        Slope = Slope,
        RoadDen = RoadDen,
        Dist2Water = Dist2Water,
        HumanMod = HumanMod,
        CanopyCover = Canopy,
        Dist2Edge = Dist2Edge,
        Landcover = landcov,
        Landcover_type = landcover_type,
        obs = obs) %>%
      mutate(
        #'  Create a study area column using the species ID info
        Area = ifelse(grepl("NE", ID), "NE", "OK"),
        Area = ifelse(grepl("MD", ID), "OK", Area),
        Area = ifelse(grepl("EA", ID), "NE", Area),
        Area = ifelse(grepl("ELK", ID), "NE", Area),
        Area = ifelse(grepl("WTD", ID), "NE", Area),
        #'  Indicate whether this location was used = 1 or available = 0
        Used = 0,
        #'  Add weights to used/available locations (used = 1, available = 5000 
        #'  per Fieberg et al. 2021)
        w = 5000
        )
    
    return(telem_covs)
    
  }
  #'  Run list of species used & available location data through function in parallel
  #'  This will take AWHILE even in parallel
  #'  Keep track of which order of selection is being extracted here!!!
  used_covs <- future_lapply(used_locs, cov_extract, future.seed = TRUE)
  md_avail_covs <- future_lapply(md_locs, cov_extract, future.seed = TRUE)
  elk_avail_covs <- future_lapply(elk_locs, cov_extract, future.seed = TRUE)
  wtd_avail_covs <- future_lapply(wtd_locs, cov_extract, future.seed = TRUE)
  coug_avail_covs <- future_lapply(coug_locs, cov_extract, future.seed = TRUE)
  wolf_avail_covs <- future_lapply(wolf_locs, cov_extract, future.seed = TRUE)
  
  #'  End time keeping
  end.time <- Sys.time()
  #'  Stop running in parallel
  parallel::stopCluster(cl)
  #'  How long did this take?
  difftime(end.time, start.time, units = "hours")
  
  #'  Make sure the used covariates and their weight are marked 1
  for(i in 1:length(used_covs)){
    used_covs[[i]]$Used <- 1
    used_covs[[i]]$w <- 1
  }
  
  #'  Correct study area for wolf data
  #'  No easy way of doing this because ID not associated with study areas
  #'  Double check lists 9 & 10 are wolf summer & winter data
  used_covs[[9]] <- mutate(used_covs[[9]],
                           Area = ifelse(grepl("W61M", ID), "OK", "NE"),  
                           Area = ifelse(grepl("W88M", ID), "OK", Area),
                           Area = ifelse(grepl("W93M", ID), "OK", Area),
                           Area = ifelse(grepl("W94M", ID), "OK", Area),
                           Area = ifelse(grepl("W110M", ID), "OK", Area))
  used_covs[[10]] <- mutate(used_covs[[10]],
                           Area = ifelse(grepl("W61M", ID), "OK", "NE"),  
                           Area = ifelse(grepl("W88M", ID), "OK", Area),
                           Area = ifelse(grepl("W93M", ID), "OK", Area),
                           Area = ifelse(grepl("W94M", ID), "OK", Area),
                           Area = ifelse(grepl("W110M", ID), "OK", Area))
  wolf_avail_covs[[1]] <- mutate(wolf_avail_covs[[1]],
                                Area = ifelse(grepl("W61M", ID), "OK", "NE"),  
                                Area = ifelse(grepl("W88M", ID), "OK", Area),
                                Area = ifelse(grepl("W93M", ID), "OK", Area),
                                Area = ifelse(grepl("W94M", ID), "OK", Area),
                                Area = ifelse(grepl("W110M", ID), "OK", Area))
  wolf_avail_covs[[2]] <- mutate(wolf_avail_covs[[2]],
                                 Area = ifelse(grepl("W61M", ID), "OK", "NE"),
                                 Area = ifelse(grepl("W88M", ID), "OK", Area),
                                 Area = ifelse(grepl("W93M", ID), "OK", Area),
                                 Area = ifelse(grepl("W94M", ID), "OK", Area),
                                 Area = ifelse(grepl("W110M", ID), "OK", Area))
  
  #'  Merge lists across lists of coordinates & covariate data
  combo_data <- function(pts, covs){
    merged <- mapply(data.frame, pts, covs, SIMPLIFY = FALSE)
    return(merged)
  }
  #'  Run used and available locations and covariates through
  used_dat <- combo_data(used_locs, used_covs)
  md_avail_dat <- combo_data(md_available, md_avail_covs)
  elk_avail_dat <- combo_data(elk_available, elk_avail_covs)
  wtd_avail_dat <- combo_data(wtd_available, wtd_avail_covs)
  coug_avail_dat <- combo_data(coug_available, coug_avail_covs)
  wolf_avail_dat <- combo_data(wolf_available, wolf_avail_covs)
  
  #'  Function to drop unneeded columns from list of used data sets
  select_cols <- function(dat) {
    used_skinny <- dat %>%
      dplyr::select(x, y, AnimalID, Season, ID.1, Season.1, Year, Elev, Slope, TPI,
                    RoadDen, Dist2Water, HumanMod, CanopyCover, Dist2Edge, 
                    Landcover, Landcover_type, 
                    obs, Area, Used, w) 
    colnames(used_skinny) <-  c("x", "y", "ID", "Season", "ID.1", "Season.1", 
                                "Year", "Elev", "Slope", "TPI", "RoadDen", 
                                "Dist2Water", "HumanMod", "CanopyCover", "Dist2Edge",
                                "Landcover", "Landcover_type", "obs", "Area", "Used", "w")
    return(used_skinny)
  }
  #'  Run the list of used locations through function
  used_dat <- lapply(used_dat, select_cols)
  
  #'  Merge used & available data per species
  md_dat_all <- rbind(used_dat[[1]], md_avail_dat[[1]], used_dat[[2]], md_avail_dat[[2]])  %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  elk_dat_all <- rbind(used_dat[[3]], elk_avail_dat[[1]], used_dat[[4]], elk_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  wtd_dat_all <- rbind(used_dat[[5]], wtd_avail_dat[[1]], used_dat[[6]], wtd_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  coug_dat_all <- rbind(used_dat[[7]], coug_avail_dat[[1]], used_dat[[8]], coug_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  wolf_dat_all <- rbind(used_dat[[9]], wolf_avail_dat[[1]], used_dat[[10]], wolf_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))

  #'  Save combined data for final RSFs
  save(md_dat_all, file = "./Outputs/RSF_pts/md_dat_all_for_pub.RData")
  save(elk_dat_all, file = "./Outputs/RSF_pts/elk_dat_all_for_pub.RData")
  save(wtd_dat_all, file = "./Outputs/RSF_pts/wtd_dat_all_for_pub.RData")
  save(coug_dat_all, file = "./Outputs/RSF_pts/coug_dat_all_for_pub.RData")
  save(wolf_dat_all, file = "./Outputs/RSF_pts/wolf_dat_all_for_pub.RData")
  
  
  #'  Next step: Run RSFs with Resource_Selection_Function_Model.R
  