# Loads in GEDI structure data and treatment and control polygons, overlays the shots with polygon attributes,
# and outputs treatment and control .csv files with information needed for analysis.
# This version for simulated GEDI footprints using real ALS ground from the UMD fork of the Hancock simulator.
# https://bitbucket.org/StevenHancock/gedisimulator/src/0309c254006bd43f3c51ba41436462316473c299/


# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 01/06/26

#library(tidyverse)
library(lubridate)
library(sf)
library(dplyr)
library(terra)
library(readr)
library(stringr)

##### Parameters #####
# Version of processing
version <- "260106"

# Selected metrics 
metricSelect = c("rhReal_25","rhReal_50","rhReal_70","rhReal_75","rhReal_95","rhReal_98","tLAI0t10mean",
                 "tLAI0t10mean_1m","tLAI_1m","FHD","FHDcanHist","FHDhist","FHD_tLAI_1m","ALS_cover",
                 "rhGauss_25","rhGauss_50","rhGauss_70","rhGauss_75","rhGauss_95","rhGauss_98",
                 "gLAI0t10mean","gLAI0t10mean_1m","FHDcanGhist","FHDcanGauss",
                 "FHD_gLAI_1m","gLAI_1m","cover")

# GEDI shot directory and files
inDir <- "E:/active/project/calfire_gedi/bailey_contract/hancock_simulation"

fileName1 <- paste0(inDir, "/treatment_simulated_waveform_metrics_tpai_fhd_6339.csv")
fileName2 <- paste0(inDir, "/treatment_simulated_waveform_metrics_tpai_fhd_6340.csv")
fileName3 <- paste0(inDir, "/control_simulated_waveform_metrics_tpai_fhd_6339.csv")
fileName4 <- paste0(inDir, "/control_simulated_waveform_metrics_tpai_fhd_6340.csv")

# File geodatabase with treatments
fgdb <- "E:/active/project/calfire_gedi/treatments/treatment_arcgis/gedi_treatments.gdb"
#fgdb <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/arcgis/gedi_treatments.gdb"

# NLCD raster (year 2018, year before GEDI analysis)
landcover <- "E:/active/project/calfire_gedi/nlcd/Annual_NLCD_LndCov_2018_CU_C1V0_CA.tif"

# ALS polygon layer
als_layer <- "E:/active/project/calfire_gedi/bailey_contract/wesm_teale_selected.shp"

# distance to buffer footprints (12.5 m without geolocational error, but 22.5 m with 10 m error)
bufdist <- 12.5

##### Functions #####

# summarize percent cover in polygons
summarizePctCover <- function(landcover,polygon){
  landcover_raster <- rast(landcover)
  extracted_values <- terra::extract(landcover_raster, vect(polygon))
  summary <- extracted_values %>%
    group_by(ID) %>% 
    summarize(
      freq_table = list(table(Band_1)), 
      majority_class = as.numeric(names(which.max(table(Band_1)))),
      majority_pct = max(table(Band_1)) / sum(table(Band_1)) * 100,
      second_class = {
        freq <- table(Band_1)
        if (length(freq) > 1) as.numeric(names(sort(freq, decreasing = TRUE)[2])) else NA
      },
      second_pct = {
        freq <- table(Band_1)
        if (length(freq) > 1) sort(freq, decreasing = TRUE)[2] / sum(freq) * 100 else NA
      },
      third_class = {
        freq <- table(Band_1)
        if (length(freq) > 2) as.numeric(names(sort(freq, decreasing = TRUE)[3])) else NA
      },
      third_pct = {
        freq <- table(Band_1)
        if (length(freq) > 2) sort(freq, decreasing = TRUE)[3] / sum(freq) * 100 else NA
      }
    )
  return(summary)
}

# Function to process treatments
prepare_treatment <- function(df,treatment_layer,landcover,bufdist){
  
  treatments <- st_read(fgdb,layer = treatment_layer)
  treatments <- st_transform(treatments, crs = 6414) # project to CA Teale if needed
  summary <- summarizePctCover(landcover,treatments)
  treatments$poly_lc_maj1 <- summary$majority_class
  treatments$poly_lc_pct1 <- summary$majority_pct
  treatments$poly_lc_maj2 <- summary$second_class
  treatments$poly_lc_pct2 <- summary$second_pct
  treatments$poly_lc_maj3 <- summary$third_class
  treatments$poly_lc_pct3 <- summary$third_pct
  
  shots <- st_as_sf(df, coords = c("Teale_X", "Teale_Y"),remove=F) 
  shots <- st_set_crs(shots, 6414)
  
  dfJoin <- st_join(shots, als_info,left = TRUE)
  dfFilter <- dfJoin %>% filter(YearALS == InfoYearALS)
  shots <- dfFilter[!duplicated(dfFilter$unique_ID), ]
  
  shotsTreatments <- st_join(shots, treatments,left = TRUE)
  shotsTreatments <- shotsTreatments[!is.na(shotsTreatments$StartYear), ]
  shotsTreatments$ShotTreatmentID <- paste0(shotsTreatments$UniqueID,"_", seq(1, by = 1, length.out = nrow(shotsTreatments)))
  
  buffers <- st_buffer(shotsTreatments, dist = bufdist) 
  buffers$Area <- st_area(buffers)
  bufferTreatmentIntersection <- st_intersection(buffers, treatments)
  bufferTreatmentIntersection$intersectionArea <- st_area(bufferTreatmentIntersection)
  shotPct <- bufferTreatmentIntersection %>%
    group_by(ShotTreatmentID) %>% 
    summarize(
      intersectionAreaSum = sum(intersectionArea),
      BufferArea = min(Area),
      .groups = "keep"
    ) %>%
    mutate(PolygonPctCovered = as.numeric(intersectionAreaSum / BufferArea * 100))
  shotsTreatments <- left_join(shotsTreatments,data.frame(shotPct),by="ShotTreatmentID")
  treatmentBoundaries <- st_boundary(treatments)
  closestTreatmentIndices <- st_nearest_feature(buffers, treatmentBoundaries) # Indices of nearest treatment polygons
  distanceToTreatmentEdge <- st_distance(buffers, treatmentBoundaries[closestTreatmentIndices, ], by_element = TRUE)
  buffers <- buffers %>%
    mutate(distanceToTreatmentEdge = as.numeric(distanceToTreatmentEdge)) %>%
    select(ShotTreatmentID,distanceToTreatmentEdge)
  
  summary <- summarizePctCover(landcover,buffers)
  buffers$shot_lc_maj1 <- summary$majority_class
  buffers$shot_lc_pct1 <- summary$majority_pct
  buffers$shot_lc_maj2 <- summary$second_class
  buffers$shot_lc_pct2 <- summary$second_pct
  buffers$shot_lc_maj3 <- summary$third_class
  buffers$shot_lc_pct3 <- summary$third_pct
  
  shotsTreatments <- left_join(shotsTreatments,data.frame(buffers),by="ShotTreatmentID")
  dfTreatments <- data.frame(shotsTreatments) 
  
  return(dfTreatments)
  
}

# Function to process
prepare_control <- function(df,control_layer,landcover,bufdist,als_info){
  
  control <- st_read(fgdb,layer = control_layer) 
  control <- st_transform(control, crs = 6414) # project to CA Teale if needed
  summary <- summarizePctCover(landcover,control)
  control$poly_lc_maj1 <- summary$majority_class
  control$poly_lc_pct1 <- summary$majority_pct
  control$poly_lc_maj2 <- summary$second_class
  control$poly_lc_pct2 <- summary$second_pct
  control$poly_lc_maj3 <- summary$third_class
  control$poly_lc_pct3 <- summary$third_pct
  
  dfFilter <- df %>% filter(Type == "control") 
  shots <- st_as_sf(dfFilter, coords = c("Teale_X", "Teale_Y"),remove=F) 
  shots <- st_set_crs(shots, 6414) 
  
  dfJoin <- st_join(shots, als_info,left = TRUE)
  dfFilter <- dfJoin %>% filter(YearALS == InfoYearALS)
  shots <- dfFilter[!duplicated(dfFilter$unique_ID), ]
  
  shotsControl <- st_join(shots, control,left = TRUE) %>% filter(!is.na(UniqueID))
  shotsControl$ShotControlID <- paste0(shotsControl$UniqueID,"_", seq(1, by = 1, length.out = nrow(shotsControl)))
  buffers <- st_buffer(shotsControl, dist = bufdist) 
  buffers$Area <- st_area(buffers)
  bufferControlIntersection <- st_intersection(buffers, control)
  bufferControlIntersection$intersectionArea <- st_area(bufferControlIntersection)
  shotPct <- bufferControlIntersection %>%
    group_by(ShotControlID) %>% 
    summarize(
      intersectionAreaSum = sum(intersectionArea),
      BufferArea = min(Area),
      .groups = "keep"
    ) %>%
    mutate(PolygonPctCovered = as.numeric(intersectionAreaSum / BufferArea * 100))
  shotsControl <- left_join(shotsControl,data.frame(shotPct),by="ShotControlID")
  controlBoundaries <- st_boundary(control)
  closestControlIndices <- st_nearest_feature(buffers, controlBoundaries) # Indices of nearest treatment polygons
  distanceToControlEdge <- st_distance(buffers, controlBoundaries[closestControlIndices, ], by_element = TRUE)
  buffers <- buffers %>%
    mutate(distanceToControlEdge = as.numeric(distanceToControlEdge)) %>%
    select(ShotControlID,distanceToControlEdge)
  summary <- summarizePctCover(landcover,buffers)
  buffers$shot_lc_maj1 <- summary$majority_class
  buffers$shot_lc_pct1 <- summary$majority_pct
  buffers$shot_lc_maj2 <- summary$second_class
  buffers$shot_lc_pct2 <- summary$second_pct
  buffers$shot_lc_maj3 <- summary$third_class
  buffers$shot_lc_pct3 <- summary$third_pct
  
  shotsControl <- left_join(shotsControl,data.frame(buffers),by="ShotControlID")
  dfControl <- data.frame(shotsControl) 
  
  return(dfControl)
  
}

############################
##### Begin processing #####
############################

# Load treatment GEDI data and select columns
gediInfo <- c("ID","YearALS","YearGEDI","true_ground","ground_slope","unique_ID","Teale_X","Teale_Y",
              "point_density","pointDense","beamDense","file_sizeMB","n_points","type",
              "gHeight")
columnSelect <- c(gediInfo,metricSelect)
df1 <- read_csv(fileName1) %>% select(!!columnSelect) %>% rename(Type = type)
df2 <- read_csv(fileName2) %>% select(!!columnSelect) %>% rename(Type = type)
df3 <- read_csv(fileName3) %>% select(!!columnSelect) %>% rename(Type = type)
df4 <- read_csv(fileName4) %>% select(!!columnSelect) %>% rename(Type = type)

dfT <- rbind(df1,df2)
dfC <- rbind(df3,df4)

# Load ALS information
als_info <- st_read(als_layer) %>% select(workunt, cllct_s,cllct_n) %>%
  rename(als_start = cllct_s, als_end = cllct_n)

als_info$InfoYearALS <- str_extract(als_info$workunt, "(19|20)[0-9]{2}")
als_info$InfoYearALS <- ifelse(
  is.na(als_info$InfoYearALS),
  2022,
  as.integer(als_info$InfoYearALS)
)
als_info <- als_info %>% select(-workunt)
als_info <- st_transform(als_info, crs = 6414) # project to CA Teale if needed

# Process treatments #
treatment_layer <- "treatment_polygons_250119"
dfTreatments <- prepare_treatment(dfT,treatment_layer,landcover,bufdist,als_info)

# Process controls
control_layer <- "control_polygons_250117"
dfControl <- prepare_control(dfC,control_layer,landcover,bufdist,als_info)

# Bind data frames
dfTreatment <- distinct(dfTreatments,.keep_all= TRUE)
dfControl <- distinct(dfControl,.keep_all= TRUE)

# Remove unwanted columns
dfTreatment1 <- dfTreatment %>% select(-c(geometry.x,geometry.y,geometry))
dfControl1 <- dfControl %>% select(-c(ORIG_FID,geometry.x,geometry.y,geometry))

# Label shots as pre- or post-treatment
dfTreatment1 <- dfTreatment1 %>%
  mutate(Timing = case_when(
    YearALS < StartYear ~ "Pre-treatment",
    YearALS > EndYear ~ "Post-treatment",
    is.na(StartYear) ~ "NA",
    TRUE ~ "During treatment" 
  ))

dfTreatment1 <- dfTreatment1 %>%
  mutate(ShotFireTiming = case_when(
    is.na(FireStartDate) ~ "No fire",
    Timing == "During treatment" ~ "During treatment",
    FireTiming == "Fire before treatment" & (YearALS > FireEndDate|is.na(FireEndDate)) ~ "Before treatment and after past fire",
    FireTiming == "Fire before treatment" & YearALS <= FireEndDate ~ "Before treatment and before/during past fire",
    FireTiming == "Fire after treatment" & YearALS < FireStartDate ~ "After treatment and before post-treatment fire",
    FireTiming == "Fire after treatment" & YearALS >= FireStartDate ~ "After treatment and after/during post-treatment fire",
    TRUE ~ "Other" 
  ))


# Make a final dataframe
dfTreatmentFinal <- as.data.frame(dfTreatment1) %>% filter(!is.na(UniqueID))
dfControlFinal <- as.data.frame(dfControl1) %>% filter(!is.na(UniqueID))

# Add aVDR calculation (Clark et al 2025; Hakkenberg et al., 2024)
dfTreatmentFinal <- dfTreatmentFinal %>%
  mutate(aVDR = if_else(rhReal_50 >= 0 & rhReal_98 >= 0,
                        1 - ((rhReal_98 - rhReal_50) / rhReal_98),
                        NA_real_),
         gaVDR = if_else(rhGauss_50 >= 0 & rhGauss_98 >= 0,
                   1 - ((rhGauss_98 - rhGauss_50) / rhGauss_98),
                   NA_real_))

dfControlFinal <- dfControlFinal %>%
  mutate(aVDR = if_else(rhReal_50 >= 0 & rhReal_98 >= 0,
                        1 - ((rhReal_98 - rhReal_50) / rhReal_98),
                        NA_real_),
         gaVDR = if_else(rhGauss_50 >= 0 & rhGauss_98 >= 0,
                         1 - ((rhGauss_98 - rhGauss_50) / rhGauss_98),
                         NA_real_))

# Add AGBD from Kellner et al., 2023
dfTreatmentFinal <- dfTreatmentFinal %>%
  mutate(tAGBD = case_when(
    shot_lc_maj1 == 41 ~ 1.052 * (-120.777 + 5.508 * sqrt(rhReal_50 + 100) + 6.808 * sqrt(rhReal_98 + 100))^2,# North America DBT
    shot_lc_maj1 == 42 ~ 1.013 * (-114.355 + 8.401 * sqrt(rhReal_70 + 100) + 3.346 * sqrt(rhReal_98 + 100))^2, # North America ENT
    shot_lc_maj1 == 43 ~ 1.013 * (-114.355 + 8.401 * sqrt(rhReal_70 + 100) + 3.346 * sqrt(rhReal_98 + 100))^2, # North America ENT (mixed forest but assume more conifer)
    shot_lc_maj1 == 52 ~ 1.118 * (-124.832 + 12.426 * sqrt(rhReal_98 + 100))^2, # GSW - shrubs, grasslands
    TRUE ~ NA
  ),gAGBD = case_when(
    shot_lc_maj1 == 41 ~ 1.052 * (-120.777 + 5.508 * sqrt(rhGauss_50 + 100) + 6.808 * sqrt(rhGauss_98 + 100))^2,# North America DBT
    shot_lc_maj1 == 42 ~ 1.013 * (-114.355 + 8.401 * sqrt(rhGauss_70 + 100) + 3.346 * sqrt(rhGauss_98 + 100))^2, # North America ENT
    shot_lc_maj1 == 43 ~ 1.013 * (-114.355 + 8.401 * sqrt(rhGauss_70 + 100) + 3.346 * sqrt(rhGauss_98 + 100))^2, # North America ENT (mixed forest but assume more conifer)
    shot_lc_maj1 == 52 ~ 1.118 * (-124.832 + 12.426 * sqrt(rhGauss_98 + 100))^2, # GSW - shrubs, grasslands
    TRUE ~ NA
  ))

dfControlFinal <- dfControlFinal %>%
  mutate(tAGBD = case_when(
    shot_lc_maj1 == 41 ~ 1.052 * (-120.777 + 5.508 * sqrt(rhReal_50 + 100) + 6.808 * sqrt(rhReal_98 + 100))^2, # North America DBT
    shot_lc_maj1 == 42 ~ 1.013 * (-114.355 + 8.401 * sqrt(rhReal_70 + 100) + 3.346 * sqrt(rhReal_98 + 100))^2, # North America ENT
    shot_lc_maj1 == 43 ~ 1.013 * (-114.355 + 8.401 * sqrt(rhReal_70 + 100) + 3.346 * sqrt(rhReal_98 + 100))^2, # North America ENT (mixed forest but assume more conifer)
    shot_lc_maj1 == 52 ~ 1.118 * (-124.832 + 12.426 * sqrt(rhReal_98 + 100))^2, # GSW - shrubs, grasslands
    TRUE ~ NA
  ),gAGBD = case_when(
    shot_lc_maj1 == 41 ~ 1.052 * (-120.777 + 5.508 * sqrt(rhGauss_50 + 100) + 6.808 * sqrt(rhGauss_98 + 100))^2,# North America DBT
    shot_lc_maj1 == 42 ~ 1.013 * (-114.355 + 8.401 * sqrt(rhGauss_70 + 100) + 3.346 * sqrt(rhGauss_98 + 100))^2, # North America ENT
    shot_lc_maj1 == 43 ~ 1.013 * (-114.355 + 8.401 * sqrt(rhGauss_70 + 100) + 3.346 * sqrt(rhGauss_98 + 100))^2, # North America ENT (mixed forest but assume more conifer)
    shot_lc_maj1 == 52 ~ 1.118 * (-124.832 + 12.426 * sqrt(rhGauss_98 + 100))^2, # GSW - shrubs, grasslands
    TRUE ~ NA
  ))

# Write out GEDI shots .csv file
write.csv(dfTreatmentFinal, paste0("G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/treatment_simulated_gedi_data_",version,".csv"), row.names = FALSE)
write.csv(dfControlFinal, paste0("G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/control_simulated_gedi_data_",version,".csv"), row.names = FALSE)

