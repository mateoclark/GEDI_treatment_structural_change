# Loads in GEDI structure data and treatment and control polygons, overlays the shots with polygon attributes,
# and outputs treatment and control .csv files with information needed for analysis.

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 07/10/2025

#library(tidyverse)
library(lubridate)
library(sf)
library(dplyr)
library(terra)

##### Parameters #####
# Version of processing
version <- "250710"

# Selected metrics 
metricSelect = c("strct_FHD","strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                 "strct_prop_int_below_10m","strct_RH_10","strct_RH_25",
                 "strct_RH_50","strct_RH_75","strct_RH_98","strct_PAI_inc_00_05","strct_PAI_inc_05_10")

# GEDI shot directory and files
inDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/GEDI_shots"

fileName1 <- paste0(inDir, "/control_zones_buffer500m_z10n_shots.csv")
fileName2 <- paste0(inDir, "/control_zones_buffer500m_z11n_shots.csv")
fileName3 <- paste0(inDir, "/Treatments_UTMz10_Only_08-18-24_shots.csv")
fileName4 <- paste0(inDir, "/Treatments_UTMz11_Only_08-18-24_shots.csv")
fileName5 <- paste0(inDir, "/FACT_FuelsHazard_TimberHarvest_AfterJuly2019_UTMz10n_shots.csv")
fileName6 <- paste0(inDir, "/FACT_FuelsHazard_TimberHarvest_AfterJuly2019_UTMz11n_shots.csv")
fileName7 <- paste0(inDir, "/FACT_Silviculture_AfterJuly2019_UTMz10n_shots.csv")
fileName8 <- paste0(inDir, "/FACT_Silviculture_AfterJuly2019_UTMz10n_shots.csv")
fileName9 <- paste0(inDir, "/c10_extra_woodlands.csv")

# File geodatabase with treatments
fgdb <- "E:/active/project/calfire_gedi/treatments/treatment_arcgis/gedi_treatments.gdb"
#fgdb <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/arcgis/gedi_treatments.gdb"

# NLCD raster (year 2018, year before GEDI analysis)
landcover <- "E:/active/project/calfire_gedi/nlcd/Annual_NLCD_LndCov_2018_CU_C1V0_CA.tif"

# distance to buffer footprints (12.5 m without geolocational error, but 22.5 m with 10 m error)
bufdist <- 22.5

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
  treatments <- st_transform(treatments, crs = 3310) # project to CA Teale if needed
  summary <- summarizePctCover(landcover,treatments)
  treatments$poly_lc_maj1 <- summary$majority_class
  treatments$poly_lc_pct1 <- summary$majority_pct
  treatments$poly_lc_maj2 <- summary$second_class
  treatments$poly_lc_pct2 <- summary$second_pct
  treatments$poly_lc_maj3 <- summary$third_class
  treatments$poly_lc_pct3 <- summary$third_pct
  
  dfFilter <- df %>% filter(Type == "Treatment")
  shots <- st_as_sf(dfFilter, coords = c("x_teale", "y_teale"),remove=F) 
  shots <- st_set_crs(shots, 3310)
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
prepare_control <- function(df,control_layer,landcover,bufdist){
  
  control <- st_read(fgdb,layer = control_layer) 
  control <- st_transform(control, crs = 3310) # project to CA Teale if needed
  summary <- summarizePctCover(landcover,control)
  control$poly_lc_maj1 <- summary$majority_class
  control$poly_lc_pct1 <- summary$majority_pct
  control$poly_lc_maj2 <- summary$second_class
  control$poly_lc_pct2 <- summary$second_pct
  control$poly_lc_maj3 <- summary$third_class
  control$poly_lc_pct3 <- summary$third_pct
  
  dfFilter <- df %>% filter(Type == "Control") 
  shots <- st_as_sf(dfFilter, coords = c("x_teale", "y_teale"),remove=F) 
  shots <- st_set_crs(shots, 3310) 
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

# Function to create a CRS string for UTM based on the zone
crsUTM <- function(zone) {
  paste0("+proj=utm +zone=", zone, " +datum=WGS84 +units=m +no_defs")
}

############################
##### Begin processing #####
############################

# Load treatment GEDI data and select columns
gediInfo <- c("GEDI_shot","GEDI_date","GEDI_utmX","GEDI_utmY","GEDI_utmZ","GEDI_beam","GEDI_elev")
columnSelect <- c(gediInfo,metricSelect)
df1 <- read.csv(fileName1) %>% select(!!columnSelect) %>% mutate(Type = "Control") 
df2 <- read.csv(fileName2) %>% select(!!columnSelect) %>% mutate(Type = "Control") 
df3 <- read.csv(fileName3) %>% select(!!columnSelect) %>% mutate(Type = "Treatment") 
df4 <- read.csv(fileName4) %>% select(!!columnSelect) %>% mutate(Type = "Treatment") 
df5 <- read.csv(fileName5) %>% select(!!columnSelect) %>% mutate(Type = "Treatment") 
df6 <- read.csv(fileName6) %>% select(!!columnSelect) %>% mutate(Type = "Treatment") 
df7 <- read.csv(fileName7) %>% select(!!columnSelect) %>% mutate(Type = "Treatment") 
df8 <- read.csv(fileName8) %>% select(!!columnSelect) %>% mutate(Type = "Treatment")
df9 <- read.csv(fileName9) %>% select(!!columnSelect) %>% mutate(Type = "Control")  # extra woodland controls in Central California

# Combine data frames
dfCombined <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9)

# Create a mean PAI 0-10m metric and remove PAI increments used to calculated it
dfCombined <- dfCombined %>% rowwise() %>% mutate(strct_mPAI_b10 = mean(c(strct_PAI_inc_00_05,strct_PAI_inc_05_10))) %>%
              select(-c("strct_PAI_inc_00_05","strct_PAI_inc_05_10"))

# Get month and year, accounting for any discrepancies in date format in csvs
dfCombined$Date <-as.Date(dfCombined$GEDI_date)
dfCombined$Month <- month(as.POSIXlt(dfCombined$Date, format="%Y/%m/%d",tz="UTC"))
dfCombined$Year <- year(as.POSIXlt(dfCombined$Date, format="%Y/%m/%d",tz="UTC"))
dfCombined <- dfCombined  %>% select(-GEDI_date) %>% filter(!is.na(Date))

# Project to California Teale Albers
# Iterate through each UTM zone and transform to California Albers
# Need a uniform projection for distance matrix
dfCombinedTeale <- do.call(rbind, lapply(split(dfCombined, dfCombined$GEDI_utmZ), function(x) {
  st_as_sf(x, coords = c("GEDI_utmX","GEDI_utmY"), crs = crsUTM(x$GEDI_utmZ[1]), remove = F) %>%
    st_transform(crs = 3310)  # California Albers (Teale Albers) NAD83 EPSG:3310
}))

# Add converted coordinates and convert back to data frame
dfCombinedTeale <- as.data.frame(dfCombinedTeale) %>% 
  mutate(x_teale = st_coordinates(geometry)[,1],
         y_teale = st_coordinates(geometry)[,2]) %>%
  select(-geometry)


# Process treatments #
treatment_layer <- "treatment_polygons_250119"
dfTreatments <- prepare_treatment(dfCombinedTeale,treatment_layer,landcover,bufdist)

# Process controls
control_layer <- "control_polygons_250117"
dfControl <- prepare_control(dfCombinedTeale,control_layer,landcover,bufdist)

# Bind data frames
dfTreatment <- distinct(dfTreatments,.keep_all= TRUE)
dfControl <- distinct(dfControl,.keep_all= TRUE)

# Remove unwanted columns
dfTreatment1 <- dfTreatment %>% select(-c(geometry.x,geometry.y,geometry))
dfControl1 <- dfControl %>% select(-c(ORIG_FID,geometry.x,geometry.y,geometry))

# Label shots as pre- or post-treatment
dfTreatment1 <- dfTreatment1 %>%
  mutate(Timing = case_when(
    Year < StartYear | (Year == StartYear & Month <= 7) ~ "Pre-treatment",
    Year > EndYear | (Year == EndYear & Month >= 7) ~ "Post-treatment",
    is.na(StartYear) ~ "NA",
    TRUE ~ "During treatment" 
  ))

# Label treatment shot timing relative to fire
dfTreatment1 <- dfTreatment1 %>%
  mutate(ShotFireTiming = case_when(
    is.na(FireStartDate) ~ "No fire",
    Timing == "During treatment" ~ "During treatment",
    FireTiming == "Fire before treatment" & (Date > FireEndDate|is.na(FireEndDate)) ~ "Before treatment and after past fire",
    FireTiming == "Fire before treatment" & Date <= FireEndDate ~ "Before treatment and before/during past fire",
    FireTiming == "Fire after treatment" & Date < FireStartDate ~ "After treatment and before post-treatment fire",
    FireTiming == "Fire after treatment" & Date >= FireStartDate ~ "After treatment and after/during post-treatment fire",
    TRUE ~ "Other" 
  ))

# Make a final dataframe
dfTreatmentFinal <- as.data.frame(dfTreatment1) %>% filter(!is.na(UniqueID))
dfControlFinal <- as.data.frame(dfControl1) %>% filter(!is.na(UniqueID))

# Write out GEDI shots .csv file
write.csv(dfTreatmentFinal, paste0("G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/treatment_gedi_data_",version,".csv"), row.names = FALSE)
write.csv(dfControlFinal, paste0("G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/control_gedi_data_",version,".csv"), row.names = FALSE)

