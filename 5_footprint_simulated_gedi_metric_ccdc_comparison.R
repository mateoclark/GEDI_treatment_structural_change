# Compares simulated GEDI to CCDC interpolated metrics

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 01/06/26

library(dplyr)
library(readr)
library(terra)
library(sf)
library(stringr)

# temporary directory
terraOptions(tempdir = "e:/temp")

    
#################################
########### Parameters ##########
#################################

# version of processing
version <- "260106"

# input directory
inDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics"

outFile <- paste0(inDir,"/data/simulated_gedi_ccdc_fusion_",version,".csv")

# input .csv file with GEDI and treatment and control information
inputTreatmentFile <- paste0(inDir,"/data/treatment_simulated_gedi_data_filtered_",version,".csv")
inputControlFile <- paste0(inDir,"/data/control_simulated_gedi_data_filtered_",version,".csv")

# Simulated GEDI lookup 
metrics = c("strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98","strct_mPAI_b10","strct_FHD","strct_cover","strct_VDR","strct_tAGBD","strct_tPAI","strct_elev",
                    "gauss_RH_25","gauss_RH_50","gauss_RH_75","gauss_RH_98","gauss_mPAI_b10","gauss_FHD","gauss_cover","gauss_VDR","gauss_tAGBD","gauss_tPAI","gauss_elev")

# CCDC years to work with
years <- c("2019", "2020", "2021", "2022")

# Region lookup
regionLookup <- data.frame(
  ccdcName = c("central","sierras","cascades","coast","klamath"),
  US_L3NAME = c("Central California Foothills and Coastal Mountains",
                "Sierra Nevada","Cascades","Coast Range","Klamath Mountains/California High North Coast Range")
)

# CCDC Metrics
metricsCCDC <- data.frame(
  metrics_raster = c("agbd","aVDR","cover","mPAI_b10","rh_98","fhd_pai_1m","pai_a0"),
  metrics_df = c("ccdc_strct_tAGBD", "ccdc_strct_VDR","ccdc_strct_cover","ccdc_strct_mPAI_b10","ccdc_strct_RH_98","ccdc_strct_FHD","ccdc_strct_tPAI")
)

# input folder with annual ccdc images
ccdcDir <- "E:/active/project/calfire_gedi/gedi_ccdc/fusion_v4"

# file geodatabase with treatment & control polygons
fgdb <- "E:/active/project/calfire_gedi/treatments/treatment_arcgis/gedi_treatments.gdb"

# feature classes
regions_layer <- "calfire_map_ecoregions_v3"

################################################
##### Load and filter simulated GEDI data ######
################################################

# load input GEDI data
df1 <- read_csv(inputTreatmentFile)

# get control footprints
df2 <- read_csv(inputControlFile) 

# Rename teale columns to fit functions
df1 <- df1 %>% select(any_of(metrics),ID,YearALS, Type, UniqueID, x_teale, y_teale, landcover_pct, als_start, als_end, ground_slope, beamDense, 
                      shot_lc_maj1,shot_lc_pct1,shot_lc_maj2,shot_lc_pct2,shot_lc_maj3,shot_lc_pct3)
df2 <- df2 %>% select(any_of(metrics),ID,YearALS, Type, UniqueID, x_teale, y_teale, landcover_pct, als_start, als_end, ground_slope, beamDense,
                      shot_lc_maj1,shot_lc_pct1,shot_lc_maj2,shot_lc_pct2,shot_lc_maj3,shot_lc_pct3)

df_combined <- rbind(df1,df2)

month_to_season <- function(m) {
  dplyr::case_when(
    m %in% c(12, 1, 2) ~ "winter",
    m %in% c(3, 4, 5)  ~ "spring",
    m %in% c(6, 7, 8)  ~ "summer",
    m %in% c(9,10,11)  ~ "fall",
    TRUE ~ NA_character_
  )
}

span_season_label <- function(start_date, end_date) {
  if (is.na(start_date) || is.na(end_date)) return(NA_character_)
  if (end_date < start_date) return(NA_character_)  # or stop("end < start")
  
  # get all months touched by the interval (inclusive)
  mseq <- seq.Date(floor_date(start_date, "month"),
                   floor_date(end_date,   "month"),
                   by = "month")
  
  seasons <- month_to_season(month(mseq))
  seasons <- seasons[!is.na(seasons)]
  
  # keep unique seasons in the order they appear
  seasons <- seasons[!duplicated(seasons)]
  
  paste(seasons, collapse = "/")
}

df_combined <- df_combined %>%
  mutate(
    season = mapply(span_season_label, als_start, als_end)
  )

################################################
##### Sample CCDC GEDI interpolated data  ######
################################################

sf_df <- st_as_sf(
  df_combined,
  coords = c("x_teale", "y_teale"),  # specify the x and y columns
  crs    = 6414                     # set coordinate reference system to EPSG:6414
)
sf_df <- st_transform(sf_df,crs = "EPSG:3310")

# Load regions
regions <- st_read(fgdb,layer = regions_layer)
regions <- st_transform(regions,crs = "EPSG:3310")
df_join <- st_join(sf_df,regions)
df_join <- left_join(df_join, regionLookup, by = "US_L3NAME")

# regions
regionsSelected <- unique(df_join$ccdcName)
regionsSelected <- regionsSelected[!regionsSelected == "x"]

statistics <- data.frame()

# Get list of .tif files
imageFiles <- list.files(ccdcDir, pattern = "ccdc.*\\.tif$", full.names = TRUE, recursive = T)
imageFiles <- imageFiles[ ! grepl("_sd_", imageFiles) ]

# Filter images to match years
pattern <- paste(years, collapse = "|")
imageFiles <- imageFiles[grepl(pattern, imageFiles)]

# Filter df_join to CCDC years, add ID
df_join <- df_join %>% filter(YearALS %in% years)
df_join$row_id <- 1:nrow(df_join)

# Loop through rasters and extract values
dr <- NULL
for (r in regionsSelected){
  cat(paste0("Region: ",r,"\n"))
  dfRegion <- df_join %>% filter(ccdcName == r)
  years_als <- intersect(unique(dfRegion$YearALS), years)
  for (y in years_als){
    cat(paste0("Year: ",y,"\n"))
    values <- NULL
    dfYear <- dfRegion %>% filter(YearALS == y)
    for (m in metricsCCDC$metrics_raster){
      cat(paste0("Metric: ",m,"\n"))
      images <- imageFiles[str_detect(imageFiles, m)]
      pattern <- paste0(r, ".*", y)
      imageSelect <- images[grep(pattern, images)]
      raster <- rast(imageSelect)
      vals <- terra::extract(raster, dfYear)
      values <- cbind(values,vals[,2])
    }
    colnames(values) <- metricsCCDC$metrics_df
    dv <- data.frame(values)
    dv$row_id <- dfYear$row_id
    dr <- rbind(dr,dv)
  }
}

# Add Teale coordinate fields to data
coords_teale <- st_coordinates(df_join)
df_join$x_teale <- coords_teale[,1]
df_join$y_teale <- coords_teale[,2]

# Join extracted values back to simulated data and scale CCDC data
df_final <- left_join(as.data.frame(df_join),dr,by="row_id") %>% select(-c(geometry,Shape_Leng, Shape_Length, Shape_Area))
df_final <- df_final %>% mutate(
  ccdc_strct_VDR = ccdc_strct_VDR/1000,
  ccdc_strct_cover = ccdc_strct_cover/100,
  ccdc_strct_mPAI_b10 = ccdc_strct_mPAI_b10/1000,
  ccdc_strct_RH_98 = ccdc_strct_RH_98/100,
  ccdc_strct_FHD = ccdc_strct_FHD/1000,
  ccdc_strct_tAGBD = ccdc_strct_tAGBD/10,
  ccdc_strct_tPAI = ccdc_strct_tPAI/1000
)

# Write out .csv file 
write.csv(df_final, outFile, row.names = FALSE)

