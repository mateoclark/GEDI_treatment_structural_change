# Calculates CCDC GEDI interpolation change in treatments

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 08/05/2025

library(dplyr)
library(terra)
library(sf)
library(stringr)

# temporary directory
terraOptions(tempdir = "e:/temp")

# version
version <- "250805"

# years to work with
years <- c("2019", "2020", "2021", "2022")

# Region lookup
regionLookup <- data.frame(
  ccdcName = c("central","sierras","cascades","coast","klamath"),
  US_L3NAME = c("Central California Foothills and Coastal Mountains",
    "Sierra Nevada","Cascades","Coast Range","Klamath Mountains/California High North Coast Range")
)

# Metrics
metrics <- c("agbd","aVDR","cover","fhd_pai_1m_a0","Gini","mPAI_b10","num_modes","pai_a0","rh_98","WI_b10m")

# input folder with annual ccdc images
ccdcDir <- "E:/active/project/calfire_gedi/gedi_ccdc/fusion_v4"

# file geodatabase with treatment & control polygons
fgdb <- "E:/active/project/calfire_gedi/treatments/treatment_arcgis/gedi_treatments.gdb"

# NLCD raster forest and shrub mask (year 2018, year before GEDI analysis)
landcover <- "E:/active/project/calfire_gedi/nlcd/Annual_NLCD_LndCov_2018_CU_C1V0_CA_forest_shrub_mask.tif"

# feature classes
treatment_layer <- "treatment_polygons_250119"
control_layer <- "control_polygons_250117"
regions_layer <- "calfire_map_ecoregions_v3"

# output directory
outDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data"

# output .csv file
#outFile <- paste0("treatment_gedi_ccdc_stats_",version,".csv")
outFile <- paste0("treatment_simulated_gedi_ccdc_stats_",version,".csv")

# treatment selection .csv file
#treatmentSelected <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/master_gedi_treatment_difference_250710.csv"
treatmentSelected <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/master_simulated_gedi_treatment_difference_250801.csv"


##############################################
# Functions
#############################################

# Calculates cell statistics within a polygoon
pixel_statistics <- function(ccdc_raster, polygon){
  extracted_values <- extract(ccdc_raster, polygon)
  statistics <- as.numeric(summary(extracted_values[,2]))[1:6]
  return(statistics)
}

# Function to crop a raster to a polygon, while checking the projections are the same

rasterCrop <- function(raster, polygon) {
  
  # Get the CRS EPSG code of the raster
  epsg_code <- crs(raster, describe = TRUE)$code
  
  # Get the CRS EPSG code of the target polygon
  target_code <- crs(polygon, describe = TRUE)$code
  
  # Initialize a variable to hold the cropped raster
  rasterCrop <- NULL
  
  # If EPSG codes are not the same, reproject to the target
  if (epsg_code != target_code) {
    # Attempt to transform and project the raster
    poly <- tryCatch({
      st_transform(polygon, paste0("EPSG:", epsg_code))
    }, error = function(e) {
      warning("Error transforming polygon to EPSG code:", epsg_code, " - ", e$message)
      return(NULL)
    })
    
    if (is.null(poly)) {
      return(NULL)  # Skip cropping if transformation failed
    }
    
    raster <- tryCatch({
      project(crop(raster, poly), paste0("EPSG:", target_code))
    }, error = function(e) {
      warning("Error cropping and projecting raster - ", e$message)
      return(NULL)
    })
    
    if (is.null(raster)) {
      return(NULL)  # Skip cropping if project failed
    }
  }
  
  # Attempt to crop the final raster to the polygon
  rasterCrop <- tryCatch({
    crop(raster, polygon)
   }, error = function(e) {
    warning("Error cropping raster to polygon - ", e$message)
    return(NULL)
  })
  
  return(rasterCrop)
}

# Function for ccdc statistics for treatments and associated control areas

ccdc_statistics <- function(fgdb, treatments, control, metrics,imageFiles){
  
  # Load land cover mask raster
  lc_mask <- rast(landcover)
  
  # Get unique treatment IDs
  treatmentID <- unique(treatments$TreatmentID)
  
  # List for change data
  treatmentccdcChange <- list()
  
  # Loop through metrics
  for (m in metrics){
      cat(paste0("Metric: ",m,"\n"))
      images<- imageFiles[str_detect(imageFiles, m)]
      
      # Loop through treatments
      for (t in treatmentID){
        
        cat(paste0("Treatment: ",t,"\n"))
        
        # Get the current treatment polygons
        treatment_polygons <- treatments %>% filter(TreatmentID == t) 
        polygon <- treatment_polygons %>% slice(1)
        
        # Extract the start and end years from the treatment polygon; if empty, make it the global min and max, respectively
        StartYear <- polygon$StartYear
        EndYear <- polygon$EndYear
        if (is.na(StartYear)) StartYear = min(treatments$StartYear,na.rm = T)
        if (is.na(EndYear)) EndYear = max(treatments$EndYear,na.rm = T)
        
        if (EndYear < 2023){
          ## Calculate ccdc change statistics and the polygon and pixel level
          pattern <- paste0(polygon$ccdcName[1], ".*", StartYear)
          imageStart <- images[grep(pattern, images)]
          rasterStart <- rast(imageStart)
          pattern <- paste0(polygon$ccdcName[1], ".*", EndYear)
          imageEnd <- images[grep(pattern, images)]
          rasterEnd <- rast(imageEnd)
          
          # ccdc change for treatments using start/end years
          rasterStartCrop <- rasterCrop(rasterStart,polygon)
          rasterEndCrop <- rasterCrop(rasterEnd,polygon)
          
          # Compare extents and resample to match if necessary
          extStart <- ext(rasterStartCrop)
          extEnd <- ext(rasterEndCrop)
          if (!identical(extStart, extEnd)) {
            rasterEndCrop <- resample(rasterEndCrop,rasterStartCrop, method = "near")
          }
          
          # Mask to forest and shrub land cover types
          lc_resamp <- resample(lc_mask,rasterStartCrop, method = "near")
          lc <- crop(lc_resamp,extStart)
          rasterStartCropMasked <- rasterStartCrop * lc
          rasterEndCropMasked <- rasterEndCrop * lc
          
          # Statistics
          statsStart <- pixel_statistics(rasterStartCropMasked,polygon)
          statsEnd <- pixel_statistics(rasterEndCropMasked,polygon)
          TreatmentChangePoly <- statsEnd - statsStart
          rasterChange <- rasterEndCrop - rasterStartCrop
          TreatmentChangePixel <- pixel_statistics(rasterChange,polygon)
          
          # Get control associated with treatment
          cpolygon <- control %>% filter(UniqueID == polygon$ControlID)
          
          # If there is a control
          if (dim(cpolygon)[1] > 0){
            # Get associated raster (may be different region)
            pattern <- paste0(cpolygon$ccdcName[1], ".*", StartYear)
            imageStart <- images[grep(pattern, images)]
            rasterStart <- rast(imageStart)
            pattern <- paste0(cpolygon$ccdcName[1], ".*", EndYear)
            imageEnd <- images[grep(pattern, images)]
            rasterEnd <- rast(imageEnd)
            
            # ccdc change for control using start/end years
            rasterStartCrop <- rasterCrop(rasterStart,cpolygon)
            rasterEndCrop <- rasterCrop(rasterEnd,cpolygon)
            
            # Compare extents and resample to match if necessary
            extStart <- ext(rasterStartCrop)
            extEnd <- ext(rasterEndCrop)
            if (!identical(extStart, extEnd)) {
              rasterEndCrop <- resample(rasterEndCrop,rasterStartCrop, method = "near")
            } 
            
            # Mask to forest and shrub land cover types
            lc_resamp <- resample(lc_mask,rasterStartCrop, method = "near")
            lc <- crop(lc_resamp,extStart)
            rasterStartCropMasked <- rasterStartCrop * lc
            rasterEndCropMasked <- rasterEndCrop * lc
            
            # Statistics
            statsStart <- pixel_statistics(rasterStartCropMasked,cpolygon)
            statsEnd <- pixel_statistics(rasterEndCropMasked,cpolygon)
            ControlChangePoly <- statsEnd - statsStart
            rasterChange <- rasterEndCrop - rasterStartCrop
            ControlChangePixel <- pixel_statistics(rasterChange,cpolygon)
            
          } else { # case with no control polygon
            ControlChangePoly <- rep(NA, 6)
            ControlChangePixel <- rep(NA, 6)
          }

          
          # Add to list
          treatmentccdcChange[[length(treatmentccdcChange) + 1]] <- c(m,t,TreatmentChangePoly,TreatmentChangePixel,
                                                                      ControlChangePoly,ControlChangePixel)
          
          remove(rasterEnd)
          remove(rasterStart)
          remove(TreatmentChangePoly)
          remove(TreatmentChangePixel)
          remove(ControlChangePoly)
          remove(ControlChangePixel)
          remove(polygon)
          remove(cpolygon)
          
        }
      } # treatment end
  } # metrics end
  
  cnames <- c("Metric","TreatmentID",
              "TPolyccdcmin","TPolyccdc1stQ","TPolyccdcmedian","TPolyccdcmean","TPolyccdc3rdQ","TPolyccdcmax",
              "TPixelccdcmin","TPixelccdc1stQ","TPixelccdcmedian","TPixelccdcmean","TPixelccdc3rdQ","TPixelccdcmax",
              "CPolyccdcmin","CPolyccdc1stQ","CPolyccdcmedian","CPolyccdcmean","CPolyccdc3rdQ","CPolyccdcmax",
              "CPixelccdcmin","CPixelccdc1stQ","CPixelccdcmedian","CPixelccdcmean","CPixelccdc3rdQ","CPixelccdcmax")
  
  treatmentStats <- data.frame()
  for (i in 1:length(treatmentccdcChange)){
    df <- as.data.frame(t(treatmentccdcChange[[i]]))
    df[, 4:ncol(df)] <- lapply(df[, 4:ncol(df)], as.numeric)
    if (dim(df)[2] == length(cnames)){
      colnames(df) <- cnames
      treatmentStats <- rbind(treatmentStats,df)
    } else {
      print(paste0("Skipped:",df$V1," ",df$V2))
    }
  }

  return(treatmentStats)
}

#################################
## Main processing
#################################

# Load treatments
treatments <- st_read(fgdb,layer = treatment_layer) %>% rename(TreatmentID = UniqueID)

# Load controls
control <- st_read(fgdb,layer = control_layer)

# Load treatment selection file
treatmentSelect <- read.csv(treatmentSelected) %>% group_by(TreatmentID) %>% slice(1) %>% select(c(TreatmentID,ControlID))

# Get unique treatment IDs
treatmentID <- unique(treatmentSelect$TreatmentID)

# Select treatments
treatments <- treatments %>% filter(TreatmentID %in% treatmentID)
treatments <- left_join(treatments,treatmentSelect,by="TreatmentID") 

# Load regions
regions <- st_read(fgdb,layer = regions_layer)

# Spatial join
treatments_join <- st_join(treatments,regions)
treatments_join <- left_join(treatments_join, regionLookup, by = "US_L3NAME")
control_join <- st_join(control,regions)
control_join <- left_join(control_join, regionLookup, by = "US_L3NAME")

# regions
regionsSelected <- unique(treatments_join$ccdcName)
regionsSelected <- regionsSelected[!regionsSelected == "x"]

statistics <- data.frame()

# Get list of .tif files
imageFiles <- list.files(ccdcDir, pattern = "ccdc", full.names = TRUE, recursive = T)

# Filter images to match years
pattern <- paste(years, collapse = "|")
imageFiles <- imageFiles[grepl(pattern, imageFiles)]
imageFiles <- imageFiles[!grepl("rf_sd",imageFiles)]

for (r in regionsSelected){
  
  cat(r)
  
  # Select treatments in region
  treatmentRegion <- treatments_join %>% filter(ccdcName == r)
  
  # Process ccdc statistics 
  treatmentStats <- ccdc_statistics(fgdb,treatmentRegion,control_join,metrics,imageFiles)
  treatmentStats$region <- r
  
  statistics <- rbind(statistics,treatmentStats)
  
} 

# Combine data and write to .csv file
out <- paste0(outDir,"/", outFile)
write.csv(statistics,out,row.names=F)
