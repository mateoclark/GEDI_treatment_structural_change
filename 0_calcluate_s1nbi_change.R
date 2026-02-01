# Calculates S1 normalized backscatter index change in treatments

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 07/31/2025

library(dplyr)
library(terra)
library(sf)
library(stringr)

# temporary directory
terraOptions(tempdir = "e:/temp")

# version
version <- "250731"

# years to work with
years <- c("2019", "2020", "2021", "2022")

# input folder with annual s1 images
s1Dir <- "E:/active/project/calfire_gedi/treatments/sentinel1"

# file geodatabase with treatment & control polygons
fgdb <- "E:/active/project/calfire_gedi/treatments/treatment_arcgis/gedi_treatments.gdb"

# NLCD raster forest and shrub mask (year 2018, year before GEDI analysis)
landcover <- "E:/active/project/calfire_gedi/nlcd/Annual_NLCD_LndCov_2018_CU_C1V0_CA_forest_shrub_mask.tif"

# output directory
outDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data"

# output .csv file
outFileVV <- paste0("treatment_gedi_s1ndbi_vv_stats_",version,".csv")
outFileVH <- paste0("treatment_gedi_s1ndbi_vh_stats_",version,".csv")

# outFileVV <- paste0("treatment_simulated_gedi_s1ndbi_vv_stats_",version,".csv")
# outFileVH <- paste0("treatment_simulated_gedi_s1ndbi_vh_stats_",version,".csv")

# treatment selection .csv file
treatmentSelected <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/master_gedi_treatment_difference_250710.csv"
#treatmentSelected <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data/master_simulated_gedi_treatment_difference_250722.csv"


##############################################
# Functions
#############################################

# Calculates cell statistics within a polygoon
pixel_statistics <- function(s1_raster, polygon){
  extracted_values <- terra::extract(s1_raster, polygon)
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

# Function for s1 statistics for treatments and associated control areas

s1_statistics <- function(fgdb, treatment_layer, control_layer, imageFiles,treatmentSelected, band){
  
  # Load land cover mask raster
  lc_mask <- rast(landcover)
  
  # Load treatments
  treatments <- st_read(fgdb,layer = treatment_layer) %>% rename(TreatmentID = UniqueID)
  
  # Load treatment selection file
  treatmentSelect <- read.csv(treatmentSelected) %>% group_by(TreatmentID) %>% slice(1) %>% select(c(TreatmentID,ControlID))
  
  # Get unique treatment IDs
  treatmentID <- unique(treatmentSelect$TreatmentID)
  
  # Join treatments and select
  treatments <- treatments %>% filter(TreatmentID %in% treatmentID)
  treatments <- left_join(treatments,treatmentSelect,by="TreatmentID") 
  
  # Load controls
  control <- st_read(fgdb,layer = control_layer)
  
  # List for change data
  treatments1Change <- list()
  
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
      
      # Get image mosaic for start and stop years for selected band
      imagesStart <- imageFiles[grep(StartYear,imageFiles)]
      rasterStart <- vrt(imagesStart,options=c("-separate","-b",band,"-vrtnodata",0))
      imagesEnd <- imageFiles[grep(EndYear,imageFiles)]
      rasterEnd <- vrt(imagesEnd,options=c("-separate","-b",band,"-vrtnodata",0))
    
      # S1 NBI change for treatments using start/end years
      rasterStartCrop <- max(rasterCrop(rasterStart,polygon),na.rm=TRUE)
      rasterEndCrop <- max(rasterCrop(rasterEnd,polygon),na.rm=TRUE)
      
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
      
      # Convert dB to linear scale
      rasterStartCropLinear <- 10^(rasterStartCropMasked/ 10)
      rasterEndCropLinear <- 10^(rasterEndCropMasked / 10)
      
      # Statistics
      statsStart <- pixel_statistics(rasterStartCropLinear,polygon)
      statsEnd <- pixel_statistics(rasterEndCropLinear,polygon)
      TreatmentChangePoly <- (statsEnd - statsStart) / (statsStart + statsEnd) # polygon level NDBI
      rasterChange <- (rasterEndCropLinear - rasterStartCropLinear) / (rasterStartCropLinear + rasterEndCropLinear)
      TreatmentChangePixel <- pixel_statistics(rasterChange,polygon) # pixel level NDBI
      
      # Get control associated with treatment
      cpolygon <- control %>% filter(UniqueID == polygon$ControlID)
      
      if (dim(cpolygon)[1] > 0){
        rasterStartCrop <- max(rasterCrop(rasterStart,cpolygon),na.rm=TRUE)
        rasterEndCrop <- max(rasterCrop(rasterEnd,cpolygon),na.rm=TRUE)
        
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
        
        # Convert dB to linear scale
        rasterStartCropLinear <- 10^(rasterStartCropMasked / 10)
        rasterEndCropLinear <- 10^(rasterEndCropMasked / 10)
        
        # Statistics
        statsStart <- pixel_statistics(rasterStartCropLinear,cpolygon)
        statsEnd <- pixel_statistics(rasterEndCropLinear,cpolygon)
        ControlChangePoly <- (statsEnd - statsStart) / (statsStart + statsEnd) # polygon level NDBI
        rasterChange <- (rasterEndCropLinear - rasterStartCropLinear) / (rasterStartCropLinear + rasterEndCropLinear)
        ControlChangePixel <- pixel_statistics(rasterChange,cpolygon) # pixel level NDBI
      } else {
        ControlChangePoly <- rep(NA, 6)
        ControlChangePixel <- rep(NA, 6)
      }
      
      # Add to list
      treatments1Change[[length(treatments1Change) + 1]] <- c("s1ndbi",t,TreatmentChangePoly,TreatmentChangePixel,
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

  cnames <- c("Metric","TreatmentID",
              "TPolys1min","TPolys11stQ","TPolys1median","TPolys1mean","TPolys13rdQ","TPolys1max",
              "TPixels1min","TPixels11stQ","TPixels1median","TPixels1mean","TPixels13rdQ","TPixels1max",
              "CPolys1min","CPolys11stQ","CPolys1median","CPolys1mean","CPolys13rdQ","CPolys1max",
              "CPixels1min","CPixels11stQ","CPixels1median","CPixels1mean","CPixels13rdQ","CPixels1max")
  
  treatmentStats <- data.frame()
  for (i in 1:length(treatments1Change)){
    df <- as.data.frame(t(treatments1Change[[i]]))
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

# Get list of .tif files
imageFiles <- list.files(s1Dir, pattern = "S1GRD_mean_VV_VH", full.names = TRUE)

# Filter images to match years
pattern <- paste(years, collapse = "|")
imageFiles <- imageFiles[grepl(pattern, imageFiles)]

# Process s1 statistics 
treatment_layer <- "treatment_polygons_250119"
control_layer <- "control_polygons_250117"
treatmentStatsVV <- s1_statistics(fgdb,treatment_layer,control_layer,imageFiles,treatmentSelected,1)
treatmentStatsVH <- s1_statistics(fgdb,treatment_layer,control_layer,imageFiles,treatmentSelected,2)

# Combine data and write to .csv file
out <- paste0(outDir,"/", outFileVV)
write.csv(treatmentStatsVV,out,row.names=F)
out <- paste0(outDir,"/", outFileVH)
write.csv(treatmentStatsVH,out,row.names=F)
