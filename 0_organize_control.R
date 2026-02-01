# Organizes control polygons

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 01/17/2025

library(sf)
library(dplyr)

##### Parameters #####
# Version of processing
version <- "250117"

# File geodatabase with treatments
fgdb <- "E:/active/project/calfire_gedi/treatments/treatment_arcgis/gedi_treatments.gdb"

# EPSG code for final layer
epsg_code = 3310 # California Teale NAD83 meters

### Functions ####

# Function to process control
prepare_control <- function(control_layer,fc){
  control <- st_read(fgdb,layer = control_layer) 
  controlTeale <- st_transform(control, crs = 3310)
  controlTeale$UniqueID <- paste0("c",fc,"_", seq(1, by = 1, length.out = nrow(controlTeale)))
  return(controlTeale)
}

##### Begin processing #####

# Process controls #
control_layer <- "control_zones_buffer500m"
fc <- "radius500m"
control<- prepare_control(control_layer,fc)

# Find distinct
dControl <- distinct(control,.keep_all= TRUE)

# write geodatabase
st_write(dControl, dsn = fgdb, layer = paste0("control_polygons_",version), append = F)

# write shapefile
controlShp <- st_zm(dControl, drop = TRUE)
outshp <- paste0("E:/active/project/calfire_gedi/treatments/shapefiles/control_polygons_",version,".shp")
st_write(controlShp, outshp, append = F,SHPT="POLYGON")
