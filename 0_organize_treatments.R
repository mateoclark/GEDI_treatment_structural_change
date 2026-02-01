# Organizes treatment polygons

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 01/19/2025

library(lubridate)
library(sf)
library(dplyr)

##### Parameters #####
# Version of processing
version <- "250119"

# File geodatabase with treatments
fgdb <- "E:/active/project/calfire_gedi/treatments/treatment_arcgis/gedi_treatments.gdb"

# EPSG code for final layer
epsg_code = 3310 # California Teale NAD83 meters

##### Lookup tables #####

# Reduce treatment names
treatmentLookup <- data.frame(
  TreatmentType = c(1,2,3,4,5,6,7,8),
  TreatmentClass = c(
    "Prescribed Fire",
    "Mechanical and Hand Fuels Reduction",
    "Timber Harvest",
    "Fuel Reduction",
    "Precommercial or Commercial Thin",
    "Clear Cut / Other Intense",
    "Uncertain",
    "N/A"
  ),
  TreatmentClassReduced = c(
    "Prescribed Fire",
    "Fuel Reduction",
    "Timber Harvest",
    "Fuel Reduction",
    "Timber Harvest",
    "Clear Cut",
    "Uncertain",
    "N/A"
  ),
  TreatmentSeverity = c(
    "Low",
    "Low",
    "Moderate",
    "Low",
    "Moderate",
    "High",
    "Uncertain",
    "N/A"
  )
)

# FACTS to Plant treatment names
factsLookup <- data.frame(
  ACTIVITY_N = c(
    "Commercial Thin",
    "Group Selection Cut (UA/RH/FH)",
    "Improvement Cut",
    "Permanent Land Clearing",
    "Precommercial Thin",
    "Prune",
    "Salvage Cut (intermediate treatment, not regeneration)",
    "Sanitation Cut",
    "Single-tree Selection Cut (UA/RH/FH)",
    "Stand Clearcut (EA/RH/FH)",
    "Stand Clearcut (w/ leave trees) (EA/RH/FH)",
    "Tree Release and Weed"
  ),
  TreatmentClass = c(
    "Precommercial or Commercial Thin", #"Commercial Thin"
    "Timber Harvest", #"Group Selection Cut (UA/RH/FH)"
    "Timber Harvest", #"Improvement Cut"
    "N/A", #"Permanent Land Clearing"
    "Precommercial or Commercial Thin", #"Precommercial Thin"
    "Mechanical and Hand Fuels Reduction", #"Prune"
    "Timber Harvest", #"Salvage Cut (intermediate treatment, not regeneration)"
    "Timber Harvest", #"Sanitation Cut"
    "Timber Harvest", #"Single-tree Selection Cut (UA/RH/FH)"
    "Clear Cut / Other Intense", #"Stand Clearcut (EA/RH/FH)"
    "Clear Cut / Other Intense", #"Stand Clearcut (w/ leave trees) (EA/RH/FH)"
    "Mechanical and Hand Fuels Reduction" #"Tree Release and Weed"
  )
)

sourceLookup <- data.frame(
  Source = seq(1:15),
  SourceType = c(
    "ITS",
    "Landtrendr",
    "CalMAPPER/Landtrendr",
    "CalMAPPER",
    "FACTS",
    "ITS/Landtrendr",
    "CalMAPPER/ITS",
    "CalMAPPER/ITS/Landtrendr",
    "ITS/FACTS",
    "CalMAPPER/FACTS",
    "CalMAPPER/ITS/FACTS",
    "FACTS/Landtrendr",
    "CalMAPPER/ITS/FACTS/Landtrendr",
    "ITS/FACTS/Landtrendr",
    "None"
  )
)

### Functions ####

# Function to process treatments
prepare_treatment <- function(treatment_layer,fc){
  treatments <- st_read(fgdb,layer = treatment_layer)
  treatmentsTeale <- st_transform(treatments, crs = 3310)
  treatmentsTeale$UniqueID <- paste0("t",fc,"_", seq(1, by = 1, length.out = nrow(treatmentsTeale)))
  return(treatmentsTeale)
}

##### Begin processing #####

# Process treatments #
treatment_layer <-"Treatments_planet_evaluation"
fc <- "planet"
TreatmentsPlanet <- prepare_treatment(treatment_layer,fc)

treatment_layer <- "FACTS_FuelsHazard_TimberHarvest_Silviculture_AfterJuly2019_updated_250117"
fc <- "facts"
TreatmentsFACTS <- prepare_treatment(treatment_layer,fc)

# Find distinct
dTreatment <- distinct(TreatmentsPlanet,.keep_all= TRUE)
dFACTSTreatment <- distinct(TreatmentsFACTS,.keep_all= TRUE)

# Harmonize data
dTreatment <- left_join(dTreatment,treatmentLookup,by = "TreatmentType")
dTreatment <- left_join(dTreatment,sourceLookup,by = "Source")
columnSelect <- c("RetentionPatch_or_Treatment","StartYear","EndYear","TreatmentType",                                    
                  "UniqueID","TreatmentClass",             
                  "TreatmentClassReduced","TreatmentSeverity","SourceType","Notes")            
dTreatment <- dTreatment %>% select(!!columnSelect)
dFACTSTreatment <- left_join(dFACTSTreatment,factsLookup,by = "ACTIVITY_N")
dFACTSTreatment <- left_join(dFACTSTreatment,treatmentLookup,by = "TreatmentClass") 
dFACTSTreatment$SourceType = "Raw FACTS"
dFACTSTreatment$Notes = paste0(dFACTSTreatment$METHOD_DES,"_",dFACTSTreatment$LAND_SUI_1,"_",dFACTSTreatment$EQUIPMENT1)
dFACTSTreatment$RetentionPatch_or_Treatment = 1 # assume FACTS have no retention patches (needs review)
dFACTSTreatment <- dFACTSTreatment %>% select(!!columnSelect)
treatments <- rbind(dTreatment,dFACTSTreatment) %>% filter(!is.na(StartYear) & !is.na(EndYear))

# find identical geometry and select entry with highest severity
treatments$TreatmentSeverity <- factor(treatments$TreatmentSeverity, levels = c("N/A","Uncertain","Low", "Moderate", "High"))
treatments <- treatments %>%
  group_by(SHAPE) %>% 
  slice_max(order_by = TreatmentSeverity) %>% 
  slice(1) %>%
  ungroup()                        

# Add fire information
fires <- st_read(fgdb,layer = "Fires2013_thru_2023_simplified") %>% 
  select(c(YEAR_,ALARM_DATE,CONT_DATE,FIRE_NAME)) %>%
  rename(
    FireYear = YEAR_,
    FireStartDate = ALARM_DATE,
    FireEndDate = CONT_DATE,
    FireName = FIRE_NAME
  )
firesTeale <- st_transform(fires, crs = epsg_code)
treatmentFire <- st_join(treatments, firesTeale,left = TRUE)

# Label treatments as having fire before, after or during treatment
treatmentFire$TreatmentStartDate <- as.Date(paste0(treatmentFire$StartYear,"-07-01"))
treatmentFire$TreatmentEndDate <- as.Date(paste0(treatmentFire$EndYear,"-07-31"))

treatmentFire <- treatmentFire %>%
  mutate(FireTiming = case_when(
    TreatmentEndDate < FireStartDate ~ "Fire after treatment",
    TreatmentStartDate > FireEndDate ~ "Fire before treatment",
    (TreatmentEndDate >= FireStartDate & TreatmentStartDate <= FireEndDate) ~ "Fire during treatment",
    (FireEndDate < as.Date("2019-01-01") | FireStartDate < as.Date("2019-01-01")) ~ "Fire before treatment",
    is.na(FireStartDate) ~ "No fire",
    TRUE ~ "Other" 
  ))

treatmentFinal <- treatmentFire %>% filter(FireTiming != "Fire during treatment" & TreatmentSeverity != "N/A") %>%
                  select(-c(RetentionPatch_or_Treatment,TreatmentType))

# write geodatabase
st_write(treatmentFinal, dsn = fgdb, layer = paste0("treatment_polygons_",version), append = F)

# write shapefile
treatmentShp <- treatmentFinal %>% 
  rename(
  Class = TreatmentClass,
  Reduced = TreatmentClassReduced,
  Severity = TreatmentSeverity,
  Source = SourceType,
  Timing = FireTiming
  ) %>% 
  select(-c(TreatmentStartDate,TreatmentEndDate,FireStartDate,FireEndDate))
treatmentShp <- st_zm(treatmentShp, drop = TRUE)
outshp <- paste0("E:/active/project/calfire_gedi/treatments/shapefiles/treatment_polygons_",version,".shp")
st_write(treatmentShp, outshp, append = F,SHPT="POLYGON")
