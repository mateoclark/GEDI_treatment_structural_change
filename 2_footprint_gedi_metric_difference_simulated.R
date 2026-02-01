# Calculates average change metric differences for treatments and nearest control area for similar dates
# This version is for simulated data

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 01/06/2026

library(dplyr)
library(readr)
library(purrr)
library(tidyr)

# Set the seed
set.seed(123)

# Simulated GEDI lookup 
metricLookup <- data.frame(
  simMetric = c("rhReal_25","rhReal_50","rhReal_75","rhReal_98","tLAI0t10mean","FHD_tLAI_1m","ALS_cover","aVDR","tAGBD","tLAI_1m","true_ground",
                "rhGauss_25","rhGauss_50","rhGauss_75","rhGauss_98","gLAI0t10mean","FHD_gLAI_1m","cover","gaVDR","gAGBD","gLAI_1m","gHeight"),
  onorbitMetric = c("strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98","strct_mPAI_b10","strct_FHD","strct_cover","strct_VDR","strct_tAGBD","strct_tPAI","strct_elev",
                    "gauss_RH_25","gauss_RH_50","gauss_RH_75","gauss_RH_98","gauss_mPAI_b10","gauss_FHD","gauss_cover","gauss_VDR","gauss_tAGBD","gauss_tPAI","gauss_elev")
)

##### Functions ######

# Function to calculate the closest point
find_closest <- function(mean_x_teale, mean_y_teale, df) {
  df$distance <- sqrt((df$mean_x_teale - mean_x_teale)^2 + (df$mean_y_teale - mean_y_teale)^2)
  closest <- df %>% filter(distance <= 10000) %>% slice_sample(n = 1) # allow up to 10 km of distance, pick one randomly
  if (dim(closest)[1] == 0) closest <- df[which.min(df$distance),] # if nothing is within 10 km, then go out to get the closest
  return(closest)
}

#################################
########### Parameters ##########
#################################

# Year matching -- Will ALS and GEDI data be matched, all else data filtered?
YearMatch <- "No"

# version of processing
version <- "260106"

# input directory
inDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics"

# minimum number of shots to include treatment in analysis
numShots <- 3

# minimum percent coverage by treatment polygon to include footprint in analysis
minCov <- 80

# minimum distance to the edge to include footprint in analysis
minDist <- 0 

# NLCD land cover classes to include
landcover_classes <- c(41,42,43,52) # 41 = deciduous forest; 42 = evergreen forest, 43 = mixed forest, 52 = shrub/scrub 

# Minimum NLCD land cover majority class percent
landcover_min <- 80

# Minimum pulse density (pulses/m2) 
minPulseDensity <- 3

# Maximum slope
max_slope <- 20

# Include footprints with these relationships with fire history
fireInclude <- c(
  "No fire",
  #"During treatment",
  "Before treatment and after past fire",
  #"Before treatment and before/during past fire",
  "After treatment and before post-treatment fire",
  #"After treatment and after/during post-treatment fire",
  "Other"
  )

# input .csv file with GEDI and treatment and control footprints
inputTreatmentFile <- paste0(inDir,"/data/treatment_simulated_gedi_data_",version,".csv")
inputControlFile <- paste0(inDir,"/data/control_simulated_gedi_data_",version,".csv")

# output .csv file with GEDI metric differences by treatment
outputFile <- paste0(inDir,"/data/master_simulated_gedi_treatment_difference_",version,".csv")

# output .csv files for filtered footprints
outputTreatmentFile <- paste0(inDir,"/data/treatment_simulated_gedi_data_filtered_",version,".csv")
outputControlFile <- paste0(inDir,"/data/control_simulated_gedi_data_filtered_",version,".csv")

#################################
##### Load and filter data ######
#################################

# load input GEDI data
df <- read_csv(inputTreatmentFile)

# get selected landcover percentages
landcover_area1 <- ifelse(df$shot_lc_maj1 %in% landcover_classes, (df$shot_lc_pct1/100) * df$BufferArea, 0)
landcover_area2 <- ifelse(df$shot_lc_maj2 %in% landcover_classes, (df$shot_lc_pct2/100) * df$BufferArea, 0)
landcover_area3 <- ifelse(df$shot_lc_maj3 %in% landcover_classes, (df$shot_lc_pct3/100) * df$BufferArea, 0)
landcover_area2[is.na(landcover_area2)] = 0
landcover_area3[is.na(landcover_area3)] = 0
df$landcover_pct <- (landcover_area1 + landcover_area2 + landcover_area3)/df$BufferArea * 100

# filter for percent coverage, distance to edge, timing, and pre-treatment height, 
# fire timing, not retention patch, valid landcover over minimum percent, pulse density and
# and below minimum slope
df1 <- df %>% filter(PolygonPctCovered >= minCov & 
                     distanceToTreatmentEdge >= minDist & 
                     Timing != "During treatment" & 
                     landcover_pct >= landcover_min &
                     ShotFireTiming %in% fireInclude &
                     beamDense >= minPulseDensity &
                     ground_slope <= max_slope
                     ) 
                  
# get GEDI structure metrics
rename_vec <- setNames(metricLookup$simMetric,metricLookup$onorbitMetric)
df1 <- df1 %>%
  rename(any_of(rename_vec))
metrics <- metricLookup$onorbitMetric

# get count of footprints
treatmentTable <- as.data.frame(table(df1$UniqueID, df1$Timing))
colnames(treatmentTable) <- c("UniqueID", "Timing", "FootprintCount")

# filter out those treatments without enough pre- and post-treatment data
treatmentKeep <- treatmentTable %>% group_by(UniqueID) %>%
  filter(all(c("Pre-treatment", "Post-treatment") %in% Timing) && 
           all(FootprintCount[Timing == "Pre-treatment"] >= numShots) && 
           all(FootprintCount[Timing == "Post-treatment"] >= numShots)) %>%
  ungroup()
df2 <- df1 %>% filter(UniqueID %in% treatmentKeep$UniqueID)

# get unique treatments
uniqueIDs <- unique(df2$UniqueID)

# get control footprints
df <- read_csv(inputControlFile)

# get selected landcover percentages
landcover_area1 <- ifelse(df$shot_lc_maj1 %in% landcover_classes, (df$shot_lc_pct1/100) * df$BufferArea, 0)
landcover_area2 <- ifelse(df$shot_lc_maj2 %in% landcover_classes, (df$shot_lc_pct2/100) * df$BufferArea, 0)
landcover_area3 <- ifelse(df$shot_lc_maj3 %in% landcover_classes, (df$shot_lc_pct3/100) * df$BufferArea, 0)
landcover_area2[is.na(landcover_area2)] = 0
landcover_area3[is.na(landcover_area3)] = 0
df$landcover_pct <- (landcover_area1 + landcover_area2 + landcover_area3)/df$BufferArea * 100

# filter valid landcover over minimum percent, pulse density, and below minimum slope
df3 <- df %>% filter(landcover_pct >= landcover_min & 
                       beamDense >= minPulseDensity &
                       ground_slope <= max_slope
                      )
  
# rename columns
df3 <- df3 %>%
  rename(any_of(rename_vec))

# Rename teale columns to fit functions
df2 <- df2 %>% rename(
  x_teale = Teale_X,
  y_teale = Teale_Y
)
df3 <- df3 %>% rename(
  x_teale = Teale_X,
  y_teale = Teale_Y
)

# Year matched data
if (YearMatch == "Yes"){
  df2 <- df2 %>% filter(YearALS == YearGEDI)
  df3 <- df3 %>% filter(YearALS == YearGEDI)
}

control <- df3 %>%
  group_by(UniqueID) %>%
  mutate(Timing = case_when(
    YearALS == min(YearALS, na.rm = TRUE) ~ "Year1",
    YearALS == max(YearALS, na.rm = TRUE) ~ "Year2",
    is.na(YearALS) ~ NA_character_,
    TRUE ~ NA_character_
  )) %>%
  ungroup()

# Find controls with adequate Year1 and Year2 shots
controlTable <- as.data.frame(table(control$UniqueID, control$Timing))
colnames(controlTable) <- c("UniqueID", "Timing", "FootprintCount")
controlKeep <- controlTable %>% group_by(UniqueID) %>%
  filter(all(FootprintCount[Timing == "Year1"] >= numShots) &&
           all(FootprintCount[Timing == "Year2"] >= numShots)) %>%
  ungroup()
controlKeep <- controlKeep %>% filter(Timing %in% c("Year1","Year2"))

control <- control %>% filter(UniqueID %in% controlKeep$UniqueID)
controlCentroid <- control %>% select(c(UniqueID,x_teale,y_teale)) %>%
  group_by(UniqueID) %>%
  summarize(
    mean_x_teale = mean(x_teale),
    mean_y_teale = mean(y_teale),
    .groups = "drop"
  )

#################################
########### Main loop ###########
#################################

# loop through treatments and calculate statistics for treatment and nearest control with data
df4 <- data.frame()
closestUsed <- c()
for (treatment in uniqueIDs){
  
  dt <- df2 %>% filter(UniqueID == treatment)
  cat(paste0("Treatment: ",treatment,"\n"))
  
  ## Calculate treatment statistics
  preTreatment <- dt %>% filter(Timing == "Pre-treatment")
  statsPre <- map_dfr(preTreatment[, metrics], ~ summary(.x), .id = "metric")
  statsPre <- statsPre[,1:7] %>% mutate(across(-metric, ~ as.numeric(.)))
  statsPre <- statsPre %>%
    pivot_longer(-metric, names_to = "Statistic", values_to = "value") %>%
    pivot_wider(names_from = metric, values_from = value)
  
  
  postTreatment <- dt %>% filter(Timing == "Post-treatment")
  statsPost <- map_dfr(postTreatment[, metrics], ~ summary(.x), .id = "metric")
  statsPost <- statsPost[,1:7] %>% mutate(across(-metric, ~ as.numeric(.)))
  statsPost <- statsPost %>%
    pivot_longer(-metric, names_to = "Statistic", values_to = "value") %>%
    pivot_wider(names_from = metric, values_from = value) 
  
  statsDelta <- statsPost
  statsDelta[,-1] <- statsPost[,-1] - statsPre[,-1]
  
  StartYear <- preTreatment$StartYear[1]
  EndYear <- postTreatment$EndYear[1]
  
  statsPre$n <- dim(preTreatment)[1]
  statsPre <- statsPre %>% rename_with(~ paste0("pre_", .x), .cols = -Statistic)
  statsPost$n <- dim(postTreatment)[1]
  statsPost <- statsPost %>% rename_with(~ paste0("post_", .x), .cols = -Statistic) %>% select(-Statistic)
  statsDelta <- statsDelta %>% rename_with(~ paste0("delta_", .x), .cols = -Statistic) %>% select(-Statistic)

  statsTreatment <- cbind(statsPre,statsPost,statsDelta)
  statsTreatment$Type <- "Treatment"
  statsTreatment$StartYear <- StartYear
  statsTreatment$EndYear <- EndYear
  statsTreatment$Mean_x_teale <- mean(dt$x_teale)
  statsTreatment$Mean_y_teale <- mean(dt$y_teale)
  statsTreatment$poly_lc_maj1 <- dt$poly_lc_maj1[1]
  statsTreatment$poly_lc_pct1 <- dt$poly_lc_pct1[1]
  statsTreatment$poly_lc_maj2 <- dt$poly_lc_maj2[1]
  statsTreatment$poly_lc_pct2 <- dt$poly_lc_pct2[1]
  statsTreatment$poly_lc_maj3 <- dt$poly_lc_maj3[1]
  statsTreatment$poly_lc_pct3 <- dt$poly_lc_pct3[1]
  statsTreatment$TreatmentClass <- dt$TreatmentClassReduced[1]
  statsTreatment$TreatmentSeverity <- dt$TreatmentSeverity[1]
  statsTreatment$TreatmentID <- treatment
  statsTreatment$preYear <- preTreatment$YearALS[1]
  statsTreatment$postYear <- postTreatment$YearALS[1]
  
  ## Calculate control statistics
  controlCentroid <- controlCentroid %>% filter(!UniqueID %in% closestUsed) # filter out those already used
  
  if (dim(controlCentroid)[1] > 0){ # if enough controls
    closestControl <- find_closest(statsTreatment$Mean_x_teale[1],statsTreatment$Mean_y_teale[1],controlCentroid)
    closestID <- closestControl$UniqueID
    closestUsed <- c(closestUsed,closestID)
    controlSelected <- df3 %>% filter(UniqueID == closestID)
    preControl <- controlSelected %>% filter(YearALS == min(controlSelected$YearALS)  & UniqueID == closestID)
    postControl <- controlSelected %>% filter(YearALS == max(controlSelected$YearALS) & UniqueID == closestID)
    
    statsPre <- map_dfr(preControl[, metrics], ~ summary(.x), .id = "metric")
    statsPre <- statsPre[,1:7] %>% mutate(across(-metric, ~ as.numeric(.)))
    statsPre <- statsPre %>%
      pivot_longer(-metric, names_to = "Statistic", values_to = "value") %>%
      pivot_wider(names_from = metric, values_from = value)
    
    statsPost <- map_dfr(postControl[, metrics], ~ summary(.x), .id = "metric")
    statsPost <- statsPost[,1:7] %>% mutate(across(-metric, ~ as.numeric(.)))
    statsPost <- statsPost %>%
      pivot_longer(-metric, names_to = "Statistic", values_to = "value") %>%
      pivot_wider(names_from = metric, values_from = value)
    
    statsDelta <- statsPost
    statsDelta[,-1] <- statsPost[,-1] - statsPre[,-1]
    
    statsPre$n <- dim(preControl)[1]
    statsPre <- statsPre %>% rename_with(~ paste0("pre_", .x), .cols = -Statistic)
    statsPost$n <- dim(postControl)[1]
    statsPost <- statsPost %>% rename_with(~ paste0("post_", .x), .cols = -Statistic) %>% select(-Statistic)
    statsDelta <- statsDelta %>% rename_with(~ paste0("delta_", .x), .cols = -Statistic) %>% select(-Statistic)
    
    statsControl <- cbind(statsPre,statsPost,statsDelta)
    statsControl$Type <- "Control"
    statsControl$StartYear <- StartYear
    statsControl$EndYear <- EndYear
    statsControl$Mean_x_teale <- closestControl$mean_x_teale
    statsControl$Mean_y_teale <- closestControl$mean_y_teale
    controlpoly <- control %>% filter(control$UniqueID == closestID) %>% slice(1) %>% select(poly_lc_maj1,poly_lc_pct1,poly_lc_maj2,poly_lc_pct2)
    statsControl$poly_lc_maj1 <- controlpoly$poly_lc_maj1
    statsControl$poly_lc_pct1 <- controlpoly$poly_lc_pct1
    statsControl$poly_lc_maj2 <- controlpoly$poly_lc_maj2
    statsControl$poly_lc_pct2 <- controlpoly$poly_lc_pct2
    statsControl$poly_lc_maj3 <- dt$poly_lc_maj3[1]
    statsControl$poly_lc_pct3 <- dt$poly_lc_pct3[1]
    statsControl$TreatmentClass <- dt$TreatmentClassReduced[1]
    statsControl$TreatmentSeverity <- dt$TreatmentSeverity[1]
    statsControl$TreatmentID <- treatment
    statsControl$preYear <- preControl$YearALS[1]
    statsControl$postYear <- postControl$YearALS[1]
    statsControl$ControlID <- closestID
    statsControl$ControlDistance <- closestControl$distance
    
    statsTreatment$ControlID <- closestID
    statsTreatment$ControlDistance <- closestControl$distance
    
    df4 <- rbind(df4,statsTreatment,statsControl)
  } else {
    statsTreatment$ControlID <- NA
    statsTreatment$ControlDistance <- NA
    df4 <- rbind(df4,statsTreatment) # else skip control if not enough
  }
}

#### Write out .csv file of difference statistics ####
write.csv(df4, outputFile, row.names = FALSE)

# final number of treatments
cat(length(unique(df4$TreatmentID)))

# Write out filtered treatment and control footprint files
ids <- unique(df4$TreatmentID)
dt <- df2 %>% filter(UniqueID %in% ids)
write.csv(dt, outputTreatmentFile, row.names = F)

ids <- unique(df4$ControlID)
dc <- df3 %>% filter(UniqueID %in% ids)
write.csv(dc, outputControlFile, row.names = F)

