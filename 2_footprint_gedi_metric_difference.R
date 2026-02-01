# Calculates average change metric differences for treatments and nearest control area for similar dates

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 07/10/2025

library(dplyr)

# Set the seed
set.seed(123)

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

# version of processing
version <- "250710"

# input directory
inDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics"

# output directory
outDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics/data"

# months for analysis
analysisMonths <- c(5:10) # May through October

# minimum number of shots to include treatment in analysis
numShots <- 3

# minimum percent coverage by treatment polygon to include footprint in analysis
minCov <- 80

# minimum distance to the edge to include footprint in analysis
minDist <- 0 

# months for analysis
analysisMonths <- c(5:10) # May through October
#analysisMonths <- c(1:12) # All months

# NLCD land cover classes to include
landcover_classes <- c(41,42,43,52) # 41 = deciduous forest; 42 = evergreen forest, 43 = mixed forest, 52 = shrub/scrub 

# Minimum NLCD land cover majority class percent
landcover_min <- 80

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


# input .csv file with GEDI and treatment and control information
inputTreatmentFile <- paste0(inDir,"/data/treatment_gedi_data_250710.csv")
inputControlFile <- paste0(inDir,"/data/control_gedi_data_250710.csv")

# output .csv file with GEDI metric differences by treatment
outputFile <- paste0(inDir,"/data/master_gedi_treatment_difference_",version,".csv")

#################################
##### Load and filter data ######
#################################

# load input GEDI data
df <- read.csv(inputTreatmentFile)

# get selected landcover percentages
landcover_area1 <- ifelse(df$shot_lc_maj1 %in% landcover_classes, (df$shot_lc_pct1/100) * df$BufferArea, 0)
landcover_area2 <- ifelse(df$shot_lc_maj2 %in% landcover_classes, (df$shot_lc_pct2/100) * df$BufferArea, 0)
landcover_area3 <- ifelse(df$shot_lc_maj3 %in% landcover_classes, (df$shot_lc_pct3/100) * df$BufferArea, 0)
landcover_area2[is.na(landcover_area2)] = 0
landcover_area3[is.na(landcover_area3)] = 0
df$landcover_pct <- (landcover_area1 + landcover_area2 + landcover_area3)/df$BufferArea * 100

# filter for percent coverage, distance to edge, timing, and pre-treatment height, months, fire timing, not retention patch
df1 <- df %>% filter(PolygonPctCovered >= minCov & 
                     distanceToTreatmentEdge >= minDist & 
                     Timing != "During treatment" & 
                     landcover_pct >= landcover_min &
                     Month %in% analysisMonths &
                     ShotFireTiming %in% fireInclude
                     ) 
                  
# get GEDI structure metrics
metrics <- colnames(df1)[grep(colnames(df1),pattern ="strct")]

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
df <- read.csv(inputControlFile)

# get selected landcover percentages
landcover_area1 <- ifelse(df$shot_lc_maj1 %in% landcover_classes, (df$shot_lc_pct1/100) * df$BufferArea, 0)
landcover_area2 <- ifelse(df$shot_lc_maj2 %in% landcover_classes, (df$shot_lc_pct2/100) * df$BufferArea, 0)
landcover_area3 <- ifelse(df$shot_lc_maj3 %in% landcover_classes, (df$shot_lc_pct3/100) * df$BufferArea, 0)
landcover_area2[is.na(landcover_area2)] = 0
landcover_area3[is.na(landcover_area3)] = 0
df$landcover_pct <- (landcover_area1 + landcover_area2 + landcover_area3)/df$BufferArea * 100

# filter for months, years and pre-treatment height
df3 <- df %>% filter(Month %in% analysisMonths &
                     Year <= max(df2$Year) &
                     landcover_pct >= landcover_min
                    ) 

#################################
########### Main loop ###########
#################################

# loop through treatments and calculate statistics for treatment and nearest control with data
df4 <- data.frame()
closestUsed <- c()
for (treatment in uniqueIDs){
  
  dt <- df2 %>% filter(UniqueID == treatment)
  region <- dt$StudyRegion[1]
  cat(paste0("Treatment: ",treatment,"\n"))
  
  ## Calculate treatment statistics
  preTreatment <- dt %>% filter(Timing == "Pre-treatment")
  statsPre <- data.frame(sapply(preTreatment[,metrics], summary))
  postTreatment <- dt %>% filter(Timing == "Post-treatment")
  statsPost <- data.frame(sapply(postTreatment[,metrics], summary))
  
  statsDelta <- statsPost - statsPre
  
  StartYear <- preTreatment$StartYear[1]
  EndYear <- postTreatment$EndYear[1]
  
  colnames(statsPre) <- paste0("pre_",colnames(statsPre))
  statsPre$pre_n <- dim(preTreatment)[1]
  colnames(statsPost) <- paste0("post_",colnames(statsPost))
  statsPost$post_n <- dim(postTreatment)[1]
  colnames(statsDelta) <- paste0("delta_",colnames(statsDelta))

  statsTreatment <- cbind(statsPre,statsPost,statsDelta)
  statsTreatment <- cbind(Statistic = rownames(statsTreatment), statsTreatment)
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
  statsTreatment$Source <- dt$Source[1]
  statsTreatment$Notes <- dt$Notes[1]
  statsTreatment$TreatmentID <- treatment
  
  ## Calculate control statistics
  
  # filter for controls with adequate pre- and post-treatment shots
  control <- df3 %>%
    mutate(Timing = case_when(
      Year < StartYear | (Year == StartYear & Month <= 7) ~ "Pre-treatment",
      Year > EndYear | (Year == EndYear & Month >= 7) ~ "Post-treatment",
      is.na(StartYear) ~ "NA",
      TRUE ~ "During treatment" 
    ))
  
  controlTable <- as.data.frame(table(control$UniqueID, control$Timing))
  colnames(controlTable) <- c("UniqueID", "Timing", "FootprintCount")
  controlKeep <- controlTable %>% group_by(UniqueID) %>%
    filter(all(FootprintCount[Timing == "Pre-treatment"] >= numShots) &&
             all(FootprintCount[Timing == "Post-treatment"] >= numShots)) %>%
    ungroup()
  controlKeep <- controlKeep %>% filter(Timing %in% c("Pre-treatment","Post-treatment"))
  
  control <- control %>% filter(UniqueID %in% controlKeep$UniqueID)
  controlCentroid <- control %>% select(c(UniqueID,x_teale,y_teale)) %>%
    group_by(UniqueID) %>%
    summarize(
      mean_x_teale = mean(x_teale),
      mean_y_teale = mean(y_teale),
      .groups = "drop"
    )
  controlCentroid <- controlCentroid %>% filter(!UniqueID %in% closestUsed) # filter out those already used
  closestControl <- find_closest(statsTreatment$Mean_x_teale[1],statsTreatment$Mean_y_teale[1],controlCentroid)
  closestID <- closestControl$UniqueID
  closestUsed <- c(closestUsed,closestID)
  preControl <- control %>% filter(Timing == "Pre-treatment" & UniqueID == closestID)
  postControl <- control %>% filter(Timing == "Post-treatment" & UniqueID == closestID)
  
  statsPre <- data.frame(sapply(preControl[,metrics], summary))
  statsPost <- data.frame(sapply(postControl[,metrics], summary))
  statsDelta <- statsPost - statsPre
  
  colnames(statsPre) <- paste0("pre_",colnames(statsPre))
  statsPre$pre_n <- dim(preControl)[1]
  colnames(statsPost) <- paste0("post_",colnames(statsPost))
  statsPost$post_n <- dim(postControl)[1]
  colnames(statsDelta) <- paste0("delta_",colnames(statsDelta))
  
  statsControl <- cbind(statsPre,statsPost,statsDelta)
  statsControl <- cbind(Statistic = rownames(statsControl), statsControl)
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
  statsControl$Source <- dt$Source[1]
  statsControl$Notes <- dt$Notes[1]
  statsControl$TreatmentID <- treatment
  statsControl$ControlID <- closestID
  statsControl$ControlDistance <- closestControl$distance
  
  statsTreatment$ControlID <- closestID
  statsTreatment$ControlDistance <- closestControl$distance
  
  df4 <- rbind(df4,statsTreatment,statsControl)
  
}

#### Write out .csv file of difference statistics ####
write.csv(df4, outputFile, row.names = FALSE)

# final number of treatments
cat(length(unique(df4$TreatmentID)))
