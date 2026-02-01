# Repeat-measures ANOVA, or mixed-effects model; for simulated GEDI data

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 06/01/2026


library(tidyverse)
library(rstatix)
library(performance)
library(nlme)
library(sp)
library(dplyr)
library(broom.mixed)

#### Parameters ####
# date for output file
version <- "260106"

# input directory
inDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics"

# input file
inputFile <- paste0(inDir,"/data/master_simulated_gedi_treatment_difference_",version,".csv")

# variable to use for grouping treatments
treatmentGroup <- "TreatmentSeverity"

# outlier filtering (above/below this value filtered)
zoutlier <- 4

#### Functions #####

# Function to rescale using z score and filter outliers with abs(z) > 4
zscore<-function(df,z=4){
  tmp <- as.matrix(apply(as.matrix(df),2,scale))
  if (ncol(tmp)==1) {tmp[which(abs(tmp)>z)] <- NA}
  if (ncol(tmp)>1) {tmp[which(abs(tmp)>z,arr.ind=T)] <- NA}
  return(tmp)
}

# function to preform spatial ANOVA
processModel <- function(d, metricSelect,outputFile){
  
  output <- data.frame()
  for (strctVar in metricSelect){
    
    # strctVar = metrics$metric[89] 
    print(strctVar)
    
    # gather data for analysis
    data <- as_tibble(d %>% filter(Metric == !!strctVar) %>% mutate_at(vars(index,Type,Timing,TreatmentGroup), as.factor))
    
    # Zscore normalization of metric values
    dataZscore <- data %>% mutate(Zscore = zscore(data$values,zoutlier)) %>% na.omit()
    
    # Small jitter to account for duplicate coordinates in duplicated controls
    dataZscore <- dataZscore %>%
      mutate(
        Mean_x_teale = jitter(Mean_x_teale, amount = 1e-6),
        Mean_y_teale = jitter(Mean_y_teale, amount = 1e-6)
      )
    
    tryCatch({
      aspatial <- nlme::lme(Zscore ~ Type * Timing, data = dataZscore, random = ~ 1|index, control = lmeControl(opt = "optim", maxIter = 1000, msMaxIter = 1000))
      spatial <- update(aspatial, correlation = corLin(form=~ Mean_x_teale + Mean_y_teale,nugget=F))
      
      aspatialResults <- tidy(aspatial)
      aspatialInterceptEst <- aspatialResults$estimate[1]
      aspatialTreatmentEst <- aspatialResults$estimate[2]
      aspatialTimeEst <- aspatialResults$estimate[3]
      aspatialTreatmentTimeEst <- aspatialResults$estimate[4]
      aspatialInterceptPval <- aspatialResults$p.value[1]
      aspatialTreatmentPval <- aspatialResults$p.value[2]
      aspatialTimePval <- aspatialResults$p.value[3]
      aspatialTreatmentTimePval <- aspatialResults$p.value[4]
      aspatialAIC <- AIC(aspatial)
      
      spatialResults <- tidy(spatial)
      spatialInterceptEst <- spatialResults$estimate[1]
      spatialTreatmentEst <- spatialResults$estimate[2]
      spatialTimeEst <- spatialResults$estimate[3]
      spatialTreatmentTimeEst <- spatialResults$estimate[4]
      spatialInterceptPval <- spatialResults$p.value[1]
      spatialTreatmentPval <- spatialResults$p.value[2]
      spatialTimePval <- spatialResults$p.value[3]
      spatialTreatmentTimePval <- spatialResults$p.value[4]
      spatialAIC <- AIC(spatial)
      
      anovaModelPval <- anova(aspatial,spatial)$`p-value`[2]
      
      n <- dim(dataZscore)[1]
      
      results <- data.frame(strctVar,n,anovaModelPval, 
                            aspatialInterceptEst,aspatialTreatmentEst,aspatialTimeEst,aspatialTreatmentTimeEst,
                            aspatialInterceptPval,aspatialTreatmentPval,aspatialTimePval,aspatialTreatmentTimePval,aspatialAIC,
                            spatialInterceptEst,spatialTreatmentEst,spatialTimeEst,spatialTreatmentTimeEst,
                            spatialInterceptPval,spatialTreatmentPval,spatialTimePval,spatialTreatmentTimePval,spatialAIC
      )
      
      colnames(results) <- c("metric","n","ModelTestPVal",
                             "aspatialInterceptEst","aspatialTreatmentEst","aspatialTimeEst","aspatialTreatmentTimeEst",
                             "aspatialInterceptPval","aspatialTreatmentPval","aspatialTimePval","aspatialTreatmentTimePval","aspatialAIC",
                             "spatialInterceptEst","spatialTreatmentEst","spatialTimeEst","spatialTreatmentTimeEst",
                             "spatialInterceptPval","spatialTreatmentPval","spatialTimePval","spatialTreatmentTimePval","spatialAIC"
      )
      
      # bind to output table
      output <- rbind(output, results)
      
    }, error = function(e) {
      # Handle error by skipping and printing a message
      cat("Error occurred: Skipping ANOVA")
      cat("  message: ", e$message, "\n", sep = "")
      cat("  call:    ", deparse(e$call), "\n", sep = "")
    })
  }
  
  # write out repeated ANOVA table
  write.csv(output, outputFile, row.names=F)
}

######################
#### Process data ####
######################

# read in data
df <- read.csv(inputFile) 
df <- df %>% rename(TreatmentGroup = !!treatmentGroup)

# selected metrics 
metricSelect <- colnames(df)[grep(colnames(df),pattern ="pre_strct")]
metricSelect <- gsub("pre_", "", metricSelect)

# filter to mean data
dfFiltered <- df %>% filter(Statistic == "Mean")
#dfFiltered <- df %>% filter(Statistic == "3rd Qu.")

# prepare data
postVars <- c(paste("post",metricSelect,sep="_"))
d1 <- dfFiltered %>% select(!!postVars,Type,TreatmentGroup,Mean_x_teale,Mean_y_teale)
colnames(d1)[1:length(metricSelect)] <- metricSelect
ncol <- dim(d1)[2]
d1<- cbind(d1[(ncol-3):ncol], stack(d1[1:(ncol-4)]))
d1$Timing <- "Post-treatment"
d1 <- tibble::rowid_to_column(d1, "index") 

preVars <- c(paste("pre",metricSelect,sep="_"))
d2 <- dfFiltered %>% select(!!preVars,Type,TreatmentGroup,Mean_x_teale,Mean_y_teale)
colnames(d2)[1:length(metricSelect)] <- metricSelect
ncol <- dim(d2)[2]
d2<- cbind(d2[(ncol-3):ncol], stack(d2[1:(ncol-4)]))
d2$Timing <- "Pre-treatment"
d2 <- tibble::rowid_to_column(d2, "index") 

d <- rbind(d1,d2)
d <- d %>% rename(Metric = ind) 

d$Type <- factor(d$Type, levels=c("Control","Treatment"))
d$Timing <- factor(d$Timing, levels=c("Pre-treatment","Post-treatment"))

# model all data
outputFile <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_All_simulated_",version,".csv")
processModel(d,metricSelect,outputFile)

# model treatment groups 
for (group in unique(d$TreatmentGroup)){
  print(group)
  outputFile <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_",group,"_simulated_",version,".csv")
  dg <- d %>% filter(TreatmentGroup == group)
  processModel(dg,metricSelect,outputFile)
}
