# Univariate spatial regression models, simulated data

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 01/06/2026


library(spatialreg)
library(sp)
library(spdep)
library(dplyr)
library(sf)
library(foreach)
library(doParallel)
library(tidyr)

#### Parameters ####

# input directory
inDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics"

# Detect the number of available cores
#numCores <- detectCores() - 1
numCores <- 10

# version date
version <- "260106"

# Number of nearest neighbors
k <- 10 

# input .csv file with GEDI metric differences
inputFile <- paste0(inDir,"/data/master_simulated_gedi_treatment_difference_260106.csv")

# input predictor files
landsat89NDVI <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_l89ndvi_stats_250731.csv")) %>%
  filter(!is.na(TPolyl89ndvimean))
landsat89NDMI <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_l89ndmi_stats_250731.csv")) %>%
  filter(!is.na(TPolyl89ndmimean))
sentinel2NDVI <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_s2ndvi_stats_250731.csv")) %>%
  filter(!is.na(TPolys2ndvimean))
sentinel2NDMI <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_s2ndmi_stats_250731.csv")) %>%
  filter(!is.na(TPolys2ndmimean))
gediCCDC <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_ccdc_stats_250805.csv")) %>%
  filter(!is.na(TPolyccdcmean))
sentinel1NDBIVV <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_s1ndbi_vv_stats_250731.csv")) %>%
  filter(!is.na(TPolys1mean)) 
sentinel1NDBIVH <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_s1ndbi_vh_stats_250731.csv")) %>%
  filter(!is.na(TPolys1mean)) 
palsar2NDBIHH <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_palsar2ndbi_hh_stats_250731.csv")) %>%
  filter(!is.na(TPolypalsar2mean)) 
palsar2NDBIHV <- read.csv(paste0(inDir,"/data/treatment_simulated_gedi_palsar2ndbi_hv_stats_250731.csv")) %>%
  filter(!is.na(TPolypalsar2mean))

# lookup table to go from CCDC output metric names to GEDI metric names
ccdcLookup <- data.frame(
  #Metric = c("agbd","cover","fhd-pai-1m","num-modes","pai","rh-98") ,
  Metric = c("aVDR","cover","fhd_pai_1m_a0","mPAI_b10","rh_98","agbd","pai_a0"),
  ccdcDeltaMetric = c("delta_ccdc_VDR","delta_ccdc_cover","delta_ccdc_FHD","delta_ccdc_mPAI_b10",
                     "delta_ccdc_RH_98","delta_ccdc_tAGBD","delta_ccdc_tPAI")
)
     
# severity land cover classes to include
landcover_classes <- c(41,42,43,52) # 41 = deciduous forest; 42 = evergreen forest, 43 = mixed forest, 52 = shrub/scrub 

#### Functions ####

# Get coefficients from model
model_coefficients <- function(model) {
  require(spatialreg)
  
  # Get model summary
  s <- summary(model)
  
  # Extract and convert coefficient matrix
  coefs <- as.data.frame(s$Coef)
  coefs$term <- rownames(coefs)
  
  # Reorder columns
  coefs <- coefs[, c("term", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
  
  # Add significance stars
  coefs$signif <- symnum(coefs$`Pr(>|z|)`,
                         corr = FALSE, na = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         symbols = c("***", "**", "*", ".", " "))
  
  rownames(coefs) <- NULL
  return(coefs)
}


# Define a function to process each predictor and metric
process_predictor_metric <- function(predictor, metricsDelta, dataFiltered, inDir) {
  output <- data.frame()
  coeffs <- data.frame()
  
  for (i in 1:length(metricsDelta)) {
    response <- metricsDelta[i]
    if (predictor == "PreTreatment") {
      p <- gsub("delta_strct","pre_strct",response)
      d <- data.frame(response = dataFiltered[[response]],
                      predictor = dataFiltered[[p]],
                      x = dataFiltered$x,
                      y = dataFiltered$y
      )
    } else if (predictor == "deltaGediCCDC") {
      p <- gsub("delta_strct","delta_ccdc",response)
      if (p %in% ccdcLookup$ccdcDeltaMetric){
        d <- data.frame(response = dataFiltered[[response]],
                        predictor = dataFiltered[[p]],
                        x = dataFiltered$x,
                        y = dataFiltered$y
        )
      } else{
        next # skip GEDI CCDC metric if it does not exist
      }
      
      
    } else {
      d <- data.frame(response = dataFiltered[[response]],
                      predictor = dataFiltered[[predictor]],
                      x = dataFiltered$x,
                      y = dataFiltered$y
      )
    }
    d <- d %>% filter(!is.infinite(response))
    d <- d %>% na.omit() 
    if ((dim(d)[1] > 10) & (sum(d$response) != 0)) {
      n <- dim(d)[1]
      coordinates(d) <- c("x","y")
      spatial_neighbors <- knearneigh(d,k = k)
      W <- knn2nb(spatial_neighbors)
      W <- nb2listw(W,style = "W")
      se_model <- errorsarlm(formula = response ~ predictor, data = d, listw = W, method="LU")
     
      ols_model <- lm(response ~ predictor, data = d)
      ols_AIC <- AIC(ols_model)
      ols_summary <- summary(ols_model)
      ols_int <- ols_summary$coefficients[1,1]
      ols_slope <- ols_summary$coefficients[2,1]
      ols_pvalslope <- ols_summary$coefficients[1,4]
      ols_pvalint <- ols_summary$coefficients[2,4]
      ols_r2 <- ols_summary$adj.r.squared
      model_summary <- summary(se_model,Nagelkerke = TRUE)
      ols_rmse <- sqrt(mean(residuals(ols_model)^2))
      ols_scaled_rmse <- ols_rmse / sd(d$response)
      
      spatial_int <-  model_summary$Coef[1,1]
      spatial_slope <-  model_summary$Coef[2,1]
      spatial_stderrint <- model_summary$Coef[1,2]
      spatial_stderrslope <- model_summary$Coef[2,2]
      spatial_pvalint <- model_summary$Coef[1,4]
      spatial_pvalslope <- model_summary$Coef[2,4]
      spatial_psuedo_r2 <- model_summary$NK
      spatial_AIC <- AIC(se_model)
      spatial_rmse <- sqrt(mean(residuals(se_model)^2))
      spatial_scaled_rmse <- spatial_rmse / sd(d$response)
     
      results <- data.frame(
        response,predictor,n,
        ols_AIC,ols_int,ols_slope,ols_pvalslope,ols_pvalint,ols_r2,ols_rmse,ols_scaled_rmse,
        spatial_AIC,spatial_int,spatial_slope,spatial_stderrint,spatial_stderrslope,spatial_pvalint,
        spatial_pvalslope,spatial_psuedo_r2,spatial_rmse,spatial_scaled_rmse)

      output <- rbind(output,results)

      remove(W)
    }
  }
  if (sum(grepl("ccdc", metricsDelta)) > 1) {
    outputFile <- paste0(inDir,"/data/univariate_spatialreg_",predictor,"_ccdc_simulated_",version,".csv")
  } else {
    outputFile <- paste0(inDir,"/data/univariate_spatialreg_",predictor,"_footprint_simulated_",version,".csv")
  }
  write.csv(output,outputFile,row.names = F)

}

###############################
######## Process data #########
###############################

# read data
data <- read.csv(inputFile)

# selected metrics 
metricSelect <- colnames(data)[grep(colnames(data),pattern ="pre_strct")]
metricSelect <- gsub("pre_", "", metricSelect)
metricSelect <- metricSelect[metricSelect != "strct_elev"]

# Pre and delta post-treatment GEDI metrics
metricsDelta <- paste0("delta_",metricSelect)
metricsPreTreatment <- paste0("pre_",metricSelect)

# filter data to include only data in treatments (not control)
dataFiltered <- data %>% filter(Statistic == "Mean" & Type == "Treatment") %>% 
                  select(any_of(c("Mean_x_teale","Mean_y_teale","TreatmentClass","TreatmentSeverity","TreatmentID",
                                  "poly_lc_maj1","poly_lc_pct1","poly_lc_maj2","poly_lc_pct2","poly_lc_maj3","poly_lc_pct3",   
                                  metricsDelta,metricsPreTreatment))) %>%
                  rename(x = Mean_x_teale,
                         y = Mean_y_teale)

# Calculate chosen severity landcover percentage
df_long <- pivot_longer(
  dataFiltered,
  cols = starts_with("poly_lc_"),
  names_to = c(".value", "group"),
  names_pattern = "poly_lc_(maj|pct)(\\d+)"
) %>% select(c(TreatmentID,maj,pct))

result <- df_long %>%
  filter(maj %in% landcover_classes) %>% # Filter for relevant classes
  group_by(TreatmentID) %>%
  summarize(forest_pct = sum(pct, na.rm = TRUE), .groups = "drop")

dataFiltered <- left_join(dataFiltered,result,by="TreatmentID")


# filter predictor datasets
landsat89NDVIFiltered <- landsat89NDVI %>% filter(TreatmentID %in% dataFiltered$TreatmentID) %>%
  select(TreatmentID,TPixell89ndvimean) %>%
  rename(deltaLandsat89NDVI = TPixell89ndvimean)
landsat89NDMIFiltered <- landsat89NDMI %>% filter(TreatmentID %in% dataFiltered$TreatmentID) %>%
  select(TreatmentID,TPixell89ndmimean) %>%
  rename(deltaLandsat89NDMI = TPixell89ndmimean)
sentinel2NDVIFiltered <- sentinel2NDVI %>% filter(TreatmentID %in% dataFiltered$TreatmentID) %>%
  select(TreatmentID,TPixels2ndvimean) %>%
  rename(deltaSentinel2NDVI = TPixels2ndvimean)
sentinel2NDMIFiltered <- sentinel2NDMI %>% filter(TreatmentID %in% dataFiltered$TreatmentID) %>%
  select(TreatmentID,TPixels2ndmimean) %>%
  rename(deltaSentinel2NDMI = TPixels2ndmimean)
gediCCDCFiltered <- gediCCDC %>% filter(TreatmentID %in% dataFiltered$TreatmentID & 
                                          !Metric %in% c("Gini","WI_b10m","pai_a0","num_modes")) %>%
  select(TreatmentID,Metric,TPixelccdcmean) %>%
  rename(deltaGediCCDC = TPixelccdcmean)
gediCCDCFiltered <- left_join(gediCCDCFiltered,ccdcLookup,by="Metric") %>% select(-Metric)
gediCCDCFiltered <- gediCCDCFiltered %>%
  distinct(TreatmentID, ccdcDeltaMetric, .keep_all = TRUE)
gediCCDCFiltered <- gediCCDCFiltered %>%
  pivot_wider(
    names_from = ccdcDeltaMetric, # Columns for the new headers
    values_from = deltaGediCCDC  # Values to fill in the table
  )
sentinel1NDBIVVFiltered <- sentinel1NDBIVV %>% filter(TreatmentID %in% dataFiltered$TreatmentID) %>%
  select(TreatmentID,TPixels1mean) %>%
  rename(deltaSentinel1NDBIVV = TPixels1mean)
sentinel1NDBIVHFiltered <- sentinel1NDBIVH %>% filter(TreatmentID %in% dataFiltered$TreatmentID) %>%
  select(TreatmentID,TPixels1mean) %>%
  rename(deltaSentinel1NDBIVH = TPixels1mean)
palsar2NDBIHHFiltered <- palsar2NDBIHH %>% filter(TreatmentID %in% dataFiltered$TreatmentID) %>%
  select(TreatmentID,TPixelpalsar2mean) %>%
  rename(deltaPalsar2NDBIHH = TPixelpalsar2mean)
palsar2NDBIHVFiltered <- palsar2NDBIHV %>% filter(TreatmentID %in% dataFiltered$TreatmentID) %>%
  select(TreatmentID,TPixelpalsar2mean) %>%
  rename(deltaPalsar2NDBIHV = TPixelpalsar2mean)


# join predictors with treatment data
dataFiltered <- left_join(dataFiltered, landsat89NDVIFiltered, by="TreatmentID")
dataFiltered <- left_join(dataFiltered, landsat89NDMIFiltered, by="TreatmentID")
dataFiltered <- left_join(dataFiltered, sentinel2NDVIFiltered, by="TreatmentID")
dataFiltered <- left_join(dataFiltered, sentinel2NDMIFiltered, by="TreatmentID")
dataFiltered <- left_join(dataFiltered, gediCCDCFiltered, by="TreatmentID")
dataFiltered <- left_join(dataFiltered, sentinel1NDBIVVFiltered, by="TreatmentID")
dataFiltered <- left_join(dataFiltered, sentinel1NDBIVHFiltered, by="TreatmentID")
dataFiltered <- left_join(dataFiltered, palsar2NDBIHHFiltered, by="TreatmentID")
dataFiltered <- left_join(dataFiltered, palsar2NDBIHVFiltered, by="TreatmentID")

# Scale CCDC to native units
dataFiltered <- dataFiltered%>% mutate(
  delta_ccdc_VDR = delta_ccdc_VDR/1000,
  delta_ccdc_cover = delta_ccdc_cover/100,
  delta_ccdc_mPAI_b10 = delta_ccdc_mPAI_b10/1000,
  delta_ccdc_RH_98 = delta_ccdc_RH_98/100,
  delta_ccdc_FHD = delta_ccdc_FHD/1000,
  delta_ccdc_tAGBD = delta_ccdc_tAGBD/10,
  delta_ccdc_tPAI = delta_ccdc_tAGBD/1000
)

#### Run models ####

# Register the parallel backend
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Run predictors with parallelization for deltaGEDI (from footprint)
predictors <- c("PreTreatment","deltaLandsat89NDVI","deltaLandsat89NDMI","deltaSentinel2NDVI","deltaSentinel2NDMI",
                "deltaGediCCDC","deltaSentinel1NDBIVV","deltaSentinel1NDBIVH","deltaPalsar2NDBIHH","deltaPalsar2NDBIHV")
foreach(predictor = predictors, .packages = c("dplyr", "sf", "sp", "spatialreg", "spdep")) %dopar% {
  process_predictor_metric(predictor, metricsDelta, dataFiltered,inDir)
}

# Run predictors with parallelization for delta CCDC GEDI (from rasters)
predictors <- c("deltaLandsat89NDVI","deltaLandsat89NDMI","deltaSentinel2NDVI","deltaSentinel2NDMI",
                "deltaSentinel1NDBIVV","deltaSentinel1NDBIVH","deltaPalsar2NDBIHH","deltaPalsar2NDBIHV")
foreach(predictor = predictors, .packages = c("dplyr", "sf", "sp", "spatialreg", "spdep")) %dopar% {
  process_predictor_metric(predictor, ccdcLookup$ccdcDeltaMetric, dataFiltered, inDir)
}

# Stop the cluster
stopCluster(cl)