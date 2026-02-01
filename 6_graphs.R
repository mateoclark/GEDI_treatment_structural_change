# Creates graphs for paper "Detecting structural change from forest fuel treatments with GEDI spaceborne lidar across California, USA"

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 01/09/2026

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)
library(forcats)
library(rstatix)
library(scales)
library(broom)
library(ggh4x)
library(viridis)
library(hexbin)
library(purrr)
library(spdep)
library(spatialreg)
library(reshape2)
library(readr)
library(lubridate)
library(ggtext)
library(fixest)
library(RColorBrewer)

#### Parameters ####
# working directory
inDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_treatment_change_in_footprint_metrics"

# lookup
metricLookup <- data.frame(
  metric = c("strct_FHD","strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
             "strct_mPAI_b10","strct_prop_int_below_10m","strct_RH_25",
             "strct_RH_50","strct_RH_75","strct_RH_98"),
  labels = c("FHD","H","nmode","aVDR","Cover","tAGBD","tPAI","mPAI0to10","REBelow10m",
             "RH25","RH50","RH75","RH98")
)

# Response labels

response_labels <- c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                     expression(mPAI["0-10m"]),expression(RE["<10m"]),
                     "RH25","RH50","RH75","RH98")
response_labels_rev <- c("RH98","RH75","RH50","RH25",
                         expression(RE["<10m"]),expression(mPAI["0-10m"]),
                         "tPAI","tAGBD","Cover","aVDR","nmode","FHD")

##########################################################################################################
#### Summary of treatments ####
##########################################################################################################

#df <- read.csv(paste0(inDir,"/data/master_gedi_treatment_difference_250710.csv")) # on-orbit
df <- read.csv(paste0(inDir,"/data/master_simulated_gedi_treatment_difference_260106.csv")) # simulated data

df1 <- df %>% filter(Statistic == "Mean" & Type == "Treatment")

print(length(unique(df1$TreatmentID))) # number of treatments
print(length(unique(df1$ControlID))) # number of controls

print(table(df1$StartYear)) # Start year count
print(table(df1$EndYear)) # End year count

print(table(df1$TreatmentSeverity)) # count of severity classes
print(table(df1$TreatmentClass)) # count of treatment class

print(summary(df1$ControlDistance)) # summary of control distance
d <- df1 %>% filter(!is.na(ControlDistance))
print(sd(d$ControlDistance))


print(summary(df1$pre_n)) # summary of pre-treatment number of footprints
print(sd(df1$pre_n))
print(summary(df1$post_n)) # summary of post-treatment number of footprints
print(sd(df1$post_n))


# Calculate forest cover (NLCD 2018) of treatments
landcover_classes <- c(41,42,43,52) # 41 = deciduous forest; 42 = evergreen forest, 43 = mixed forest, 52 = shrub/scrub 

df2 <- pivot_longer(
  df1,
  cols = starts_with("poly_lc_"),
  names_to = c(".value", "group"),
  names_pattern = "poly_lc_(maj|pct)(\\d+)"
) %>% select(c(TreatmentID,maj,pct))


df3 <- df2 %>%
  filter(maj %in% landcover_classes) %>% # Filter for relevant classes
  group_by(TreatmentID) %>%
  summarize(forest_pct = sum(pct, na.rm = TRUE), .groups = "drop")
df3$forest_pct <- ifelse(df3$forest_pct > 100, 100, df3$forest_pct)

df4 <- left_join(df1,df3,by="TreatmentID")
print(summary(df4$forest_pct))
print(sd(df4$forest_pct))

# Calculate forest cover (NLCD 2018) of controls
df1 <- df %>% filter(Statistic == "Mean" & Type == "Control")

df2 <- pivot_longer(
  df1,
  cols = starts_with("poly_lc_"),
  names_to = c(".value", "group"),
  names_pattern = "poly_lc_(maj|pct)(\\d+)"
) %>% select(c(TreatmentID,maj,pct))


df3 <- df2 %>%
  filter(maj %in% landcover_classes) %>% # Filter for relevant classes
  group_by(TreatmentID) %>%
  summarize(forest_pct = sum(pct, na.rm = TRUE), .groups = "drop")
df3$forest_pct <- ifelse(df3$forest_pct > 100, 100, df3$forest_pct)

df4 <- left_join(df1,df3,by="TreatmentID")
print(summary(df4$forest_pct))
print(sd(df4$forest_pct))

##########################################################################################################
#### Plot NLME Spatial Repeat-measures ANOVA for On-orbit GEDI ####
##########################################################################################################

# version date
version <- "250710"

# Levels and labels
response_levels <- c("strct_FHD","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                     "strct_mPAI_b10","strct_prop_int_below_10m",
                     "strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98")

pvalue = 0.01
plot_list = list()

# y-axis limits
ylimits = c(-1.3,1.3)

# All treatment types

file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_All_",version,".csv")
df <- read.csv(file) %>% filter(!metric %in% c("strct_RH_10", "strct_H"))

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df1 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df1) <- c("metric","ind","pval","pval-type","zscore","effect")
df1$significant <- ifelse(df1$pval <= pvalue,"Sig.","Not Sig.")
df1$effectFactor <- factor(df1$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Treatment", "Timing","Interaction"))
df1$metricFactor <- factor(df1$metric, 
                           levels = response_levels)
pal <- RColorBrewer::brewer.pal(3, "Set2")

if (length(unique(df1$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot1 <- ggplot(df1, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(
    position = "dodge", stat = "identity",
    pattern_density = 0.1,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_spacing = 0.02,
  ) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "A. All treatment types", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.key.width = unit(.1, "in"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#dev.new();print(plot1)
plot_list[[1]] = plot1

# Low impact
file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_Low_",version,".csv")
df <- read.csv(file) %>% filter(!metric %in% c("strct_RH_10", "strct_H"))

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df2 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df2) <- c("metric","ind","pval","pval-type","zscore","effect")
df2$significant <- ifelse(df2$pval <= pvalue,"Sig.","Not Sig.")
df2$effectFactor <- factor(df2$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df2$metricFactor <- factor(df2$metric, 
                           levels = response_levels)
if (length(unique(df2$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot2 <- ggplot(df2, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "B. Low impact", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text( size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot2)
plot_list[[2]] = plot2

# Moderate impact
file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_Moderate_",version,".csv")
df <- read.csv(file) %>% filter(!metric %in% c("strct_RH_10", "strct_H"))

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df3 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df3) <- c("metric","ind","pval","pval-type","zscore","effect")
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")
df3$effectFactor <- factor(df3$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df3$metricFactor <- factor(df3$metric, 
                           levels = response_levels)
if (length(unique(df3$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot3 <- ggplot(df3, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "C. Moderate impact", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot3)
plot_list[[3]] = plot3

# High impact
file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_High_",version,".csv")
df <- read.csv(file) %>% filter(!metric %in% c("strct_RH_10", "strct_H"))

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df4 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df4) <- c("metric","ind","pval","pval-type","zscore","effect")
df4$significant <- ifelse(df4$pval <= pvalue,"Sig.","Not Sig.")
df4$effectFactor <- factor(df4$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df4$metricFactor <- factor(df4$metric, 
                           levels = response_levels)
if (length(unique(df4$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot4 <- ggplot(df4, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "D. High impact") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot4)
plot_list[[4]] = plot4

# Uncertain impact
file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_Uncertain_",version,".csv")
df <- read.csv(file) %>% filter(!metric %in% c("strct_RH_10", "strct_H"))

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df5 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df5) <- c("metric","ind","pval","pval-type","zscore","effect")
df5$significant <- ifelse(df5$pval <= pvalue,"Sig.","Not Sig.")
df5$effectFactor <- factor(df5$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df5$metricFactor <- factor(df5$metric, 
                           levels = response_levels)
if (length(unique(df5$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot5 <- ggplot(df5, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "E. Not classified") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot4)
plot_list[[5]] = plot5

combined1 <- ggarrange(plotlist = plot_list,
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 2, nrow = 3,
                      align = "hv")

dev.new();print(combined1)

outFile <- paste0(inDir,"/figures/nlme_repeated_anova_spatial_zscore_",version,".png")
ggsave(outFile, width = 8, height = 7, units = "in", dpi=600)

##########################################################################################################
#### Plot NLME Spatial Repeat-measures ANOVA for Simulated GEDI ####
##########################################################################################################

# version date
version <- "260106"

# Levels and labels
response_levels_sim <- c("strct_FHD","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                     "strct_mPAI_b10","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98")
response_labels_sim <- expression("FHD", "aVDR", "Cover", "tAGBD", "tPAI", mPAI["0-10m"],
                                  "RH25", "RH50", "RH75", "RH98")

pvalue = 0.01
plot_list = list()

# y-axis limits
ylimits = c(-1.7,1.7)

# All treatment types

file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_All_simulated_",version,".csv")
df <- read.csv(file) %>% filter(metric != "strct_elev")

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df1 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df1) <- c("metric","ind","pval","pval-type","zscore","effect")
df1$significant <- ifelse(df1$pval <= pvalue,"Sig.","Not Sig.")
df1$effectFactor <- factor(df1$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Treatment", "Timing","Interaction"))
df1$metricFactor <- factor(df1$metric, 
                           levels = response_levels_sim)
pal <- RColorBrewer::brewer.pal(3, "Set2")

if (length(unique(df1$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot1 <- ggplot(df1, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(
    position = "dodge", stat = "identity",
    pattern_density = 0.1,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_spacing = 0.02,
  ) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "A. All treatment types", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels_sim) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.key.width = unit(.1, "in"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#dev.new();print(plot1)
plot_list[[1]] = plot1

# Low impact
file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_Low_simulated_",version,".csv")
df <- read.csv(file) %>% filter(metric != "strct_elev")

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df2 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df2) <- c("metric","ind","pval","pval-type","zscore","effect")
df2$significant <- ifelse(df2$pval <= pvalue,"Sig.","Not Sig.")
df2$effectFactor <- factor(df2$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df2$metricFactor <- factor(df2$metric, 
                           levels = response_levels_sim)
if (length(unique(df2$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot2 <- ggplot(df2, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "B. Low impact", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels_sim) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text( size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot2)
plot_list[[2]] = plot2

# Moderate impact
file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_Moderate_simulated_",version,".csv")
df <- read.csv(file) %>% filter(metric != "strct_elev")

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df3 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df3) <- c("metric","ind","pval","pval-type","zscore","effect")
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")
df3$effectFactor <- factor(df3$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df3$metricFactor <- factor(df3$metric, 
                           levels = response_levels_sim)
if (length(unique(df3$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot3 <- ggplot(df3, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "C. Moderate impact", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels_sim) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot3)
plot_list[[3]] = plot3

# High impact
file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_High_simulated_",version,".csv")
df <- read.csv(file) %>% filter(metric != "strct_elev")

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df4 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df4) <- c("metric","ind","pval","pval-type","zscore","effect")
df4$significant <- ifelse(df4$pval <= pvalue,"Sig.","Not Sig.")
df4$effectFactor <- factor(df4$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df4$metricFactor <- factor(df4$metric, 
                           levels = response_levels_sim)
if (length(unique(df4$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot4 <- ggplot(df4, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "D. High impact", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels_sim) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot4)
plot_list[[4]] = plot4


# Uncertain impact
file <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_Uncertain_simulated_",version,".csv")
df <- read.csv(file) %>% filter(metric != "strct_elev")

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df5 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df5) <- c("metric","ind","pval","pval-type","zscore","effect")
df5$significant <- ifelse(df5$pval <= pvalue,"Sig.","Not Sig.")
df5$effectFactor <- factor(df5$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df5$metricFactor <- factor(df5$metric, 
                           levels = response_levels_sim)
if (length(unique(df5$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot5 <- ggplot(df5, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  geom_vline(
    xintercept = seq(
      from = 1.5,
      to   = length(levels(df1$metricFactor)) - 0.5,
      by   = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  labs(title = "E. Not classified") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = response_labels_sim) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot5)
plot_list[[5]] = plot5

combined1 <- ggarrange(plotlist = plot_list,
                       common.legend = TRUE,
                       legend = "bottom",
                       ncol = 2, nrow = 3,
                       align = "hv")

dev.new();print(combined1)

outFile <- paste0(inDir,"/figures/nlme_repeated_anova_spatial_zscore_simulated_",version,".png")
ggsave(outFile, width = 8, height = 7, units = "in", dpi=600)


###########################################################################################################
### Spatial regression of change in pre- and post-treatment univariate predictor r2 and slopes for on-orbit
### data with all treatment types. This graph using change in average per-pixel change for treatment polygons
##########################################################################################################

# version date
version <- "250805"

# Levels and labels
predictor_levels <- c("deltaLandsat89NDVI","deltaLandsat89NDMI","deltaSentinel2NDVI","deltaSentinel2NDMI",  
                      "deltaPalsar2NDBIHH","deltaPalsar2NDBIHV","deltaSentinel1NDBIVV","deltaSentinel1NDBIVH","deltaGediCCDC")
predictor_labels <- c("L8/9 \u0394NDVI", "L8/9 \u0394NDMI",
                      "S2 \u0394NDVI", "S2 \u0394NDMI", "P2-HH NDBI", "P2-HV NDBI", 
                      "S1-VV NDBI", "S1-VH NDBI","\u0394GEDI Fusion")
response_levels <- c("delta_strct_FHD","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD","delta_strct_tPAI",
                     "delta_strct_mPAI_b10","delta_strct_prop_int_below_10m",
                     "delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75","delta_strct_RH_98")

# get regression results tables
Landsat89NDVI <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaLandsat89NDVI_footprint_",version,".csv"))
Landsat89NDMI <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaLandsat89NDMI_footprint_",version,".csv"))
Sentinel2NDVI <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel2NDVI_footprint_",version,".csv"))
Sentinel2NDMI <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel2NDMI_footprint_",version,".csv"))
Palsar2NDBIHH <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaPalsar2NDBIHH_footprint_",version,".csv"))
Palsar2NDBIHV <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaPalsar2NDBIHV_footprint_",version,".csv"))
Sentinel1NDBIVV <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel1NDBIVV_footprint_",version,".csv"))
Sentinel1NDBIVH <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel1NDBIVH_footprint_",version,".csv"))
GEDICCDC <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaGediCCDC_footprint_",version,".csv"))


df_orb <- rbind(Landsat89NDVI,Landsat89NDMI,Sentinel2NDVI,Sentinel2NDMI,
            Palsar2NDBIHH,Palsar2NDBIHV,Sentinel1NDBIVV,Sentinel1NDBIVH,GEDICCDC) %>% 
      filter(response != "delta_strct_H")

df_orb <- df_orb %>% filter(response %in% paste0("delta_",metricLookup$metric))

df_orb$metricFactor <- factor(df_orb$response, 
                          levels = response_levels)

pal <- RColorBrewer::brewer.pal(10, "Paired")
#pal1 <- c(pal[3:10],pal[1])

# Plot of r2
df1 <- df_orb %>% select(ols_r2,spatial_psuedo_r2,predictor,metricFactor)
df1$modelFactor <- factor(df1$predictor, 
                          levels = predictor_levels,
                          labels = predictor_labels)
plotR2 <- ggplot(df1,aes(x=metricFactor,y=spatial_psuedo_r2,fill = modelFactor)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_vline(
    xintercept = seq(
      1.5,
      length(levels(df1$metricFactor)) - 0.5,
      by = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.05)) +
  scale_fill_manual(values=pal) +
  ylab(expression("Pseudo-R "^2)) +
  scale_x_discrete(labels = response_labels) +
  #guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(.1, "in"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

#dev.new();print(plotR2)

# Plot of scaled RMSE
df2 <- df_orb %>% select(ols_scaled_rmse,spatial_scaled_rmse,predictor,metricFactor)
df2$modelFactor <- factor(df2$predictor, 
                          levels = predictor_levels,
                          labels = predictor_labels)
plotRMSE <- ggplot(df2,aes(x=metricFactor,y=spatial_scaled_rmse,fill = modelFactor)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_vline(
    xintercept = seq(
      1.5,
      length(levels(df1$metricFactor)) - 0.5,
      by = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0.7, 1.0, by = 0.05)) +
  coord_cartesian(ylim = c(0.7, 1.0)) +
  scale_fill_manual(values=pal) +
  ylab("Scaled RMSE") +
  scale_x_discrete(labels = response_labels) + 
  #guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(.1, "in"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

#dev.new();print(plotRMSE)

# combine plots into final figure
combined2 <- ggarrange(plotR2,plotRMSE,
                       labels = c("A", "B"),
                       common.legend =TRUE,
                       legend = "bottom",
                       ncol = 1, nrow = 2)

dev.new();print(combined2)

outFile <- paste0(inDir,"/figures/regression_models_univariate_r2_rmse_",version,".png")
ggsave(outFile, width = 8, height = 6, units = "in", dpi=600)

outFile <- paste0(inDir,"/data/regression_models_univariate_statistics_",version,".csv")
write.csv(df,outFile, row.names = F)

# Pivot spatial pseudo R2 values by predictors
df_spatial_r2 <- df_orb %>%
  select(response, predictor, spatial_psuedo_r2) %>%
  pivot_wider(names_from = predictor, values_from = spatial_psuedo_r2)
summary(df_spatial_r2)

# Pivot spatial RMSE values by predictors
df_spatial_rmse <- df_orb %>%
  select(response, predictor, spatial_scaled_rmse) %>%
  pivot_wider(names_from = predictor, values_from = spatial_scaled_rmse)
summary(df_spatial_rmse)

# Pivot OLS R2 values by predictors
df_ols_r2 <- df_orb %>%
  select(response, predictor, ols_r2) %>%
  pivot_wider(names_from = predictor, values_from = ols_r2)
summary(df_ols_r2)

# NDVI comparison
mean(df_spatial_r2$deltaSentinel2NDVI)
mean(df_spatial_r2$deltaLandsat89NDVI)
wilcox.test(df_spatial_r2$deltaSentinel2NDVI, df_spatial_r2$deltaLandsat89NDVI, paired = TRUE)

# NDMI comparison
mean(df_spatial_r2$deltaSentinel2NDMI)
mean(df_spatial_r2$deltaLandsat89NDMI)
wilcox.test(df_spatial_r2$deltaSentinel2NDMI, df_spatial_r2$deltaLandsat89NDMI, paired = TRUE)

# Summarize spatial_psuedo_r2 for each predictor
summary_stats_orb <- df_orb %>%
  group_by(predictor) %>%
  summarize(
    mean_r2 = mean(spatial_psuedo_r2, na.rm = TRUE),
    sd_r2 = sd(spatial_psuedo_r2, na.rm = TRUE),
    min_r2 = min(spatial_psuedo_r2, na.rm = TRUE),
    max_r2 = max(spatial_psuedo_r2, na.rm = TRUE),
    mean_rmse = mean(spatial_scaled_rmse, na.rm = TRUE),
    sd_rmse = sd(spatial_scaled_rmse, na.rm = TRUE),
    n = sum(!is.na(spatial_psuedo_r2)),
    .groups = "drop"
  )

# View the results
print(summary_stats_orb)

# Test CCDC GEDI interpolation against other predictors

data_filtered <- df_orb %>%
  filter(!is.na(spatial_psuedo_r2), !is.na(predictor))
gedi_vals <- data_filtered %>%
  filter(predictor == "deltaGediCCDC") %>%
  pull(spatial_psuedo_r2)
results <- data_filtered %>%
  filter(predictor != "deltaGediCCDC" & !response %in% c("delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75")) %>%
  group_by(predictor) %>%
  summarize(
    p_value = wilcox.test(gedi_vals, spatial_psuedo_r2, alternative = "greater")$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "holm"))

# View results
print(results)

sar_orb <- df_orb %>% filter(predictor %in% c("deltaSentinel1NDBIVV","deltaSentinel1NDBIVH","deltaPalsar2NDBIHH","deltaPalsar2NDBIHV")) %>% select(predictor, response, spatial_psuedo_r2)
ccdc_orb <- df_orb %>% filter(predictor == "deltaGediCCDC") %>% select(predictor, response, ols_r2, spatial_psuedo_r2)

biomass_orb <- df_orb %>% filter(response == "delta_strct_tAGBD") %>% select(predictor, response, spatial_psuedo_r2)
cover_orb <- df_orb %>% filter(response == "delta_strct_cover") %>% select(predictor, response, spatial_psuedo_r2)

###########################################################################################################
### Spatial regression of ALS simulation change in pre- and post-treatment univariate predictor r2 and slopes
### all treatments
### This graph using change in average per-pixel change for treatment polygons
##########################################################################################################

# version date
version <- "260106"

# Levels and labels
predictor_levels <- c("deltaLandsat89NDVI","deltaLandsat89NDMI","deltaSentinel2NDVI","deltaSentinel2NDMI",  
                      "deltaPalsar2NDBIHH","deltaPalsar2NDBIHV","deltaSentinel1NDBIVV","deltaSentinel1NDBIVH","deltaGediCCDC")
predictor_labels <- c("L8/9 \u0394NDVI", "L8/9 \u0394NDMI",
                      "S2 \u0394NDVI", "S2 \u0394NDMI", "P2-HH NDBI", "P2-HV NDBI", 
                      "S1-VV NDBI", "S1-VH NDBI","\u0394GEDI Fusion")
response_levels_sim <- c("delta_strct_tAGBD","delta_strct_FHD","delta_strct_VDR","delta_strct_cover", "delta_strct_tPAI",
                         "delta_strct_mPAI_b10","delta_strct_RH_25","delta_strct_RH_50",
                         "delta_strct_RH_75","delta_strct_RH_98")
response_labels_sim <- c("tAGBD","FHD","aVDR","Cover","tPAI", expression(mPAI["0-10m"]),"RH25","RH50","RH75","RH98")
response_labels_sim_rev <- c("RH98","RH75","RH50","RH25",
                             expression(mPAI["0-10m"]),"tPAI", "Cover","aVDR","FHD","tAGBD")

# get regression results tables
Landsat89NDVI <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaLandsat89NDVI_footprint_simulated_",version,".csv"))
Landsat89NDMI <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaLandsat89NDMI_footprint_simulated_",version,".csv"))
Sentinel2NDVI <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel2NDVI_footprint_simulated_",version,".csv"))
Sentinel2NDMI <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel2NDMI_footprint_simulated_",version,".csv"))
Palsar2NDBIHH <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaPalsar2NDBIHH_footprint_simulated_",version,".csv"))
Palsar2NDBIHV <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaPalsar2NDBIHV_footprint_simulated_",version,".csv"))
Sentinel1NDBIVV <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel1NDBIVV_footprint_simulated_",version,".csv"))
Sentinel1NDBIVH <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel1NDBIVH_footprint_simulated_",version,".csv"))
GEDICCDC <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaGediCCDC_footprint_simulated_",version,".csv"))

df_sim <- rbind(Landsat89NDVI,Landsat89NDMI,Sentinel2NDVI,Sentinel2NDMI,
                 Palsar2NDBIHH,Palsar2NDBIHV,Sentinel1NDBIVV,Sentinel1NDBIVH,GEDICCDC)

df_sim <- df_sim  %>% filter(response %in% paste0("delta_",metricLookup$metric))

df_sim$metricFactor <- factor(df_sim$response, 
                               levels = response_levels_sim)


pal <- RColorBrewer::brewer.pal(10, "Paired")

# Plot of r2
df1 <- df_sim  %>% select(ols_r2,spatial_psuedo_r2,predictor,metricFactor)
df1$modelFactor <- factor(df1$predictor, 
                          levels = predictor_levels,
                          labels = predictor_labels)
plotR2 <- ggplot(df1,aes(x=metricFactor,y=spatial_psuedo_r2,fill = modelFactor)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_vline(
    xintercept = seq(
      1.5,
      length(levels(df1$metricFactor)) - 0.5,
      by = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1)) +
  scale_fill_manual(values=pal) +
  ylab(expression("Pseudo-R "^2)) +
  scale_x_discrete(labels = response_labels_sim) +
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(.1, "in"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

#dev.new();print(plotR2)

# Plot of scaled RMSE
df2 <- df_sim  %>% select(ols_scaled_rmse,spatial_scaled_rmse,predictor,metricFactor)
df2$modelFactor <- factor(df2$predictor, 
                          levels = predictor_levels,
                          labels = predictor_labels)
plotRMSE <- ggplot(df2,aes(x=metricFactor,y=spatial_scaled_rmse,fill = modelFactor)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_vline(
    xintercept = seq(
      1.5,
      length(levels(df1$metricFactor)) - 0.5,
      by = 1
    ),
    color = "grey80",
    linewidth = 0.3
  ) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0.5, 1.0, by = 0.05)) +
  coord_cartesian(ylim = c(0.5, 1.0)) +
  scale_fill_manual(values=pal) +
  ylab("Scaled RMSE") +
  scale_x_discrete(labels = response_labels_sim) + 
  #guides(fill = guide_legend(ncol = 4, byrow = TRUE)) +
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(.1, "in"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )


# combine plots into final figure
combined2 <- ggarrange(plotR2,plotRMSE,
                       labels = c("A", "B"),
                       common.legend =TRUE,
                       legend = "bottom",
                       ncol = 1, nrow = 2)

dev.new();print(combined2)

outFile <- paste0(inDir,"/figures/regression_models_univariate_r2_rmse_simulated_",version,".png")
ggsave(outFile, width = 8, height = 6, units = "in", dpi=600)

# Pivot pseudo R2 values by predictors
df_wide <- df_sim  %>%
  select(response, predictor, spatial_psuedo_r2) %>%
  pivot_wider(names_from = predictor, values_from = spatial_psuedo_r2)
summary(df_wide)

# NDVI comparison
mean(df_wide$deltaSentinel2NDVI)
mean(df_wide$deltaLandsat89NDVI)
wilcox.test(df_wide$deltaSentinel2NDVI, df_wide$deltaLandsat89NDVI, paired = TRUE)

# NDMI comparison
mean(df_wide$deltaSentinel2NDMI)
mean(df_wide$deltaLandsat89NDMI)
wilcox.test(df_wide$deltaSentinel2NDMI, df_wide$deltaLandsat89NDMI, paired = TRUE)

# Summarize spatial_psuedo_r2 for each predictor
summary_stats_sim <- df_sim  %>%
  group_by(predictor) %>%
  summarize(
    mean_r2 = mean(spatial_psuedo_r2, na.rm = TRUE),
    sd_r2 = sd(spatial_psuedo_r2, na.rm = TRUE),
    min_r2 = min(spatial_psuedo_r2, na.rm = TRUE),
    max_r2 = max(spatial_psuedo_r2, na.rm = TRUE),
    mean_rmse = mean(spatial_scaled_rmse, na.rm = TRUE),
    sd_rmse = sd(spatial_scaled_rmse, na.rm = TRUE),
    n = sum(!is.na(spatial_psuedo_r2)),
    .groups = "drop"
  )

# View the results
print(summary_stats_sim)

# Write out results
outFile <- paste0(inDir,"/data/regression_models_univariate_statistics_simulated_",version,".csv")
write.csv(summary_stats_sim,outFile, row.names = F)

# Test CCDC GEDI interpolation against other predictors

data_filtered <- df_sim  %>%
  filter(!is.na(spatial_psuedo_r2), !is.na(predictor))
gedi_vals <- data_filtered %>%
  filter(predictor == "deltaGediCCDC") %>%
  pull(spatial_psuedo_r2)
results <- data_filtered %>%
  filter(predictor != "deltaGediCCDC") %>%
  group_by(predictor) %>%
  summarize(
    p_value = wilcox.test(gedi_vals, spatial_psuedo_r2, alternative = "less")$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "holm"))

# View results
print(results)


sar_sim <- df_sim %>% filter((predictor == "deltaSentinel1NDBIVV" | predictor == "deltaSentinel1NDBIVH") & (response == "delta_strct_FHD" | response == "delta_strct_RH_98")) %>% select(predictor, response, spatial_psuedo_r2)
ccdc_sim <- df_sim %>% filter(predictor == "deltaGediCCDC") %>% select(predictor, response, ols_r2, spatial_psuedo_r2)

biomass_sim <- df_sim %>% filter(response == "delta_strct_tAGBD") %>% select(predictor, response, spatial_psuedo_r2)
cover_sim <- df_sim %>% filter(response == "delta_strct_cover") %>% select(predictor, response, spatial_psuedo_r2)

###############################################################################
## Simulated vs on-orbit regression comparison table
###############################################################################

## Table for pseudo-R2

summary_stats_sim <- df_sim  %>%
  group_by(predictor) %>%
  summarize(
    mean_r2 = mean(spatial_psuedo_r2, na.rm = TRUE),
    sd_r2 = sd(spatial_psuedo_r2, na.rm = TRUE),
    min_r2 = min(spatial_psuedo_r2, na.rm = TRUE),
    max_r2 = max(spatial_psuedo_r2, na.rm = TRUE),
    mean_rmse = mean(spatial_scaled_rmse, na.rm = TRUE),
    sd_rmse = sd(spatial_scaled_rmse, na.rm = TRUE),
    n = sum(!is.na(spatial_psuedo_r2)),
    .groups = "drop"
  )

# Compute standard error for both datasets
summary_stats_sim_r2 <- summary_stats_sim %>%
  mutate(se_sim = sd_r2 / sqrt(n)) %>%
  select(predictor, mean_r2_sim = mean_r2, se_sim)

summary_stats_orb <- df_orb %>%
  filter(response %in% unique(df_sim$response)) %>% # filter to those variables in simulated
  group_by(predictor) %>%
  summarize(
    mean_r2 = mean(spatial_psuedo_r2, na.rm = TRUE),
    sd_r2 = sd(spatial_psuedo_r2, na.rm = TRUE),
    min_r2 = min(spatial_psuedo_r2, na.rm = TRUE),
    max_r2 = max(spatial_psuedo_r2, na.rm = TRUE),
    mean_rmse = mean(spatial_scaled_rmse, na.rm = TRUE),
    sd_rmse = sd(spatial_scaled_rmse, na.rm = TRUE),
    n = sum(!is.na(spatial_psuedo_r2)),
    .groups = "drop"
  )

summary_stats_orb_r2 <- summary_stats_orb %>%
  mutate(se_orb = sd_r2 / sqrt(n)) %>%
  select(predictor, mean_r2_orb = mean_r2, se_orb)

# Join and compute delta R?
summary_combined <- summary_stats_sim_r2 %>%
  inner_join(summary_stats_orb_r2, by = "predictor") %>%
  mutate(
    delta_r2 = mean_r2_sim - mean_r2_orb,
    sim_label = sprintf("%.2f ± %.2f", mean_r2_sim, se_sim),
    orb_label = sprintf("%.2f ± %.2f", mean_r2_orb, se_orb),
    delta_label = sprintf("+%.2f", delta_r2)
  )

# Desired order of predictors
predictor_order <- c(
  "deltaLandsat89NDVI", "deltaLandsat89NDMI",
  "deltaSentinel2NDVI", "deltaSentinel2NDMI",
  "deltaPalsar2NDBIHH", "deltaPalsar2NDBIHV",
  "deltaSentinel1NDBIVV", "deltaSentinel1NDBIVH",
  "deltaGediCCDC"
)

# Apply ordering and clean labels
final_table_r2 <- summary_combined %>%
  mutate(predictor = factor(predictor, levels = predictor_order)) %>%
  arrange(predictor) %>%
  select(
    Predictor = predictor,
    `Simulated R2 ± SE` = sim_label,
    `On-orbit R2 ± SE` = orb_label,
    `R2 (Sim - Orbit)` = delta_label
  )

# View table
print(final_table_r2, n = Inf)
write.csv(final_table_r2, paste0(inDir,"/data/regression_models_univariate_r2_260106.csv"))

## Table for RMSE

# Compute standard error for both datasets
summary_stats_sim_rmse <- summary_stats_sim %>%
  mutate(se_sim = sd_rmse / sqrt(n)) %>%
  select(predictor, mean_rmse_sim = mean_rmse, se_sim)

summary_stats_orb_rmse <- summary_stats_orb %>%
  mutate(se_orb = sd_rmse / sqrt(n)) %>%
  select(predictor, mean_rmse_orb = mean_rmse, se_orb)

# Join and compute delta R²
summary_combined <- summary_stats_sim_rmse %>%
  inner_join(summary_stats_orb_rmse, by = "predictor") %>%
  mutate(
    delta_rmse = mean_rmse_sim - mean_rmse_orb,
    sim_label = sprintf("%.2f ± %.2f", mean_rmse_sim, se_sim),
    orb_label = sprintf("%.2f ± %.2f", mean_rmse_orb, se_orb),
    delta_label = sprintf("%.2f", delta_rmse)
  )

# Apply ordering and clean labels
final_table_rmse <- summary_combined %>%
  mutate(predictor = factor(predictor, levels = predictor_order)) %>%
  arrange(predictor) %>%
  select(
    Predictor = predictor,
    `Simulated RMSE ± SE` = sim_label,
    `On-orbit RMSE ± SE` = orb_label,
    `RMSE (Sim - Orbit)` = delta_label
  )

# View table
print(final_table_rmse, n = Inf)
write.csv(final_table_rmse, paste0(inDir,"/data/regression_models_univariate_rmse_260106.csv"))

###############################################################################
### Plot of spatial regression plot slopes 
###############################################################################

extract_slopes <- function(coef_table, predictor_var = "predictor", response_var = NA, class_var = NA) {
  
  # Clean term names
  coef_table$term <- as.character(coef_table$term)

  # Filter to class and response
  coef_table <- coef_table %>% filter(response == response_var & class == class_var)
  
  # Base (reference) class: just the main effect of the predictor
  base_row <- coef_table %>% filter(term == predictor_var)
  if (nrow(base_row) == 0) stop("Predictor term not found in table.")
  
  base_slope <- base_row$Estimate
  base_se <- base_row$Std..Error
  
  # Extract all interaction terms for the predictor
  interaction_terms <- grep(paste0("^", predictor_var, ":"), coef_table$term, value = TRUE)
  
  # Parse class levels from interaction terms (after colon)
  class_levels <- gsub(paste0("^", predictor_var, ":"), "", interaction_terms)
  
  # Reference class is what's *not* in the interaction terms
  all_classes <- unique(gsub(paste0("^", predictor_var, ":"), "", coef_table$term[grep(paste0("^", predictor_var, ":class"), coef_table$term)]))
  reference_class <- "Reference"
  
  # Build result list
  results <- list(
    data.frame(
      term = predictor_var,
      slope = base_slope,
      se = base_se,
      class = reference_class,
      pval = base_row$Pr...z..,
      signif = base_row$signif,
      stringsAsFactors = FALSE
    )
  )
  
  # Process each interaction
  for (t in interaction_terms) {
    class_name <- sub(paste0("^", predictor_var, ":"), "", t)
    int_row <- coef_table %>% filter(term == t)
    
    slope <- base_slope + int_row$Estimate
    slope_se <- sqrt(base_se^2 + int_row$Std..Error^2)  # conservative approx
    
    results[[class_name]] <- data.frame(
      term = t,
      slope = slope,
      se = slope_se,
      class = class_name,
      pval = base_row$Pr...z..,
      signif = int_row$signif,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine and add metadata
  final <- bind_rows(results)
  final$response <- response_var
  final$predictor <- predictor_var
  final$class <- class_var
  return(final)
}

# version date
version <- "250722"

# Pvalue for significance
pvalue = 0.01

plot_list = list()

### Plot for Pre-treatment Structure

# Levels 
response_levels <- c("delta_strct_FHD","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD","delta_strct_tPAI",
                     "delta_strct_mPAI_b10","delta_strct_prop_int_below_10m",
                     "delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75","delta_strct_RH_98")

# get regression results tables 
df <- read.csv(paste0(inDir,"/data/univariate_spatialreg_PreTreatment_footprint_",version,".csv"))
df_coeff <- read.csv(paste0(inDir,"/data/univariate_spatialreg_PreTreatment_footprint_zscore_coeff_",version,".csv"))

# filter
df <- df %>% filter(response %in% response_levels)
df_coeff <- df_coeff %>% filter(response %in% response_levels)

# Plot of slopes treatment impact
df1 <- df %>% select(spatial_slope,spatial_pvalslope, response) %>%
              rename(slope = spatial_slope,
                     pval = spatial_pvalslope) %>%
              mutate(term = "baseline") %>%
              select(term, slope, pval, response)
df2 <- data.frame()
for (r in response_levels){
  d <- extract_slopes(df_coeff, predictor_var = "predictor", response_var = r, class_var = "TreatmentSeverity")
  df2 <- rbind(df2,d)
}
df2 <- df2 %>% select(-c(predictor,se, class,signif))
df3 <- rbind(df1,df2)
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")

df3$modelFactor <- factor(df3$term, 
                          levels = c("baseline","predictor:class1Low","predictor:class1Moderate","predictor","predictor:class1Uncertain"),
                          labels = c("Baseline","Low","Moderate","High","Uncertain"))
df3$responseFactor <- factor(df3$response, 
                             levels = rev(response_levels))
df3 <- df3 %>%
  mutate(size_group = ifelse(modelFactor == "Baseline", "Baseline", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:4])

p <- ggplot(df3,aes(x=slope,y=responseFactor,shape = significant,color = modelFactor, size = size_group)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +
  theme_bw() +
  scale_size_manual(values = c("Baseline" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = response_labels_rev) +
  labs(x = expression("Effect on " * Delta * "GEDI On-orbit"),
       title = "Pre-treatment Structure") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

plot_list[[1]] = p
#dev.new();print(p)

# Plot of slopes for 2018 forest cover
df1 <- df %>% select(spatial_slope,spatial_pvalslope, response) %>%
  rename(slope = spatial_slope,
         pval = spatial_pvalslope) %>%
  mutate(term = "baseline") %>%
  select(term, slope, pval, response)
df2 <- data.frame()
for (r in response_levels){
  d <- extract_slopes(df_coeff, predictor_var = "predictor", response_var = r, class_var = "forestCover")
  df2 <- rbind(df2,d)
}
df2 <- df2 %>% select(-c(predictor,se, class,signif))
df3 <- rbind(df1,df2)
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")

df3$modelFactor <- factor(df3$term, 
                          levels = c("baseline","predictor:class3Forest cover 25-95%","predictor"),
                          labels = c("Baseline","Forest Cover 25-95%","Forest Cover >= 95%"))
df3$responseFactor <- factor(df3$response, 
                             levels = rev(response_levels))
df3 <- df3 %>%
  mutate(size_group = ifelse(modelFactor == "Baseline", "Baseline", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:4])

p <- ggplot(df3,aes(x=slope,y=responseFactor,shape = significant,color = modelFactor, size = size_group)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +
  theme_bw() +
  scale_size_manual(values = c("Baseline" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = response_labels_rev) +
  labs(x = expression("Effect on " * Delta * "GEDI On-orbit"),
       title = "Pre-treatment Structure") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

plot_list[[2]] = p
#dev.new();print(p)

### Plot for fusion GEDI Structure

# Levels 
response_levels <- c("delta_strct_FHD","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD",
                     "delta_strct_tPAI","delta_strct_mPAI_b10","delta_strct_prop_int_below_10m","delta_strct_RH_98")
response_labels_ccdc_rev <- expression("RH98", RE["<10m"], mPAI["0-10m"], 
                                       "tPAI", "tAGBD", "Cover", "aVDR", "nmode", "FHD")

# get regression results tables 
df <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaGediCCDC_footprint_zscore_",version,".csv"))
df_coeff <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaGediCCDC_footprint_zscore_coeff_",version,".csv"))

# filter
df <- df %>% filter(response %in% response_levels)
df_coeff <- df_coeff %>% filter(response %in% response_levels)

# Plot of slopes treatment impact
df1 <- df %>% select(spatial_slope,spatial_pvalslope, response) %>%
  rename(slope = spatial_slope,
         pval = spatial_pvalslope) %>%
  mutate(term = "baseline") %>%
  select(term, slope, pval, response)
df2 <- data.frame()
for (r in response_levels){
  d <- extract_slopes(df_coeff, predictor_var = "predictor", response_var = r, class_var = "TreatmentSeverity")
  df2 <- rbind(df2,d)
}
df2 <- df2 %>% select(-c(predictor,se, class,signif))
df3 <- rbind(df1,df2)
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")

df3$modelFactor <- factor(df3$term, 
                          levels = c("baseline","predictor:class1Low","predictor:class1Moderate","predictor","predictor:class1Uncertain"),
                          labels = c("Baseline","Low","Moderate","High","Uncertain"))
df3$responseFactor <- factor(df3$response, 
                         levels = rev(response_levels))
df3 <- df3 %>%
  mutate(size_group = ifelse(modelFactor == "Baseline", "Baseline", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:4])

p <- ggplot(df3,aes(x=slope,y=responseFactor,shape = significant,color = modelFactor, size = size_group)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +
  theme_bw() +
  scale_size_manual(values = c("Baseline" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = response_labels_ccdc_rev) +
  labs(x = expression("Effect on " * Delta * "GEDI On-orbit"),
       title = "\u0394GEDI Fusion") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

plot_list[[3]] = p
#dev.new();print(p)

# Plot of slopes 2018 forest cover
df1 <- df %>% select(spatial_slope,spatial_pvalslope, response) %>%
  rename(slope = spatial_slope,
         pval = spatial_pvalslope) %>%
  mutate(term = "baseline") %>%
  select(term, slope, pval, response)
df2 <- data.frame()
for (r in response_levels){
  d <- extract_slopes(df_coeff, predictor_var = "predictor", response_var = r, class_var = "forestCover")
  df2 <- rbind(df2,d)
}
df2 <- df2 %>% select(-c(predictor,se, class,signif))
df3 <- rbind(df1,df2)
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")

df3$modelFactor <- factor(df3$term, 
                          levels = c("baseline","predictor:class3Forest cover 25-95%","predictor"),
                          labels = c("Baseline","Forest Cover 25-95%","Forest Cover >= 95%"))
df3$responseFactor <- factor(df3$response, 
                             levels = rev(response_levels))
df3 <- df3 %>%
  mutate(size_group = ifelse(modelFactor == "Baseline", "Baseline", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:4])

p <- ggplot(df3,aes(x=slope,y=responseFactor,shape = significant,color = modelFactor, size = size_group)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +
  theme_bw() +
  scale_size_manual(values = c("Baseline" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = response_labels_ccdc_rev) +
  labs(x = expression("Effect on " * Delta * "GEDI On-orbit"),
       title = "\u0394GEDI Fusion") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

plot_list[[4]] = p
#dev.new();print(p)

### Plot for Sentinel-2 NDVI

# Levels 
response_levels <- c("delta_strct_FHD","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD","delta_strct_tPAI",
                     "delta_strct_mPAI_b10","delta_strct_prop_int_below_10m",
                     "delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75","delta_strct_RH_98")

# get regression results tables 
df <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel2NDVI_footprint_zscore_",version,".csv"))
df_coeff <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel2NDVI_footprint_zscore_coeff_",version,".csv"))

# filter
df <- df %>% filter(response %in% response_levels)
df_coeff <- df_coeff %>% filter(response %in% response_levels)

# Plot of slopes treatment impact
df1 <- df %>% select(spatial_slope,spatial_pvalslope, response) %>%
  rename(slope = spatial_slope,
         pval = spatial_pvalslope) %>%
  mutate(term = "baseline") %>%
  select(term, slope, pval, response)
df2 <- data.frame()
for (r in response_levels){
  d <- extract_slopes(df_coeff, predictor_var = "predictor", response_var = r, class_var = "TreatmentSeverity")
  df2 <- rbind(df2,d)
}
df2 <- df2 %>% select(-c(predictor,se, class,signif))
df3 <- rbind(df1,df2)
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")

df3$modelFactor <- factor(df3$term, 
                          levels = c("baseline","predictor:class1Low","predictor:class1Moderate","predictor","predictor:class1Uncertain"),
                          labels = c("Baseline","Low","Moderate","High","Uncertain"))
df3$responseFactor <- factor(df3$response, 
                             levels = rev(response_levels))
df3 <- df3 %>%
  mutate(size_group = ifelse(modelFactor == "Baseline", "Baseline", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:4])

p <- ggplot(df3,aes(x=slope,y=responseFactor,shape = significant,color = modelFactor, size = size_group)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +
  theme_bw() +
  scale_size_manual(values = c("Baseline" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = response_labels_rev) +
  labs(x = expression("Effect on " * Delta * "GEDI On-orbit"),
       title = "Sentinel-2 \u0394NDVI") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

plot_list[[5]] = p
#dev.new();print(p)

# Plot of slopes for 2018 forest cover
df1 <- df %>% select(spatial_slope,spatial_pvalslope, response) %>%
  rename(slope = spatial_slope,
         pval = spatial_pvalslope) %>%
  mutate(term = "baseline") %>%
  select(term, slope, pval, response)
df2 <- data.frame()
for (r in response_levels){
  d <- extract_slopes(df_coeff, predictor_var = "predictor", response_var = r, class_var = "forestCover")
  df2 <- rbind(df2,d)
}
df2 <- df2 %>% select(-c(predictor,se, class,signif))
df3 <- rbind(df1,df2)
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")

df3$modelFactor <- factor(df3$term, 
                          levels = c("baseline","predictor:class3Forest cover 25-95%","predictor"),
                          labels = c("Baseline","Forest Cover 25-95%","Forest Cover >= 95%"))
df3$responseFactor <- factor(df3$response, 
                             levels = rev(response_levels))
df3 <- df3 %>%
  mutate(size_group = ifelse(modelFactor == "Baseline", "Baseline", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:4])

p <- ggplot(df3,aes(x=slope,y=responseFactor,shape = significant,color = modelFactor, size = size_group)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +
  theme_bw() +
  scale_size_manual(values = c("Baseline" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = response_labels_rev) +
  labs(x = expression("Effect on " * Delta * "GEDI On-orbit"),
       title = "Sentinel-2 \u0394NDVI") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

plot_list[[6]] = p
#dev.new();print(p)

### Plot for Sentinel-2 NDMI

# Levels 
response_levels <- c("delta_strct_FHD","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD","delta_strct_tPAI",
                     "delta_strct_mPAI_b10","delta_strct_prop_int_below_10m",
                     "delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75","delta_strct_RH_98")

# get regression results tables 
df <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel2NDMI_footprint_zscore_",version,".csv"))
df_coeff <- read.csv(paste0(inDir,"/data/univariate_spatialreg_deltaSentinel2NDMI_footprint_zscore_coeff_",version,".csv"))

# filter
df <- df %>% filter(response %in% response_levels)
df_coeff <- df_coeff %>% filter(response %in% response_levels)

# Plot of slopes treatment impact
df1 <- df %>% select(spatial_slope,spatial_pvalslope, response) %>%
  rename(slope = spatial_slope,
         pval = spatial_pvalslope) %>%
  mutate(term = "baseline") %>%
  select(term, slope, pval, response)
df2 <- data.frame()
for (r in response_levels){
  d <- extract_slopes(df_coeff, predictor_var = "predictor", response_var = r, class_var = "TreatmentSeverity")
  df2 <- rbind(df2,d)
}
df2 <- df2 %>% select(-c(predictor,se, class,signif))
df3 <- rbind(df1,df2)
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")

df3$modelFactor <- factor(df3$term, 
                          levels = c("baseline","predictor:class1Low","predictor:class1Moderate","predictor","predictor:class1Uncertain"),
                          labels = c("Baseline","Low","Moderate","High","Uncertain"))
df3$responseFactor <- factor(df3$response, 
                             levels = rev(response_levels))
df3 <- df3 %>%
  mutate(size_group = ifelse(modelFactor == "Baseline", "Baseline", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:4])

p <- ggplot(df3,aes(x=slope,y=responseFactor,shape = significant,color = modelFactor, size = size_group)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +
  theme_bw() +
  scale_size_manual(values = c("Baseline" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = response_labels_rev) +
  labs(x = expression("Effect on " * Delta * "GEDI On-orbit"),
       title = "Sentinel-2 \u0394NDMI") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

plot_list[[7]] = p
#dev.new();print(p)

# Plot of slopes for 2018 forest cover
df1 <- df %>% select(spatial_slope,spatial_pvalslope, response) %>%
  rename(slope = spatial_slope,
         pval = spatial_pvalslope) %>%
  mutate(term = "baseline") %>%
  select(term, slope, pval, response)
df2 <- data.frame()
for (r in response_levels){
  d <- extract_slopes(df_coeff, predictor_var = "predictor", response_var = r, class_var = "forestCover")
  df2 <- rbind(df2,d)
}
df2 <- df2 %>% select(-c(predictor,se, class,signif))
df3 <- rbind(df1,df2)
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")

df3$modelFactor <- factor(df3$term, 
                          levels = c("baseline","predictor:class3Forest cover 25-95%","predictor"),
                          labels = c("Baseline","Forest Cover 25-95%","Forest Cover >= 95%"))
df3$responseFactor <- factor(df3$response, 
                             levels = rev(response_levels))
df3 <- df3 %>%
  mutate(size_group = ifelse(modelFactor == "Baseline", "Baseline", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:4])

p <- ggplot(df3,aes(x=slope,y=responseFactor,shape = significant,color = modelFactor, size = size_group)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2)) +
  theme_bw() +
  scale_size_manual(values = c("Baseline" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = response_labels_rev) +
  labs(x = expression("Effect on " * Delta * "GEDI On-orbit"),
       title = "Sentinel-2 \u0394NDMI") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

plot_list[[8]] = p
#dev.new();print(p)

combined3 <- ggarrange(plotlist = c(plot_list[1],plot_list[5],plot_list[3]),
                       labels = c("A", "B", "C"),
                       common.legend = TRUE,
                       legend = "bottom",
                       ncol = 3, nrow = 1,
                       align = "hv")

dev.new();print(combined3)

outFile <- paste0(inDir,"/figures/regression_models_multivariate_statistics_treatment_intensity_zscore_",version,".png")
ggsave(outFile, width = 8, height = 4, units = "in", dpi=600)

combined4 <- ggarrange(plotlist = c(plot_list[2],plot_list[6],plot_list[4]),
                       labels = c("A", "B", "C"),
                       common.legend = TRUE,
                       legend = "bottom",
                       ncol = 3, nrow = 1,
                       align = "hv")

dev.new();print(combined4)

outFile <- paste0(inDir,"/figures/regression_models_multivariate_statistics_nlcd2018_forest_cover_zscore_",version,".png")
ggsave(outFile, width = 8, height = 4, units = "in", dpi=600)


###############################################################################
### Scatterplot of on-orbit vs. simulated GEDI
###############################################################################

# Levels and labels
response_levels <- c("strct_FHD","strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                     "strct_mPAI_b10","strct_prop_int_below_10m","strct_RH_10",
                     "strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98")
response_levels_sim <- c("strct_FHD","strct_VDR","strct_cover","strct_tAGBD", "strct_tPAI",
                         "strct_mPAI_b10","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98")
response_intersect <- intersect(response_levels, response_levels_sim)
response_not_common <- union(setdiff(response_levels, response_levels_sim), setdiff(response_levels_sim,response_levels))

dsim <- read.csv(paste0(inDir,"/data/master_simulated_gedi_treatment_difference_260106.csv"))

dsim <- dsim[, !grepl("elev|gauss", names(dsim))] # remove elevation and Gaussian-ground metrics

#dsim <- dsim[, !grepl("elev|strct", names(dsim))] # remove elevation and structure
#names(dsim) <- gsub("gauss", "strct", names(dsim))

d <- dsim %>% filter(Statistic == "Mean")
summary_df <- d %>%
  group_by(Type, preYear) %>%
  summarise(total_pre_n = sum(pre_n, na.rm = TRUE),
            count = n()) %>%
  arrange(Type, preYear)

print(summary_df)

dorb <- read.csv(paste0(inDir,"/data/master_gedi_treatment_difference_250710.csv"))

commonIDs <- intersect(dsim$TreatmentID,dorb$TreatmentID)
dsim_treatment <- dsim %>% filter(TreatmentID %in% commonIDs & Statistic == "Mean" & Type == "Treatment") %>%
  mutate(GEDI_type = "Simulated", UniqueID = TreatmentID) %>%
  select(-matches(paste(response_not_common, collapse = "|"))) %>%
  select(-c("TreatmentID","ControlID","preYear","postYear"))
dorb_treatment <- dorb %>% filter(TreatmentID %in% commonIDs & Statistic == "Mean" & Type == "Treatment") %>%
  mutate(GEDI_type = "Onorbit", UniqueID = TreatmentID) %>%
  select(-matches(paste(response_not_common, collapse = "|"))) %>%
  select(-c("Source","Notes","TreatmentID","ControlID"))

commonIDs <- intersect(dsim$ControlID,dorb$ControlID)
dsim_control <- dsim %>% filter(ControlID %in% commonIDs & Statistic == "Mean" & Type == "Control") %>%
  mutate(GEDI_type = "Simulated", UniqueID = ControlID) %>%
  select(-matches(paste(response_not_common, collapse = "|"))) %>%
  select(-c("TreatmentID","ControlID","preYear","postYear"))
dorb_control <- dorb %>% filter(ControlID %in% commonIDs & Statistic == "Mean" & Type == "Control") %>%
  mutate(GEDI_type = "Onorbit", UniqueID = ControlID) %>%
  select(-matches(paste(response_not_common, collapse = "|"))) %>%
  select(-c("Source","Notes","TreatmentID","ControlID"))

df <- rbind(dsim_treatment, dsim_control, dorb_treatment, dorb_control)

strct_cols <- grep("strct", names(df), value = TRUE)
id_cols <- c("GEDI_type", "Type", "UniqueID")
df_strct <- df %>%
  select(all_of(c(strct_cols, id_cols)))

df_long <- df_strct %>%
  pivot_longer(cols = -all_of(id_cols), names_to = "metric", values_to = "value") %>%
  extract(metric, into = c("stage", "metric_base"), regex = "^(pre|post|delta)_(.+)$") %>%
  filter(!is.na(stage))

df_wide <- df_long %>%
  pivot_wider(names_from = GEDI_type, values_from = value) %>%
  drop_na(Simulated, Onorbit)

axis_ranges <- df_wide %>%
  pivot_longer(c(Simulated, Onorbit), names_to = "source", values_to = "value") %>%
  group_by(metric_base, stage) %>%
  summarise(min_val = min(value, na.rm = TRUE),
            max_val = max(value, na.rm = TRUE), .groups = "drop")

df_plot <- df_wide %>%
  left_join(axis_ranges, by = c("metric_base", "stage"))

coords <- df %>%
  group_by(UniqueID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(UniqueID, Mean_x_teale, Mean_y_teale)
df_plot <- left_join(df_plot, coords, by = "UniqueID")

stats_spdep <- df_plot %>%
  group_by(metric_base, stage) %>%
  group_modify(~ {
    dat <- .x
    # build spatial weights on jittered x_teale/y_teale
    coords <- cbind(dat$Mean_x_teale, dat$Mean_y_teale)
    knn_nb <- knearneigh(coords, k = 10)        
    nb     <- knn2nb(knn_nb, row.names = dat$shot_id)
    lw     <- nb2listw(nb, style = "W", zero.policy = TRUE)
    
    # fit a spatial error (SEM) model
    sem <- errorsarlm(Onorbit ~ Simulated,
                      data        = dat,
                      listw       = lw,
                      method="LU",
                      zero.policy = TRUE)
    
    model_summary <- summary(sem,Nagelkerke = TRUE)
    
    # fit OLS model
    ols_model <- lm(Onorbit ~ Simulated, data = dat)
    ols_summary <- summary(ols_model)
    
    # 3) extract metrics
    tibble(
      spatial_int =  model_summary$Coef[1,1],
      spatial_slope =  model_summary$Coef[2,1],
      spatial_stderrint = model_summary$Coef[1,2],
      spatial_stderrslope = model_summary$Coef[2,2],
      spatial_pvalint = model_summary$Coef[1,4],
      spatial_pvalslope = model_summary$Coef[2,4],
      spatial_psuedo_r2 = model_summary$NK,
      ols_int = ols_summary$coefficients[1,1],
      ols_slope = ols_summary$coefficients[2,1],
      ols_pvalslope = ols_summary$coefficients[1,4],
      ols_pvalint = ols_summary$coefficients[2,4],
      ols_r2 = ols_summary$adj.r.squared
    )
  }) %>%
  ungroup()

stats_df <- stats_spdep %>% 
  mutate(
    sem_stars = case_when(
      spatial_pvalslope < 0.001 ~ "***",
      spatial_pvalslope < 0.01  ~ "**",
      spatial_pvalslope < 0.05  ~ "*",
      TRUE ~ ""
    ),
    ols_stars = case_when(
      ols_pvalslope < 0.001 ~ "***",
      ols_pvalslope < 0.01  ~ "**",
      ols_pvalslope < 0.05  ~ "*",
      TRUE ~ ""
    ),
    sem_label = sprintf(
      "Pseudo-R^2 == %.2f*'%s'",
      spatial_psuedo_r2,
      sem_stars
    ),
    ols_label = sprintf(
      "OLS~R^2 == %.2f*'%s'",
      ols_r2,
      ols_stars
    ),
    combo = paste(sem_label, ols_label, sep = "*'\n'*")
  )

df_plot$stage <- factor(df_plot$stage, levels = c("pre", "post", "delta"))

stats_df$stage <- factor(stats_df$stage, levels = c("pre", "post", "delta"))

stage_labels <- c(
  "pre" = "Year 1",
  "post" = "Year 2",
  "delta" = "Change"
)

response_levels <- c(
  "strct_FHD"       = "FHD",
  "strct_VDR"       = "aVDR",
  "strct_cover"     = "Cover",
  "strct_tAGBD"     = "tAGBD",
  "strct_tPAI"     = "tPAI",
  "strct_mPAI_b10"  = "mPAI['0-10m']",
  "strct_RH_25"     = "RH25",
  "strct_RH_50"     = "RH50",
  "strct_RH_75"     = "RH75",
  "strct_RH_98"     = "RH98"
)
df_plot$metric_base <- factor(df_plot$metric_base, levels = names(response_levels))
stats_df$metric_base <- factor(stats_df$metric_base, levels = names(response_levels))

p1 <- ggplot(df_plot, aes(x = Simulated, y = Onorbit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_abline(
    data      = stats_df,
    aes(slope = spatial_slope, intercept = spatial_int),
    linetype  = "solid",
    color     = "black"
  ) +
  geom_point(aes(color = Type), alpha = 0.7) +
  geom_blank(aes(x = min_val, y = min_val)) +
  geom_blank(aes(x = max_val, y = max_val)) +
  geom_text(
    data = stats_df,
    aes(x = Inf, y = -Inf, label = sem_label),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = -0.5, size = 3,
    parse = TRUE
  ) +
  facet_grid2(
    metric_base ~ stage,
    scales = "free",
    independent = "x",
    labeller = labeller(
      metric_base = as_labeller(response_levels, label_parsed),
      stage = as_labeller(stage_labels)
    )
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  labs(
    x = "Simulated GEDI",
    y = "On-orbit GEDI"
  )

dev.new();print(p1)

outFile <- paste0(inDir,"/figures/scatterplot_onorbit_simulated_260106.png")
#outFile <- paste0(inDir,"/figures/scatterplot_onorbit_simulated_gaussian_260101.png")
ggsave(outFile, width = 8, height = 9, units = "in", dpi=600)


## Version with separated Control and Treatment facets

df_plot <- df_plot %>%
  mutate(stage = case_when(
    stage == "delta" & Type == "Treatment" ~ "Change (Treatment)",
    stage == "delta" & Type == "Control" ~ "Change (Control)",
    TRUE ~ stage
  ))

df_plot$stage <- factor(df_plot$stage, levels = c("pre", "post", "Change (Control)", "Change (Treatment)"))
stage_labels <- c(
  "pre" = "Year 1",
  "post" = "Year 2",
  "Change (Control)" = "Change (Control)",
  "Change (Treatment)" = "Change (Treatment)"
)

axis_limits <- df_plot %>%
  pivot_longer(cols = c(Simulated, Onorbit), names_to = "source", values_to = "value") %>%
  group_by(metric_base, stage) %>%
  summarise(
    axis_min = min(value, na.rm = TRUE),
    axis_max = max(value, na.rm = TRUE),
    .groups = "drop"
  )
df_plot <- df_plot %>%
  left_join(axis_limits, by = c("metric_base", "stage"))

r2_labels <- df_plot %>%
  group_by(metric_base, stage) %>%
  summarise({
    mod <- lm(Onorbit ~ Simulated)
    r2 <- summary(mod)$r.squared
    pval <- glance(mod)$p.value
    stars <- case_when(
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      TRUE ~ ""
    )
    label <- paste0("R2 = ", round(r2, 2), " ", stars)
    tibble(label = label)
  }, .groups = "drop")

p2 <- ggplot(df_plot, aes(x = Simulated, y = Onorbit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_point(aes(color = Type), alpha = 0.7) +
  geom_blank(aes(x = axis_min, y = axis_min)) +
  geom_blank(aes(x = axis_max, y = axis_max)) +
  geom_text(data = r2_labels, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3.5) +
  facet_grid2(
    metric_base ~ stage,
    scales = "free",
    independent = "x",
    labeller = labeller(
      metric_base = as_labeller(response_levels, label_parsed),
      stage = stage_labels
    )) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  labs(
    x = "Simulated GEDI",
    y = "On-orbit GEDI"
  )


dev.new();print(p2)

outFile <- paste0(inDir,"/figures/scatterplot_onorbit_simulated_separate_change_260106.png")
#outFile <- paste0(inDir,"/figures/scatterplot_onorbit_simulated_separate_change_gaussian_260101.png")
ggsave(outFile, width = 8, height = 8, units = "in", dpi=600)

# summary statistics

summary_by_stage <- stats_df %>%
  group_by(stage) %>%
  summarise(
    # your existing stats
    across(
      spatial_psuedo_r2,
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x,   na.rm = TRUE),
        min  = ~min(.x,  na.rm = TRUE),
        max  = ~max(.x,  na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    min_metric = metric_base[which.min(spatial_psuedo_r2)],
    max_metric = metric_base[which.max(spatial_psuedo_r2)]
  )

print(summary_by_stage)
write.csv(summary_by_stage,paste0(inDir,"/data/simulated_onorbit_r2_summary.csv"), row.names = F)

table(df_plot$Type)/length(unique(df_plot$metric_base))/3

###############################################################################
### Scatterplot of CCDC fusion vs. simulated GEDI
###############################################################################

ground <- "real"
#ground <- "gauss"

winter <- "with_winter"
#winter <- "no_winter"

# NLCD land cover classes to include: 41 = deciduous forest; 42 = evergreen forest, 43 = mixed forest, 52 = shrub/scrub 
landcover_classes <- c(42) # just conifer 
# landcover_classes <- c(41,42,43,52) # all forest and shrub

version <- "260106"

df <- read_csv(paste0(inDir,"/data/simulated_gedi_ccdc_fusion_",version,".csv"))
df$shot_id <- paste0(df$x_teale,"_",df$y_teale)

# Use real or Gaussian ground
if (ground == "real"){ # use real ground
  df <- df[, !grepl("gauss", names(df), ignore.case = TRUE)] # drop Gaussian-fitted metrics
  names(df) <- gsub("^strct", "sim_strct", names(df))
} else { # use gaussian ground
  df <- df[, !grepl("^strct", names(df), ignore.case = TRUE)] # drop real ground metrics
  names(df) <- gsub("gauss", "sim_strct", names(df))
}

# filter out acquisitions with winter
if (winter == "no_winter"){
  df <- df %>% filter(!grepl("winter", season))
}

# filter for desired landcover classes based on majority
df <- df %>% filter(shot_lc_maj1 %in% landcover_classes)

df_long <- df %>%
  pivot_longer(
    cols = matches("^(sim|ccdc)_strct_"), 
    names_to  = c("source", "metric"),
    names_pattern = "(sim|ccdc)_strct_(.*)",
    values_to = "value"
  ) %>%
  filter(metric %in% c("RH_98", "mPAI_b10", "cover", "VDR", "FHD","tAGBD","tPAI")) %>%
  pivot_wider(names_from = source, values_from = value)

df_plot <- df_long %>%
  filter(!is.na(sim), !is.na(ccdc)) %>%
  mutate(
    x_jit = jitter(x_teale, amount = 0.001),
    y_jit = jitter(y_teale, amount = 0.001)
  )

# filters to footprints that have both pre and post
df_plot2 <- df_plot %>%
  group_by(shot_id) %>%
  filter(n_distinct(YearALS) > 1) %>%
  mutate(
    stage = case_when(
      YearALS == min(YearALS) ~ "Pre",
      YearALS == max(YearALS) ~ "Post",
      TRUE                   ~ NA_character_
    )
  ) %>%
  ungroup()
df_plot2$shot_id <- 1:nrow(df_plot2)

# Outlier removal helper (robust z using MAD) ---
robust_keep <- function(x, zmax = 4) {
  m <- median(x, na.rm = TRUE)
  s <- mad(x, constant = 1.4826, na.rm = TRUE)  # ~robust SD
  if (is.na(s) || s == 0) return(rep(TRUE, length(x)))
  abs((x - m) / s) <= zmax
}

# Choose threshold (3.5-5 is common; 4 is a decent default)
zmax <- 4

df_plot2 <- df_plot2 %>%
  group_by(metric) %>%
  filter(
    robust_keep(sim, zmax = zmax),
    robust_keep(ccdc, zmax = zmax)
  ) %>%
  ungroup()

# year filter, if desired
#df_plot2 <- df_plot2 %>% filter(YearALS == 2020)

stats_spdep <- df_plot2 %>%
  group_by(metric) %>%
  group_modify(~ {
    dat <- .x
    # build spatial weights on jittered x_teale/y_teale
    coords <- cbind(dat$x_jit, dat$y_jit)
    knn_nb <- knearneigh(coords, k = 50)        
    nb     <- knn2nb(knn_nb, row.names = dat$shot_id)
    lw     <- nb2listw(nb, style = "W", zero.policy = TRUE)
    
    # fit a spatial error (SEM) model
    sem <- errorsarlm(ccdc ~ sim,
                      data        = dat,
                      listw       = lw,
                      method="LU",
                      zero.policy = TRUE)
    
    model_summary <- summary(sem,Nagelkerke = TRUE)
    
    # fit OLS model
    ols_model <- lm(ccdc ~ sim, data = dat)
    ols_summary <- summary(ols_model)
    
    # 3) extract metrics
    tibble(
      spatial_int =  model_summary$Coef[1,1],
      spatial_slope =  model_summary$Coef[2,1],
      spatial_stderrint = model_summary$Coef[1,2],
      spatial_stderrslope = model_summary$Coef[2,2],
      spatial_pvalint = model_summary$Coef[1,4],
      spatial_pvalslope = model_summary$Coef[2,4],
      spatial_psuedo_r2 = model_summary$NK,
      ols_int = ols_summary$coefficients[1,1],
      ols_slope = ols_summary$coefficients[2,1],
      ols_pvalslope = ols_summary$coefficients[1,4],
      ols_pvalint = ols_summary$coefficients[2,4],
      ols_r2 = ols_summary$adj.r.squared
    )
  }) %>%
  ungroup()

stats_df <- stats_spdep %>%
  mutate(
    sem_stars = case_when(
      spatial_pvalslope < 0.001 ~ "***",
      spatial_pvalslope < 0.01  ~ "**",
      spatial_pvalslope < 0.05  ~ "*",
      TRUE ~ ""
    ),
    ols_stars = case_when(
      ols_pvalslope < 0.001 ~ "***",
      ols_pvalslope < 0.01  ~ "**",
      ols_pvalslope < 0.05  ~ "*",
      TRUE ~ ""
    ),
    strip_label = sprintf(
      "Pseudo-R<sup>2</sup> = %.2f%s<br>OLS R<sup>2</sup> = %.2f%s",
      spatial_psuedo_r2, sem_stars,
      ols_r2, ols_stars
    )
  )


save(stats_df , file = paste0(inDir,"/data/models_simulated_ccdc_fusion_pre-post_filtered_",winter,"_",ground,"_",version,".RData"))

metric_labs <- c(
  cover     = "Cover",
  VDR       = "aVDR",
  RH_98     = "RH98",
  mPAI_b10  = "mPAI0-10m",
  FHD  = "FHD",
  tAGBD = "tAGBD",
  tPAI = "tPAI"
)

lims_df <- df_plot2 %>%
  group_by(metric) %>%
  summarise(
    lo = min(c(sim, ccdc), na.rm = TRUE),
    hi = max(c(sim, ccdc), na.rm = TRUE),
    .groups = "drop"
  )

anchor_df <- lims_df %>%
  pivot_longer(c(lo, hi), names_to = "which", values_to = "v") %>%
  transmute(metric, sim = v, ccdc = v)

p_hex <- ggplot(df_plot2, aes(sim, ccdc)) +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(name = "Count") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_abline(data = stats_df, aes(slope = spatial_slope, intercept = spatial_int),
              linetype = "dashed", color = "red") +
  geom_abline(data = stats_df, aes(slope = ols_slope, intercept = ols_int),
              linetype = "dashed", color = "cyan") +
  geom_blank(data = anchor_df) +   # forces x/y limits to match within each facet
  facet_wrap(
    ~ metric,
    scales = "free",
    labeller = labeller(metric = function(m) {
      paste0("<b>", metric_labs[m], "</b><br>",
             stats_df$strip_label[match(m, stats_df$metric)])
    })
  ) +
  labs(x = "Simulated GEDI", y = "Fusion GEDI") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = ggtext::element_markdown(lineheight = 1.1, size = 10),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA)
  )

dev.new();print(p_hex)

outFile <- paste0(inDir,"/figures/scatterplot_simulated_ccdc_fusion_pre-post_filtered_",winter,"_",ground,"_",version,".png")
ggsave(outFile, width = 8, height = 7, units = "in", dpi=600)

# R2 summary
stats_yearly <- stats_df %>% select(metric, ols_r2, spatial_psuedo_r2)

# Percent of ALS footprints by ecoregion
d <- df_plot %>% filter(metric == "RH_98")
table(d$US_L3NAME)/dim(d)[1]

# 
footprints_per_polygon_stats <- d %>%
  group_by(Type, UniqueID) %>%
  summarise(n = n(), .groups = "drop") %>%        # Step 1: count per UniqueID
  group_by(Type) %>%
  summarise(
    mean_count = mean(n),
    sd_count = sd(n),
    .groups = "drop"
  )

###############################################################################
### Scatterplot of CCDC fusion vs. simulated GEDI change
###############################################################################

ground <- "real"
#ground <- "gauss"

winter <- "with_winter"
#winter <- "no_winter"

# NLCD land cover classes to include: 41 = deciduous forest; 42 = evergreen forest, 43 = mixed forest, 52 = shrub/scrub 
landcover_classes <- c(42) # just conifer 
# landcover_classes <- c(41,42,43,52) # all forest and shrub

version <- "260106"

df <- read_csv(paste0(inDir,"/data/simulated_gedi_ccdc_fusion_",version,".csv"))
df$shot_id <- paste0(df$x_teale,"_",df$y_teale)

# Use real or Gaussian ground
if (ground == "real"){ # use real ground
  df <- df[, !grepl("gauss", names(df), ignore.case = TRUE)] # drop Gaussian-fitted metrics
  names(df) <- gsub("^strct", "sim_strct", names(df))
} else { # use gaussian ground
  df <- df[, !grepl("^strct", names(df), ignore.case = TRUE)] # drop real ground metrics
  names(df) <- gsub("gauss", "sim_strct", names(df))
}

# filter out acquisitions with winter
if (winter == "no_winter"){
  df <- df %>% filter(!grepl("winter", season))
}

# filter for desired landcover classes based on majority
df <- df %>% filter(shot_lc_maj1 %in% landcover_classes)
  
df_long <- df %>%
  pivot_longer(
    cols = matches("^(sim|ccdc)_strct_"), 
    names_to  = c("source", "metric"),
    names_pattern = "(sim|ccdc)_strct_(.*)",
    values_to = "value"
  ) %>%
  #filter(metric %in% c("RH_98", "mPAI_b10", "cover", "VDR", "FHD","tAGBD","tPAI")) %>%
  filter(metric %in% c("RH_98", "mPAI_b10", "cover", "VDR", "FHD","tAGBD","tPAI") & Type == "treatment") %>% # just use treatments for change
  pivot_wider(names_from = source, values_from = value) %>% 
  select(-c(ID, UniqueID, row_id, x_teale, y_teale, ccdcName, US_L3NAME))

df_plot <- df_long %>%
  filter(!is.na(sim), !is.na(ccdc))

df_plot2 <- df_plot %>%
  group_by(shot_id) %>%
  filter(n_distinct(YearALS) > 1) %>%
  mutate(
    stage = case_when(
      YearALS == min(YearALS) ~ "Pre",
      YearALS == max(YearALS) ~ "Post",
      TRUE                   ~ NA_character_
    )
  ) %>%
  ungroup()

# Outlier removal helper (robust z using MAD) ---
robust_keep <- function(x, zmax = 4) {
  m <- median(x, na.rm = TRUE)
  s <- mad(x, constant = 1.4826, na.rm = TRUE)  # ~robust SD
  if (is.na(s) || s == 0) return(rep(TRUE, length(x)))
  abs((x - m) / s) <= zmax
}

# Choose threshold (3.5-5 is common; 4 is a decent default)
zmax <- 4

df_plot2 <- df_plot2 %>%
  group_by(metric) %>%
  filter(
    robust_keep(sim, zmax = zmax),
    robust_keep(ccdc, zmax = zmax)
  ) %>%
  ungroup()

# Calculate change
df_change <- df_plot2 %>%
  filter(!is.na(stage)) %>%
  group_by(Type, shot_id, metric) %>%
  filter(n_distinct(stage) == 2) %>%
  summarise(
    change_sim  = sim [stage == "Post"]  - sim [stage == "Pre"],
    change_ccdc = ccdc[stage == "Post"]  - ccdc[stage == "Pre"],
    .groups = "drop"
  ) %>%
  separate(
    shot_id,
    into    = c("x_teale", "y_teale"),
    sep     = "_",
    convert = TRUE  # turns them into numeric
  )

stats_spdep <- df_change %>%
  group_by(metric) %>%
  group_modify(~ {
    dat <- .x
    # build spatial weights on x_teale/y_teale
    coords <- cbind(dat$x_teale, dat$y_teale)
    knn_nb <- knearneigh(coords, k = 50)        
    nb     <- knn2nb(knn_nb, row.names = dat$shot_id)
    lw     <- nb2listw(nb, style = "W", zero.policy = TRUE)
    
    # fit a spatial error (SEM) model
    sem <- errorsarlm(change_ccdc ~ change_sim,
                      data        = dat,
                      listw       = lw,
                      method="LU",
                      zero.policy = TRUE)
    
    model_summary <- summary(sem,Nagelkerke = TRUE)
    
    # fit OLS model
    ols_model <- lm(change_ccdc ~ change_sim, data = dat)
    ols_summary <- summary(ols_model)
    
    # 3) extract metrics
    tibble(
      spatial_int =  model_summary$Coef[1,1],
      spatial_slope =  model_summary$Coef[2,1],
      spatial_stderrint = model_summary$Coef[1,2],
      spatial_stderrslope = model_summary$Coef[2,2],
      spatial_pvalint = model_summary$Coef[1,4],
      spatial_pvalslope = model_summary$Coef[2,4],
      spatial_psuedo_r2 = model_summary$NK,
      ols_int = ols_summary$coefficients[1,1],
      ols_slope = ols_summary$coefficients[2,1],
      ols_pvalslope = ols_summary$coefficients[1,4],
      ols_pvalint = ols_summary$coefficients[2,4],
      ols_r2 = ols_summary$adj.r.squared
    )
  }) %>%
  ungroup()

stats_df <- stats_spdep %>%
  mutate(
    sem_stars = case_when(
      spatial_pvalslope < 0.001 ~ "***",
      spatial_pvalslope < 0.01  ~ "**",
      spatial_pvalslope < 0.05  ~ "*",
      TRUE ~ ""
    ),
    ols_stars = case_when(
      ols_pvalslope < 0.001 ~ "***",
      ols_pvalslope < 0.01  ~ "**",
      ols_pvalslope < 0.05  ~ "*",
      TRUE ~ ""
    ),
    strip_label = sprintf(
      "Pseudo-R<sup>2</sup> = %.2f%s<br>OLS R<sup>2</sup> = %.2f%s",
      spatial_psuedo_r2, sem_stars,
      ols_r2, ols_stars
    )
  )

save(stats_df , file = paste0(inDir,"/data/models_simulated_ccdc_fusion_change_filtered_",winter,"_",ground,"_",version,".RData"))

metric_labs <- c(
  cover     = "Cover",
  VDR       = "aVDR",
  RH_98     = "RH98",
  mPAI_b10  = "mPAI0-10m",
  FHD  = "FHD",
  tAGBD = "tAGBD",
  tPAI = "tPAI"
)

lims_df <- df_change %>%
  group_by(metric) %>%
  summarise(
    lo = min(c(change_sim, change_ccdc), na.rm = TRUE),
    hi = max(c(change_sim, change_ccdc), na.rm = TRUE),
    .groups = "drop"
  )

anchor_df <- lims_df %>%
  pivot_longer(c(lo, hi), names_to = "which", values_to = "v") %>%
  transmute(metric, change_sim = v, change_ccdc = v)

p_hex <- ggplot(df_change, aes(change_sim, change_ccdc)) +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(name = "Count") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_abline(data = stats_df, aes(slope = spatial_slope, intercept = spatial_int),
              linetype = "dashed", color = "red") +
  geom_abline(data = stats_df, aes(slope = ols_slope, intercept = ols_int),
              linetype = "dashed", color = "cyan") +
  geom_blank(data = anchor_df) +   # forces x/y limits to match within each facet
  facet_wrap(
    ~ metric,
    scales = "free",
    labeller = labeller(metric = function(m) {
      paste0("<b>", metric_labs[m], "</b><br>",
             stats_df$strip_label[match(m, stats_df$metric)])
    })
  ) +
  labs(x = "Simulated GEDI", y = "Fusion GEDI") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = ggtext::element_markdown(lineheight = 1.1, size = 10),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA)
  )


dev.new();print(p_hex)

outFile <- paste0(inDir,"/figures/scatterplot_simulated_ccdc_fusion_change_filtered_",winter,"_",ground,"_",version,".png")
ggsave(outFile, width = 8, height = 7, units = "in", dpi=600)

stats_change <- stats_df %>% select(metric, ols_r2, spatial_psuedo_r2)
stats_change 

###############################################################################
## Correlation pre- & post-treatment ALS metrics in Control Areas - Polygon Scale
###############################################################################

version <- "260110"

response_levels_sim <- c("strct_FHD","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                         "strct_mPAI_b10","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98", "strct_elev")
response_labels_sim <- c("FHD","aVDR","Cover","tAGBD","tPAI",expression(mPAI["0-10m"]),"RH25","RH50","RH75","RH98","elev")


df <- read.csv(paste0(inDir,"/data/master_simulated_gedi_treatment_difference_",version,".csv"))

df <- df[, !grepl("gauss", names(df), ignore.case = TRUE)] # drop Gaussian-fitted metrics

# Filter to Control rows
control_df <- df %>% filter(Type == "Control" & Statistic == "Mean")

# Select only columns starting with 'pre_' or 'post_'
metric_df <- control_df %>% select(starts_with("pre_"), starts_with("post_")) %>% 
  select(-c(pre_n, post_n))

# Reorder
metric_df_ordered <- metric_df %>% select(all_of(paste0("pre_",response_levels_sim)),all_of(paste0("post_",response_levels_sim)))

# Compute the correlation matrix
cor_matrix <- cor(metric_df_ordered, use = "pairwise.complete.obs")

# Reshape the correlation matrix for plotting
cor_melted <- melt(cor_matrix)

# Create pre and post labels
labels <- c(paste0("Year1~", response_labels_sim), paste0("Year2~", response_labels_sim))

# Create the plot
p <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  # scale_fill_distiller(palette = "Reds", direction = 1, name = "Pearson r", 
  #                      limits = c(0.5, 1), oob = scales::squish) +
  scale_fill_viridis_c(option = "magma", direction = 1, name = "Pearson r") +
  geom_vline(xintercept = 11.5, color = "blue", linewidth = 1) +
  geom_hline(yintercept = 11.5, color = "blue", linewidth = 1) +
  coord_fixed() +
  scale_x_discrete(labels = parse(text = labels)) +
  scale_y_discrete(labels = parse(text = labels)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Correlation of Simulated Bi-temporal\nGEDI Metrics (Control Areas, Polygon Scale)")

dev.new();print(p)
outFile <- paste0(inDir,"/figures/correlation_matrix_simulated_GEDI_pre-post_control_areas_polygon_",version,".png")
ggsave(outFile, width = 6, height = 6, units = "in", dpi=600,bg = "white")

# Create vector of Year1 and Year2 column names
year1_cols <- paste0("pre_", response_levels_sim )
year2_cols <- paste0("post_", response_levels_sim )

# Initialize vector to store correlations
cor_values <- numeric(length(response_levels_sim))

# Loop to calculate correlations between matching Year1 and Year2 metrics
for (i in seq_along(response_levels_sim)) {
  x <- control_df[[year1_cols[i]]]
  y <- control_df[[year2_cols[i]]]
  cor_values[i] <- cor(x, y, use = "pairwise.complete.obs")
}

# Assemble into data frame for plotting
cor_df <- data.frame(
  Metric = response_levels_sim,
  Pearson_r = cor_values
)

# Relabel for plot
cor_df$Metric <- factor(cor_df$Metric, levels = response_levels_sim, labels = response_labels_sim)

# Plot as bar chart
p <- ggplot(cor_df, aes(x = Metric, y = Pearson_r)) +
  geom_col(fill = "steelblue") +
  ylim(0, 1) +
  geom_text(aes(label = round(Pearson_r, 2)), vjust = -0.5) +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = parse(text = response_labels_sim)) +
  labs(
    title = "Bi-temporal Correlation of GEDI Metrics\nControl Areas, Polygon Scale\n2018–2020 vs 2022",
    y = "Pearson r",
    x = ""
  ) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.new();print(p)
outFile <- paste0(inDir,"/figures/correlation_barchart_simulated_GEDI_pre-post_control_areas_polygon_",version,".png")
ggsave(outFile, width = 5, height = 5, units = "in", dpi=600,bg = "white")

# Model of impact of Year 1 on bi-temporal relationships
df_wide <- control_df %>% select(starts_with("pre"), starts_with("post")) %>% 
  select(-c(pre_n, post_n))

df_balanced <- df_wide %>%
  group_by(preYear) %>%
  slice_sample(n = min(table(df_wide$preYear))) %>%
  ungroup()

results <- response_levels_sim %>%
  map_df(~{
    p <- paste0("pre_",.x)
    r <- paste0("post_",.x)
    
    df_model <- df_wide %>%
      select(Predictor = all_of(p),
             Response = all_of(r),
             Year = preYear) %>%
      filter(!is.na(Predictor) & !is.na(Response) & !is.na(Year))
    
    model <- lm(Response ~ Predictor * Year, data = df_model)
    tidy(model) %>%
      filter(term == "Predictor:Year") %>%  # interaction effect
      mutate(metric = .x)
  })

# View interaction term effects
results %>%
  select(metric, estimate, p.value)
results

# Pre-treatment footprint counts by Type and Year
pre_summary <- df %>%
  filter(preYear %in% c(2018, 2019, 2020)) %>%
  group_by(Type, preYear) %>%
  summarise(total_pre_footprints = sum(pre_n), .groups = "drop")
pre_summary

# Post-treatment footprint counts by Type and Year (only 2022)
post_summary <- df %>%
  filter(postYear == 2022) %>%
  group_by(Type, postYear) %>%
  summarise(total_post_footprints = sum(post_n), .groups = "drop")
post_summary

# Combine into one table
summary_table <- bind_rows(
  pre_summary %>% rename(Year = preYear, Footprints = total_pre_footprints),
  post_summary %>% rename(Year = postYear, Footprints = total_post_footprints)
) %>%
  arrange(Type, Year)

summary_table

###############################################################################
## Correlation pre- & post-treatment ALS metrics with both Gaussian-fitted and real ground
###############################################################################

version <- "260106"

response_levels_sim <- c("strct_FHD","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                         "strct_mPAI_b10","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98", "strct_elev")
response_labels_sim <- c("FHD","aVDR","Cover","tAGBD","tPAI",expression(mPAI["0-10m"]),"RH25","RH50","RH75","RH98","elev")

# Get footprint level data in control areas; filtered version used in the analysis
df <- read_csv(paste0(inDir,"/data/control_simulated_gedi_data_filtered_",version,".csv"))

df$shot_id <- sub(".*_(X[0-9]+_Y[0-9]+)$", "\\1", df$unique_ID)

# filters to footprints that have both pre and post
df1 <- df %>%
  filter(!is.na(YearALS)) %>%  # Remove rows where YearALS is completely NA 
  group_by(shot_id) %>%
  filter(n_distinct(YearALS) > 1) %>%  # Only keep groups with >1 unique non-NA year
  mutate(
    min_year = min(YearALS),
    max_year = max(YearALS),
    stage = case_when(
      YearALS == min_year ~ "Year1",
      YearALS == max_year ~ "Year2",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-min_year, -max_year) %>%  # Clean up helper columns
  ungroup()



# Base metric names (without prefix)
metric_base <- sub("^strct_", "", response_levels_sim)
metrics <- c(paste0("strct_",metric_base),paste0("gauss_",metric_base))
df2 <- df1 %>% select(YearALS,stage, shot_id, all_of(metrics))

df_filter <- df2 %>%
  group_by(shot_id, stage) %>%
  slice(1) %>%   # keeps the first row per group
  ungroup()

df_wide <- df_filter %>%
  pivot_wider(
    names_from = stage,
    values_from = all_of(c(metrics, "YearALS")),
    names_sep = "_"
  )

# Helper: compute Year1 vs Year2 correlations for strct + gauss, for one df_wide subset
calc_cor_by_method <- function(df_sub, metric_base) {
  
  # Build correlation rows for both methods and all metrics
  bind_rows(lapply(c("strct", "gauss"), function(meth) {
    
    # Build the column names we need
    y1 <- paste0(meth, "_", metric_base, "_Year1")
    y2 <- paste0(meth, "_", metric_base, "_Year2")
    
    # Keep only metrics that actually exist in the data
    ok <- y1 %in% names(df_sub) & y2 %in% names(df_sub)
    y1 <- y1[ok]; y2 <- y2[ok]
    mbase_ok <- metric_base[ok]
    
    # Compute correlations metric-by-metric
    r <- mapply(function(a, b) cor(df_sub[[a]], df_sub[[b]], use = "pairwise.complete.obs"),
                y1, y2)
    
    tibble(
      metric = mbase_ok,
      method = meth,
      Pearson_r = as.numeric(r)
    )
  }))
}

cor_2018 <- df_wide %>%
  filter(YearALS_Year1 == 2018) %>%
  calc_cor_by_method(metric_base) %>%
  mutate(panel = "2018 vs 2022")

cor_2019 <- df_wide %>%
  filter(YearALS_Year1 == 2019) %>%
  calc_cor_by_method(metric_base) %>%
  mutate(panel = "2019 vs 2022")

cor_2020 <- df_wide %>%
  filter(YearALS_Year1 == 2020) %>%
  calc_cor_by_method(metric_base) %>%
  mutate(panel = "2020 vs 2022")

cor_all <- df_wide %>%
  calc_cor_by_method(metric_base) %>%
  mutate(panel = "2018–2020 vs 2022")

cor_df <- bind_rows(cor_2018, cor_2019, cor_2020, cor_all) %>%
  mutate(
    metric = factor(metric, levels = metric_base),
    method = factor(method, levels = c("strct", "gauss"))
  )

# Dodge spec reused for bars + text (critical)
pd <- position_dodge2(width = 0.8, preserve = "single")

p_cor <- ggplot(cor_df, aes(x = metric, y = Pearson_r, fill = method)) +
  geom_col(
    position = pd,
    width = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", Pearson_r),
      # stagger labels by method
      y = Pearson_r + ifelse(method == "strct", 0.025, 0.008)
    ),
    position = pd,
    size = 2
  ) +
  facet_wrap(~ panel, ncol = 2) +
  scale_x_discrete(labels = parse(text = response_labels_sim)) +
  scale_fill_manual(
    values = c(
      "strct" = brewer.pal(8, "Set2")[3],
      "gauss" = brewer.pal(8, "Set2")[4]
    ),
    labels = c(
      "strct" = "ALS (real)",
      "gauss" = "Gaussian-fitted"
    )
  ) +
  coord_cartesian(
    ylim = c(0, 1.08),   # headroom for labels
    clip = "off"
  ) +
  labs(
    y = "Pearson r",
    x = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 12)  # extra top margin
  )

dev.new(); print(p_cor)

outFile <- paste0(inDir, "/figures/correlation_barchart_strct_vs_gauss_pre-post_control_areas_footprint_", version, ".png")
ggsave(outFile, p_cor, width = 8, height = 6, units = "in", dpi = 600, bg = "white")


###############################################################################
## Box plot of polygon-level change on-orbit data
###############################################################################

# Levels and labels
response_levels <- c("strct_FHD","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                     "strct_mPAI_b10","strct_prop_int_below_10m",
                     "strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98")
response_labels <- c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                     expression(mPAI["0-10m"]),expression(RE["<10m"]),
                     "RH25","RH50","RH75","RH98")

dorb <- read.csv(paste0(inDir,"/data/master_gedi_treatment_difference_250710.csv"))
dorb <- dorb %>% filter(Statistic == "Mean") 

df <- dorb %>% mutate(TreatmentSeverity = if_else(Type == "Control", "Control", TreatmentSeverity))
strct_cols <- grep("strct", names(df), value = TRUE)
id_cols <- c("TreatmentSeverity")
df_strct <- df %>%
  select(all_of(c(strct_cols, id_cols)))

df_long <- df_strct %>%
  pivot_longer(cols = -all_of(id_cols), names_to = "metric", values_to = "value") %>%
  extract(metric, into = c("stage", "metric_base"), regex = "^(pre|post|delta)_(.+)$") %>%
  filter(!is.na(stage) & stage != "delta" & !metric_base %in% c("strct_RH_10", "strct_H"))
df_long$stage <- factor(df_long$stage, levels = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment"))
df_long$metric <- factor(df_long$metric_base, levels = response_levels, labels = response_labels)
df_long$TreatmentSeverity <- factor(df_long$TreatmentSeverity, 
                                    levels = c("Control", "Low", "Moderate", "High", "Uncertain"),
                                    labels = c("Control", "Low", "Moderate", "High", "Not Classified"))

df_long <- df_long %>%
  mutate(
    value = if_else(metric_base == "strct_cover", value * 100, value)
  )

# outlier clipping
df_long <- df_long %>%
  group_by(metric) %>%
  mutate(
    lower = quantile(value, 0.01, na.rm = TRUE),
    upper = quantile(value, 0.99, na.rm = TRUE),
    value_clipped = pmin(pmax(value, lower), upper)
  ) %>%
  ungroup()


p <- ggplot(df_long, aes(x = TreatmentSeverity, y = value_clipped, fill = stage)) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.7, 
               position = position_dodge(width = 0.75),
               linewidth = 0.3) +
  scale_fill_manual(
    values = c("Pre-treatment" = "#8DA0CB", "Post-treatment" = "#E78AC3")
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4, labeller = label_parsed) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    x = "Treatment Impact",
    y = "Value",
    fill = "Stage"
  )

dev.new();print(p)
outFile <- paste0(inDir,"/figures/metrics_barchart_onorbit_250804.png")
ggsave(outFile, width = 7, height = 5, units = "in", dpi=600,bg = "white")

###############################################################################
## Box plot of polygon-level change simulated data
###############################################################################

# Levels and labels
response_levels <- c("strct_FHD","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                     "strct_mPAI_b10","strct_prop_int_below_10m",
                     "strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98")
response_labels <- c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                     expression(mPAI["0-10m"]),expression(RE["<10m"]),
                     "RH25","RH50","RH75","RH98")

response_levels_sim <- c("strct_FHD","strct_VDR","strct_cover","strct_tAGBD",
                         "strct_mPAI_b10","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98")
response_labels_sim <- c("FHD","aVDR","Cover","tAGBD",expression(mPAI["0-10m"]),"RH25","RH50","RH75","RH98")

dsim <- read.csv(paste0(inDir,"/data/master_simulated_gedi_treatment_difference_260106.csv"))
dsim <- dsim %>% filter(Statistic == "Mean")
dsim <- dsim[, !grepl("gauss|elev", names(dsim), ignore.case = TRUE)] # drop Gaussian-fitted metrics and elevation

df <- dsim %>% mutate(TreatmentSeverity = if_else(Type == "Control", "Control", TreatmentSeverity))
strct_cols <- grep("strct", names(df), value = TRUE)
id_cols <- c("TreatmentSeverity")
df_strct <- df %>%
  select(all_of(c(strct_cols, id_cols)))

df_long <- df_strct %>%
  pivot_longer(cols = -all_of(id_cols), names_to = "metric", values_to = "value") %>%
  extract(metric, into = c("stage", "metric_base"), regex = "^(pre|post|delta)_(.+)$") %>%
  filter(!is.na(stage) & stage != "delta" & !metric_base %in% c("strct_RH_10", "strct_H"))
df_long$stage <- factor(df_long$stage, levels = c("pre", "post"), labels = c("Pre-treatment", "Post-treatment"))
df_long$metric <- factor(df_long$metric_base, levels = response_levels, labels = response_labels)
df_long$TreatmentSeverity <- factor(df_long$TreatmentSeverity, 
                                    levels = c("Control", "Low", "Moderate", "High", "Uncertain"),
                                    labels = c("Control", "Low", "Moderate", "High", "Not Classified"))

df_long <- df_long %>%
  mutate(
    value = if_else(metric_base == "strct_cover", value * 100, value)
  )

# outlier clipping
df_long <- df_long %>%
  group_by(metric) %>%
  mutate(
    lower = quantile(value, 0.01, na.rm = TRUE),
    upper = quantile(value, 0.99, na.rm = TRUE),
    value_clipped = pmin(pmax(value, lower), upper)
  ) %>%
  ungroup()

p <- ggplot(df_long, aes(x = TreatmentSeverity, y = value_clipped, fill = stage)) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.7, 
               position = position_dodge(width = 0.75),
               linewidth = 0.3) +
  scale_fill_manual(
    values = c("Pre-treatment" = "#8DA0CB", "Post-treatment" = "#E78AC3")
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2, labeller = label_parsed) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    x = "Treatment Impact",
    y = "Value",
    fill = "Stage"
  )

dev.new();print(p)
outFile <- paste0(inDir,"/figures/metrics_barchart_simulated_260106.png")
ggsave(outFile, width = 6, height = 6, units = "in", dpi=600,bg = "white")

###############################################################################
## CCDC fusion Map Statistics - on-orbit held-out test data
###############################################################################

df <- read.csv(paste0(inDir,"/data/calfire_gedi_ccdc_eval_combined.csv"))

# Function to calculate standard error
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# Summarize R2adj and RMSE by metric
summary_df <- df %>%
  group_by(metric) %>%
  summarise(
    R2adj_mean = mean(R2adj, na.rm = TRUE),
    R2adj_min  = min(R2adj, na.rm = TRUE),
    R2adj_max  = max(R2adj, na.rm = TRUE),
    R2adj_se   = se(R2adj),
    RMSE_mean  = mean(RMSE, na.rm = TRUE),
    RMSE_min   = min(RMSE, na.rm = TRUE),
    RMSE_max   = max(RMSE, na.rm = TRUE),
    RMSE_se    = se(RMSE),
    .groups = "drop"
  )

# View result
print(summary_df)

outFile <- paste0(inDir,"/data/calfire_gedi_ccdc_eval_combined_summary.csv")
write.csv(summary_df, outFile, row.names = F)

###############################################################################
## Pulse and point density summary for ALS data used in simulation
###############################################################################

version <- "260106"

dc <- read_csv(paste0(inDir,"/data/control_simulated_gedi_data_filtered_",version,".csv"))
dt <- read_csv(paste0(inDir,"/data/treatment_simulated_gedi_data_filtered_",version,".csv"))

dt1 <- dt %>% select(pointDense,beamDense, YearALS, Type, als_start, als_end)
dc1 <- dc %>% select(pointDense,beamDense, YearALS, Type, als_start, als_end)
dcombo <- rbind(dt1,dc1)

# Point and pulse density by year
dcombo %>%
  group_by(YearALS) %>%
  summarise(
    n_total = n(),
    n_treatment = sum(Type == "treatment", na.rm = TRUE),
    n_control   = sum(Type == "control",   na.rm = TRUE),
    mean_pointDense = mean(pointDense, na.rm = TRUE),
    sd_pointDense   = sd(pointDense,   na.rm = TRUE),
    mean_beamDense  = mean(beamDense,  na.rm = TRUE),
    sd_beamDense    = sd(beamDense,    na.rm = TRUE),
    .groups = "drop"
  )


# Point and pulse density -- all years
dcombo %>%
  summarise(
    n_total = n(),
    n_treatment = sum(Type == "treatment", na.rm = TRUE),
    n_control   = sum(Type == "control",   na.rm = TRUE),
    mean_pointDense = mean(pointDense, na.rm = TRUE),
    sd_pointDense   = sd(pointDense,   na.rm = TRUE),
    mean_beamDense  = mean(beamDense,  na.rm = TRUE),
    sd_beamDense    = sd(beamDense,    na.rm = TRUE),
    .groups = "drop"
  )


# Show percent by season of acquisition

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

dcombo <- dcombo %>%
  mutate(
    season = mapply(span_season_label, als_start, als_end)
  )

dcombo <- dcombo %>%
  mutate(
    als_start_year = year(als_start),
    als_end_year = year(als_end)
  )

dcombo %>%
  count(season, YearALS, als_start_year, als_end_year) %>%
  mutate(
    percent = 100 * n / sum(n)
  )

###############################################################################
## Ground-finding effects on ALS-simulated metrics
###############################################################################

version <- "260106"

response_levels_sim <- c("FHD","VDR","cover","tAGBD","tPAI","mPAI_b10","RH_25","RH_50","RH_75","RH_98", "elev")
response_labels_sim <- c("FHD","aVDR","Cover","tAGBD","tPAI",expression(mPAI["0-10m"]),"RH25","RH50","RH75","RH98","elev")

# read filtered footprint data (note: has max slope constraints)
dc <- read_csv(paste0(inDir,"/data/control_simulated_gedi_data_filtered_",version,".csv"))
dt <- read_csv(paste0(inDir,"/data/treatment_simulated_gedi_data_filtered_",version,".csv"))

dt1 <- dt %>% select(contains("gauss_"),contains("strct_"),unique_ID,YearALS) %>% select(-rhGauss_70,-rhGauss_95)
dc1 <- dc %>% select(contains("gauss_"),contains("strct_"),unique_ID,YearALS) %>% select(-rhGauss_70,-rhGauss_95)
dcombo <- rbind(dt1,dc1) # treatments & control
dcombo <- rbind(dt1) # treatments only

dcombo$shot_id <- sub(".*_(X[0-9]+_Y[0-9]+)$", "\\1", dcombo$unique_ID)
dcombo <- dcombo %>% select(-unique_ID)

# filters to footprints that have both pre and post
df1 <- dcombo %>%
  filter(!is.na(YearALS)) %>%  # Remove rows where YearALS is completely NA 
  group_by(shot_id) %>%
  filter(n_distinct(YearALS) > 1) %>%  # Only keep groups with >1 unique non-NA year
  mutate(
    min_year = min(YearALS),
    max_year = max(YearALS),
    stage = case_when(
      YearALS == min_year ~ "Year1",
      YearALS == max_year ~ "Year2",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-min_year, -max_year) %>%  # Clean up helper columns
  ungroup()

df_filter <- df1 %>%
  group_by(shot_id, stage) %>%
  slice(1) %>%   # keeps the first row per group
  ungroup()

metrics <- c(paste0("strct_",response_levels_sim),paste0("gauss_",response_levels_sim))

df_wide <- df_filter %>%
  pivot_wider(
    names_from = stage,
    values_from = all_of(c(metrics, "YearALS")),
    names_sep = "_"
  )

df_delta <- df_wide %>%
  pivot_longer(
    cols = matches("^(strct|gauss)_.*_Year[12]$"),
    names_to = c("ground_method", "metric", "year"),
    names_pattern = "(strct|gauss)_(.*)_Year([12])",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = year,
    values_from = value,
    names_prefix = "Year"
  ) %>%
  mutate(
    delta_metric = Year2 - Year1
  )

results <- df_delta %>%
  filter(!is.na(delta_metric)) %>%
  group_by(metric) %>%
  summarise(
    fit = list(feols(delta_metric ~ ground_method | shot_id, data = cur_data())),
    .groups = "drop"
  ) %>%
  mutate(tid = map(fit, broom::tidy)) %>%
  select(metric, tid) %>%
  unnest(tid) %>%
  filter(term == "ground_methodstrct")


df_delta %>%
  filter(metric %in% c("RH_25","RH_50","RH_75","RH_98")) %>%
  group_by(metric, ground_method) %>%
  summarise(sd_delta = sd(delta_metric, na.rm = TRUE), .groups="drop")

effect_sizes <- df_delta %>%
  select(metric, shot_id, ground_method, delta_metric) %>%
  pivot_wider(names_from = ground_method, values_from = delta_metric) %>%
  mutate(D = strct - gauss) %>%
  group_by(metric) %>%
  summarise(
    mean_D = mean(D, na.rm=TRUE),
    sd_D   = sd(D, na.rm=TRUE),
    d      = mean_D / sd_D,
    n      = sum(!is.na(D)),
    
    .groups="drop"
  )

effect_sizes <- df_delta %>%
  filter(!is.na(delta_metric)) %>%
  group_by(metric) %>%
  mutate(
    mean_abs_delta = mean(abs(delta_metric), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(metric, shot_id, ground_method, delta_metric, mean_abs_delta) %>%
  pivot_wider(names_from = ground_method, values_from = delta_metric) %>%
  mutate(
    D = strct - gauss
  ) %>%
  group_by(metric) %>%
  summarise(
    mean_D = mean(D, na.rm = TRUE),
    sd_D   = sd(D, na.rm = TRUE),
    d      = mean_D / sd_D,
    n      = sum(!is.na(D)),
    rel_effect = abs(mean_D) / first(mean_abs_delta),
    .groups = "drop"
  )

effect_sizes %>%
  arrange(desc(rel_effect))