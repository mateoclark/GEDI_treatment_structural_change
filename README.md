# GEDI_treatment_structural_change

Code and data for paper "Changes in GEDI-based measures of forest structure after large California wildfires relative to pre-fire conditions"

Clark, M. L., Bailey, T., Burns, P., Gutman, T., Haines, M., Wotring, I., Hakkenberg, C., & Goetz, S. (in review). Detecting structural change from forest fuel treatments with GEDI spaceborne lidar across California, USA. 

Matthew L. Clark, Ph.D.
Department of Geography, Environment, and Planning
Sonoma State University, California USA
matthew.clark@sonoma.edu

**Code**

*0_calcluate_ccdc_change.R*
- Calculates CCDC GEDI interpolation change in treatments

*0_calcluate_l89ndmi_change.R*
- Calculates Landsat 8/9 NDMI change in treatments

*0_calcluate_l89ndvi_change.R*
- Calculates Landsat 8/9 NDVI change in treatments

*0_calcluate_palsar2nbi_change.R*
- Calculates PALSAR2 normalized backscatter index change in treatments

*0_calcluate_s1nbi_change.R*
- Calculates Sentinel-1 normalized backscatter index change in treatments

*0_calcluate_s2ndmi_change.R*
- Calculates Sentinel-2 NDMI change in treatments

*0_calcluate_s2ndvi_change.R*
- Calculates Sentinel-2 NDVI change in treatments

*0_organize_control.R*
Organizes control polygons

*0_organize_treatments.R*
Organizes treatment polygons

*1_prepare_master_data_file.R*  
*1_prepare_master_data_file_simulated.R*
- Loads in GEDI structure data and treatment and control polygons, overlays the shots with polygon attributes,
and outputs treatment and control .csv files with information needed for analysis. On-orbit and simulated versions.

*2_footprint_gedi_metric_difference.R*  
*2_footprint_gedi_metric_difference_simulated.R*
- Calculates average change metric differences for treatments and nearest control area for similar dates. On-orbit and simulated versions.

*3_gedi_repeated_measures_ANOVA.R*  
*3_gedi_repeated_measures_ANOVA_simulated.R*
- Repeat-measures ANOVA, or mixed-effects model. On-orbit and simulated versions.

*4_spatial_regression_univariate.R*  
*4_spatial_regression_univariate_simulated.R*
- Univariate spatial regression models. On-orbit and simulated versions.

*5_footprint_simulated_gedi_metric_ccdc_comparison.R*
- Compares simulated GEDI to CCDC interpolated metrics

*6_graphs.R*
- Graphs for manuscript

**Data**

To do
