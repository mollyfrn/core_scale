# core scale: a subset of the core-transient project 
Master's thesis research: data, methods, results, and figures


This repo contains the final workflow and scripts of the core-transient scale analysis project. The workflow and order of scripts is as follows: 

1) pre_analysis_tidyload.R - This script downloads Breeding Bird Survey (BBS) 50-stop resolution data through ecoretriever. It also pares the data down to only routes with 15 consecutive years (2000-2014) of data available. 

2) 
..a) core_scale_analyses.R - This script takes cleaned BBS 50-stop data from the previous script and renders the proportion of core species for ecological assemblages across 69 spatial scales. It also generates descriptive parameters characterizing the relationship between the proportion of core species and scale for all sampled assemblages. 
..b) core_scale_analyses50.R - This supplemental script takes cleaned BBS 50-stop data from the previous script and renders the proportion of core species for ecological assemblages across 69 spatial scales, using a cutoff threshold of 50% (e.g. 2/4) for determing the proportion of core species in a given assemblage.It also generates descriptive parameters characterizing the relationship between the proportion of core species and scale for all sampled assemblages. 
..c) core_scale_analyses80.R - This supplemental script takes cleaned BBS 50-stop data from the previous script and renders the proportion of core species for ecological assemblages across 69 spatial scales, using a cutoff threshold of 80% (e.g. 4/5) for determing the proportion of core species in a given assemblage.It also generates descriptive parameters characterizing the relationship between the proportion of core species and scale for all sampled assemblages. 

3) envs_processing.R - This script downloads and cleans elevation and NDVI raster data for use as proxies of environmental heterogeneity, based on site locations and the samples scales. 

4) coef_hetero_analyses.R - This script tests the correlation of the core-scale relationship parameters and environmental heterogeneity. 

5) figures_script.R - This script uses intermediate data files and analyses generated by the previous scripts to create descriptive figures. 


