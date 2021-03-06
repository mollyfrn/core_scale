#Percent-Core Occupancy-scale analysis: alternative to mean occupancy for ease of reader + metric 

####Occ-scale analysis####

##feat. alternative to curve-fitting parameters: slope, intercept, and x value @ scale of 3 aggregated routes.
# author: Molly F. Jenkins
# date: 06/27/2017


# setwd("C:/git/core_scale")
#'#' Please download and install the following packages:
library(raster)
library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(fields)
library(tidyverse)
library(nlme)
library(gridExtra)
library(wesanderson)
library(stats)
library(viridis)

# To run this script, you need temperature, precip, etc data, 
# which are currently stored in the following directories off of github: 

# Data directories
tempdatadir = '//bioark.ad.unc.edu/HurlbertLab/GIS/ClimateData/BIOCLIM_meanTemp/'
precipdata = '//bioark.ad.unc.edu/HurlbertLab/GIS/ClimateData/2-25-2011/prec/'
ndvidata = "//bioark.ad.unc.edu/HurlbertLab/GIS/MODIS NDVI/"
BBS = '//bioark.ad.unc.edu/HurlbertLab/Jenkins/Intermediate scripts/BBS scaled/'


####Occ-line equation for pct Core scale relationship####

#Occupancy-scale analysis

# author: Molly F. Jenkins
# date: 06/27/2017

####Summary####
#This script takes a dataset of BBS routes with occupancy values calculated for every scale, 
#scales ranging in size from 1/10th of a BBS route to 66 aggregated BBS routes. 
#Using this data we analyze the relationship between temporal occupancy of communities and scale. 
#We characterize this relationship in a series of simple linear models between occupancy and scale 
#(using log(area) as our variable of spatial scale). From these models and our existing data, 
#for every focal route we traced how occupancy changes across scale, from smallest to largest scale. 
#We characterized these changes through a series of variables: 
#1)the occupancy value at the minimum scale for a focal route 
#2)the occupancy value at the maximum scale for a focal route 
#3)the slope of the line linking the minimum and maximum values 
#the steepness (or flatness) of this line corresponds with the rate of accumulation of core species 
#for a given community 
#4)the scale at which mean occupancy first reaches or surpasses 50% 
#At what scale, at what area in km must we reach in order to reliably see consistent occupancy 
#PREDICTION: in areas of fairly uniform habitat type, this scale should be lower 
#PREDICTION: in areas of fairly high habitat heterogeneity, this scale should be higher 
#5)and we look at the straightness or curvature of the actual data 
#as compared to the data derived from our model.
#PREDICTION: focal routes with larger "curvy" values will occur in regions of greater habitat heterogeneity, 
#and deviance from the line will correspond with the greater environmental variance
#associated with the location of that focal route.
#i.e. curvy is a proxy for AUC, where greater deviance of the predicted vals vs the actual vals will result 
#in larger overall curvy values and larger area under the curve

#We explore the variation in this relationship 
#and attempt to characterize whether or not it is best explained by habitat heterogeneity 
#using the variation present across several environmental variables as proxies for habitat heterogeneity. 

#### CURRENT ISSUES #### #revised 09/10/2017

#1) Have below-scale duplicates w/diff starting locations i.e. 
#"scale" = Factor w/83 levels "10-1", "10-2", "10-3", "10-4", "10-5", "5-1", "5-2", "5-3", "5-4", "5-5", "5-6", "5-7", "5-8", etc. 
#no duplicates for 1:66 obvi, these were calculated correctly 
#can I have multi values for a single scale or should I calc avg occ and pctCore and pctTran across the segments w/in a route? 
#i.e. occ/15 at a scale, but scales designated by segments and occ calcd initially based on the starting stop # 
# was told this summer that that was alright and that we were pointedly NOT aggregating across the lower scales, and that's fine 
#BUT: it DOES mean we will have pseudo-duplicates at the lower scales w/diff starting points of segments, and so 
#multiple plotting points for the lower scales instead of a single representative point. This necessitates the "min" qualifier
#in calculating some of the coefficients.  

#2) In calculating the AUC proxy "curvy" values, should I be squaring the differences before summing them? 
#Because we are interested in the magnitude of deviance from the observed values in instances of greater heterogeneity, 
#and not the specific direction - simply adding them could cancel out or diminish the magnitude or gulf of deviance 
#in areas where the deviance is large in a negative direction at some scales, and large in a positive direction in others. 
#These stateroutes would express similar "curviness" values as areas of high homogeneity, where this is little net difference 
#and a smaller magnitude of deviance from the observed values. Because of this, I think we should be squaring and summing. 

#3) Curviness is a great proxy for calculating the AUC per focal route, but it doesn't allow us to see how that changes across scales. 
#Considering having a second df output from the creation of coefs that allows us to see how those individual deviances change 
#esp across scales. 


####Env analysis####
#Background data prep is just envs_processing script
core_coefs = read.csv("intermed/core_coefs.csv", header = TRUE) 
env_all = read.csv("intermed/env_all.csv", header = TRUE) 

core_env_coefs = env_all %>%
  inner_join(core_coefs, by = "stateroute") %>%
  dplyr::select(stateroute, scale, area, ndvi.var, elev.var, PCA.min, PCA.max, PCA.slope, PCA.mid, PCA.curvature,
         PCN.min, PCN.max, PCN.slope, PCN.mid, PCN.curvature) %>% 
  filter(scale > 2) #ensure env analyses only happening with 2+ agg routes
write.csv(core_env_coefs, "intermed/core_env_coefs.csv", row.names = FALSE)
#updated 05/09


####Coef & habitat heterogeneity models####
core_env_coefs = read.csv("intermed/core_env_coefs.csv", header = TRUE) #same # of cols as old 

#check out cov matrix to inform model generation and predictions:
covmatrix = round(cor(core_env_coefs[, 4:ncol(core_env_coefs)]), 2) #round to 2 decimal digits
covmatrix = as.data.frame(covmatrix)
write.csv(covmatrix, "intermed/core_covmatrix.csv", row.names = FALSE)
#mean and var - interpret how covary with coefs and direction - i.e. ndvi mean covaries positively with pmin, 
#but elev mean covaries negatively with pmin 
#and both variances covary negatively with pmin 
#ndvi mean covaries negatively with pthresh, elev mean positively
#no other variables have this divergent relationship in directionality of covariance, 
#ndvi var and elev var always in unison when strong
#hab het variance measures: pslope and pthresh positive, pmin and pmax both negative covariance


####UP NEXT: run series of models####
# nested loop for examining variation in coefs/fitted curves explained by env heterogeneity 

core_rsqrd_hetero = data.frame(dep = character(), ind = character(), 
                          r2 = numeric(), adjr = numeric(), corr_r = numeric(), uppr = numeric(), lowr = numeric())

for (d in 4:5) { #adjust columns appropriately -> make sure correct order of ind and dep vars!
  for (i in 6:15) {
    tempmod = lm(core_env_coefs[,d] ~ core_env_coefs[,i])
    tempcor = cor.test(core_env_coefs[,d], core_env_coefs[,i], method = "pearson")
    
    
    tempdf = data.frame(dep = names(core_env_coefs)[d], 
                        ind = names(core_env_coefs)[i], 
                        r2 = summary(tempmod)$r.squared, 
                        adjr = summary(tempmod)$adj.r.squared, 
                        corr_r = as.numeric(tempcor$estimate), 
                        uppr = as.numeric(tempcor$conf.int[2]), 
                        lowr = as.numeric(tempcor$conf.int[1]))

    core_rsqrd_hetero = rbind(core_rsqrd_hetero, tempdf)
  }
}


write.csv(core_rsqrd_hetero, "intermed/core_rsqrd_hetero.csv", row.names = FALSE) 
#updated 05/07 using corrected hab_het vals, only variances characterizing sites


####Visually Characterizing measures of habitat heterogeneity####
core_rsqrd_hetero = read.csv("intermed/core_rsqrd_hetero.csv", header = TRUE)
# hab_het = read.csv("scripts/R-scripts/scale_analysis/hab_het.csv", header = TRUE)

r_plot = ggplot(data = core_rsqrd_hetero, aes(y = corr_r))+geom_col(aes(x=dep))+facet_wrap(~ind)+
  theme_bw()
r_plot 

#scale on x and r on y, panel by coef of interest, line color by var measure

# goal plot -> ggplot(envcoefs, aes(x = scale, y = corr_r))+geom_line(aes(color = dep))+facet_wrap(~ind)
#I want a corr_r value for every dep and ind variable at every scale, for every focal
#for every scale, for every focal route - will have a LOT - maybe just do a subset for meeting 

#the correlation coefficients themselves won't change, bc representative of the overall 
#occ-scale relationship, that's fine - the hab_het vals will change though bc measures 
#at each scale 
#starting at scale of 1 since that's lowest res we have for habhet across scales, 


#rerun previous dep/ind loop with new mods
core_env_coefs = read.csv("intermed/core_env_coefs.csv", header = TRUE)

core_scales_hetero = data.frame(dep = character(), ind = character(), 
                           r2 = numeric(), adjr = numeric(), corr_r = numeric(), 
                           uppr = numeric(), lowr = numeric(), scale = numeric())
scales = unique(core_env_coefs$scale)


for (s in scales) {
  env_coefs2 = core_env_coefs %>% 
    filter(scale == s)
  for (d in 4:5) { #adjust columns appropriately -> make sure correct order of ind and dep vars!
    for (i in 6:15) {
      tempmod = lm(env_coefs2[,d] ~ env_coefs2[,i])
      tempcor = cor.test(env_coefs2[,d], env_coefs2[,i], method = "pearson")
      
      
      tempdf = data.frame(dep = names(env_coefs2)[d], 
                          ind = names(env_coefs2)[i], 
                          r2 = summary(tempmod)$r.squared, 
                          adjr = summary(tempmod)$adj.r.squared, 
                          corr_r = as.numeric(tempcor$estimate),
                          uppr = as.numeric(tempcor$conf.int[2]),
                          lowr = as.numeric(tempcor$conf.int[1]),
                          scale = s)
      
     core_scales_hetero = rbind(core_scales_hetero, tempdf)
    }
  }
}

write.csv(core_scales_hetero, "intermed/core_scales_hetero.csv", row.names = FALSE) 
#updated 05/07 using corrected hab_het vals, only variances characterizing sites

#INTERPRETATION: 
#elevation and ndvi at the highest scales explain more variation in the pmin and pthresh values 
#compared to any other measures of habitat heterogeneity 
#However, elevation explains more variation than ndvi at the top scales - AND 
#the gulf between ndvi and elevation in terms of explanatory power widens at the larger scales 
# - the two measures are closer together in importance at the lowest scales. 
#NDVI variance at local scale actually has an 'outlier' point further above even local elev variation 
#depending on which coefficient metric that corresponds to, it may be that ndvi explains more 
#local variation in heterogeneity and variation at the lower scales of the occ-scale relationship 
#while elevational heterogeneity explains variation at the higher scales of the occ-scale relationship?  
