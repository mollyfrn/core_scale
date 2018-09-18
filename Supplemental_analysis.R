#Molly F. Jenkins 
#updated 09/18/2018
####Supplemental analysis detailing species temporal occupancy distributions####

#Set working directory to core_scale folder thru github i.e. setwd("C:/git/core_scale")
#years are 2000-2014 because provides 1005 routes

## Please download and install the following packages:
library(raster)
library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(dplyr)
library(fields)
library(rdataretriever)
library(stats)
library(viridis)
library(gridExtra)

####Bringing in BBS50 stop data and prepping it for scale partitioning####
#check name of bbs 50 stop data and correct query if necessary 
rdataretriever::datasets()
bbs = rdataretriever::fetch('breed-bird-survey-50stop') 
bbs50 = bbs$breed_bird_survey_50stop_counts  #raw counts
bbs50 = bbs50 %>%
  mutate(stateroute = as.integer(statenum*1000 + route)) #create unique state + route combo ID
write.csv(bbs50, "intermed/bbs50.csv", row.names = FALSE)  

#smaller data file for testing (REMOVE WHEN FINISHED EDITS)
#bbs50 = read.csv("fifty1.csv", header = TRUE)
#routes = read.csv("routes.csv", header = TRUE)

# Get subset of BBS routes (just routes) btw 2000-2014 surveyed in EVERY year
routes = bbs$breed_bird_survey_50stop_routes
routes = routes %>%
  mutate(stateroute = statenum*1000 + route)
#file of ALL BBS routes for N. Am. with Lat Longs

require(dplyr)
good_rtes = bbs50 %>% 
  filter(year >= 2000, year < 2015) %>% 
  select(year, stateroute) %>%
  unique() %>%    
  group_by(stateroute) %>%  
  count(stateroute) %>% 
  filter(n == 15) 
write.csv(good_rtes, "intermed/good_rtes.csv", row.names = FALSE)

# Subset the full BBS dataset to the routes with latlong data present 
require(dplyr)
fifty_allyears = bbs50 %>% 
  filter(year >= 2000, year < 2015) %>% 
  filter(stateroute %in% good_rtes$stateroute)
write.csv(fifty_allyears,"intermed/fifty_allyears.csv", row.names = FALSE)


# merge lat longs from routes file to the list of "good" routes (2001-2015 present all years)
require(dplyr)
good_rtes2 = routes %>% 
  filter(stateroute %in% good_rtes$stateroute) %>% 
  dplyr::select(stateroute, latitude, longitude)
write.csv(good_rtes2, "intermed/good_rtes2.csv", row.names = FALSE)
#updated 05/03/2018

#################################################################
####Temporal occupancy for each spp below, at, and above scale of a bbs route


####Data prep for calculating occupancy at and below the scale of a BBS route####
occ_counts2 = function(countData, countColumns, scale) {
  bbssub = countData[, c("stateroute", "year", "aou", countColumns)] 
  bbssub$groupCount = rowSums(bbssub[, countColumns]) 
  bbsu = unique(bbssub[bbssub[, "groupCount"]!= 0, 
                       c("stateroute", "year", "aou", "groupCount")])
  return(bbsu)
}

#read in output of pre_analysis_tidyload.R
fifty_allyears = read.csv("intermed/fifty_allyears.csv", header = TRUE) 
fifty_bestAous = fifty_allyears %>% 
  filter(aou > 2880 & !(aou >= 3650 & aou <= 3810) &
           !(aou >= 3900 & aou <= 3910) & 
           !(aou >= 4160 & aou <= 4210) & aou != 7010) #leaving out less reliable data


#loop to establish segment breakpoints for occ
c_scales = c(5, 10, 25, 50) 
output = c()

for (scale in c_scales) {
  numGroups = floor(50/scale)
  for (g in 1:numGroups) {
    groupedCols = paste("stop", ((g-1)*scale + 1):(g*scale), sep = "")
    temp = occ_counts2(fifty_bestAous, groupedCols, scale)
    temp$scale = scale
    temp$seg = g #added segment specifier 
    output = rbind(output, temp) 
  }
}

bbs_below_guide = data.frame(output)
write.csv(bbs_below_guide, "intermed/bbs_below_guide.csv", row.names = FALSE)

####Data prep for calculating occupancy above the scale of a BBS route####
#Route grouping - scales from 2 aggregated BBS routes and up
#grouping, ranking stateroutes by distance from next nearest neighbor 

occ_counts2 = function(countData, countColumns, scale) {
  bbssub = countData[, c("stateroute", "year", "aou", countColumns)] 
  bbssub$groupCount = rowSums(bbssub[, countColumns]) 
  bbsu = unique(bbssub[bbssub[, "groupCount"]!= 0, 
                       c("stateroute", "year", "aou", "groupCount")])
  return(bbsu)
}


fifty_allyears = read.csv("intermed/fifty_allyears.csv", header = TRUE) 
fifty_bestAous = fifty_allyears %>% 
  filter(aou > 2880 & !(aou >= 3650 & aou <= 3810) & 
           !(aou >= 3900 & aou <= 3910) & 
           !(aou >= 4160 & aou <= 4210) & aou != 7010) 
#leaving out owls, waterbirds as less reliable data

good_rtes2 = read.csv("intermed/good_rtes2.csv", header = TRUE) 
require(fields)
distances = rdist.earth(matrix(c(good_rtes2$longitude, 
                                 good_rtes2$latitude), ncol=2),
                        matrix(c(good_rtes2$longitude, 
                                 good_rtes2$latitude), ncol=2),
                        miles=FALSE, R=6371)
dist.df = data.frame(rte1 = rep(good_rtes2$stateroute, each = nrow(good_rtes2)),
                     rte2 = rep(good_rtes2$stateroute, times = nrow(good_rtes2)),
                     dist = as.vector(distances))
write.csv(dist.df, "intermed/dist_df.csv", row.names = FALSE) 


c_scales = c(50)
output = c()

for (scale in c_scales) {
  numGroups = floor(50/scale)
  for (g in 1:numGroups) {
    groupedCols = paste("stop", ((g-1)*scale + 1):(g*scale), sep = "")
    temp = occ_counts2(fifty_bestAous, groupedCols, scale)
    output = rbind(output, temp) 
  }
}

bbs_above_guide = data.frame(output)
write.csv(bbs_above_guide, "intermed/bbs_above_guide.csv", row.names = FALSE)


####Paring bbs_below down to 983 routes####
dist.df = read.csv("intermed/dist_df.csv", header = TRUE)
bbs_below_guide = read.csv("intermed/bbs_below_guide.csv", header = TRUE)


#filter out to only routes that are up to 1000km radius away from each other 
far = dist.df %>% 
  arrange(rte1, dist) %>% 
  group_by(rte1) %>% 
  slice(66)

hist(far$dist)
far2 = far %>% 
  filter(dist < 1000)

bbs_below_guide = bbs_below_guide %>% 
  filter(stateroute %in% far2$rte1)

#test should return [1] 1 2; 25 stop scale should only have two segments per route
test = bbs_below_guide %>% 
  filter(scale == "25")
unique(test$seg)

###################################################################

