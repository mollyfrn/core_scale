#Molly F. Jenkins 
#updated 05/03/2018
####Pre-analysis prep and partitioning of BBS50stop data####

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

####Bringing in BBS50 stop data and prepping it for scale partitioning####
#check name of bbs 50 stop data and correct query if necessary 
rdataretriever::datasets()
bbs = rdataretriever::fetch('breed-bird-survey-50stop') 
bbs50 = bbs$breed_bird_survey_50stop_counts  #raw counts
bbs50 = bbs50 %>%
  mutate(stateroute = as.integer(statenum*1000 + route)) #create unique state + route combo ID
write.csv(bbs50, "intermed/bbs50.csv", row.names = FALSE)  


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

