# Analysis of the proportion of core species in communities across different scales
# author: Molly F. Jenkins
# updated: 05/03/2017

# Summary: The following script demonstrates the first original cutoffs of >2/3 and <1/3, 
# which can be viewed in the core_scale_analyses.R script. 
#This script explores cutoffs of: 2/3 $ 1/3 (original occ-scale-processing and creation of bbs_allscales.csv data)
#The supplemental cutoff scripts can be viewed in allscales_coefs50.R and allscales_coefs80.R
###################################################################################################################

# setwd("C:/git/core_scale")
# Please download and install the following packages:
library(raster)
library(tidyverse)
library(fields)
library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(nlme)
library(gridExtra)
library(wesanderson)
library(stats)
library(viridis)


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


#group by scale and segment and THEN take means of segments
require(tidyverse)
uniqrtes = unique(bbs_below_guide$stateroute) 
scale = unique(bbs_below_guide$scale)
rte_segments = unique(bbs_below_guide$seg)

output = data.frame(focalrte = NULL,
                    scale = NULL, 
                    meanOcc = NULL, 
                    pctCore = NULL,
                    pctTran = NULL,
                    aveN = NULL)


for (r in uniqrtes) { #for each focal route
  for (nu in scale) { #for each level of scale aggregated to each focal route
      focal_c = bbs_below_guide %>% 
      filter(stateroute == r & scale == nu)
    
    for (s in unique(focal_c$seg)) {
        focal_clustr = focal_c %>% 
        filter(seg == s) 
      
      abun.summ = focal_clustr %>% 
        group_by(year) %>%  
        summarize(totalN = sum(groupCount)) %>%
        summarize(aveN = mean(totalN), 
                  stateroute = r)
      
      occ.summ = focal_clustr %>% 
        dplyr::select(year, aou) %>% 
        distinct() %>% 
        count(aou) %>% #how many times does that AOU show up in that clustr that year 
        mutate(occ = n/15, scale = nu) %>% 
        summarize(focalrte = r, 
                  scale = nu, 
                  meanOcc = mean(occ), 
                  pctCore = sum(occ > 2/3)/length(occ),
                  pctTran = sum(occ <= 1/3)/length(occ)) 
      
      
      occ2 = occ.summ %>% 
        group_by(focalrte) %>% 
        summarize(scale = nu, 
                  meanOcc = mean(meanOcc), #mean of means, community mean  
                  pctCore = mean(pctCore), 
                  pctTran = mean(pctTran)) %>%
        left_join(abun.summ, by = c('focalrte' = 'stateroute'))
      
      output = rbind(output, occ2)
      print(paste("Focal rte", r, "#' rtes sampled", nu)) #for viewing progress
      
      
    } #segment loop 
  } #n loop
  
} #r loop

bbs_below = as.data.frame(output)
write.csv(bbs_below, "intermed/bbs_below_new.csv", row.names = FALSE)

#average for unique scale-route combo (currently duplicates based on segments)
bbs_below_avgs = bbs_below %>% 
  group_by(focalrte, scale) %>% 
  summarize(meanOcc = mean(meanOcc), 
            pctCore = mean(pctCore), 
            pctTran = mean(pctTran), 
            aveN = mean(aveN))

write.csv(bbs_below_avgs, "intermed/bbs_below_avgs.csv", row.names = FALSE)
#updated 05/03/2018


####Above-scale and merging code####
dist.df = read.csv("intermed/dist_df.csv", header = TRUE)
bbs_above_guide = read.csv("intermed/bbs_above_guide.csv", header = TRUE)


#filter out to only routes that are up to 1000km radius away from each other
far = dist.df %>% 
  arrange(rte1, dist) %>% 
  group_by(rte1) %>% 
  slice(66)
hist(far$dist)
far2 = far %>% 
  filter(dist < 1000)

bbs_above_guide = bbs_above_guide %>% 
  filter(stateroute %in% far2$rte1)


uniqrtes = unique(bbs_above_guide$stateroute) 
numrtes = 2:66 
output = data.frame(focalrte = NULL,
                    scale = NULL, 
                    meanOcc = NULL, 
                    pctCore = NULL,
                    pctTran = NULL,
                    maxdist = NULL,
                    aveN = NULL)

for (r in uniqrtes) { #for each focal route
  for (nu in numrtes) { #for each level of scale aggregated to each focal route
    
    tmp_rte_group = dist.df %>% #changes with size of nu but caps at 66
      filter(rte1 == r) %>% 
      top_n(66, desc(dist)) %>% #fixed ordering by including arrange parameter 
      arrange(dist) %>%
      slice(1:nu) %>% 
      dplyr::select(everything()) %>% data.frame()
    
    #this chunk takes varying list from above and uses it to subset the bbs data 
    #so that occ can be calculated for AOUs anywhere within the cluster 
    focal_clustr = bbs_above_guide %>% 
      filter(stateroute %in% tmp_rte_group$rte2) 
    
    abun.summ = focal_clustr %>% #abundance
      group_by(year) %>%  
      summarize(totalN = sum(groupCount)) %>%
      summarize(aveN = mean(totalN), 
                stateroute = r)
    
    occ.summ = focal_clustr %>% 
      dplyr::select(year, aou) %>% 
      distinct() %>% 
      count(aou) %>% #how many times does that AOU show up in that clustr that year 
      mutate(occ = n/15, scale = nu) %>% 
      summarize(focalrte = r, 
                scale = nu, 
                meanOcc = mean(occ), 
                pctCore = sum(occ > 2/3)/length(occ),
                pctTran = sum(occ <= 1/3)/length(occ), 
                maxdist = max(tmp_rte_group$dist)) 
    
    occ2 = occ.summ %>% 
      group_by(focalrte) %>% 
      summarize(scale = nu, 
                meanOcc = mean(meanOcc), #mean of means, community mean  
                pctCore = mean(pctCore), 
                pctTran = mean(pctTran), 
                maxdist = mean(maxdist)) %>%
      left_join(abun.summ, by = c('focalrte' = 'stateroute'))
    
    output = rbind(output, occ2)
    print(paste("Focal rte", r, "#' rtes sampled", nu)) #for viewing progress

  } #n loop
  
} #r loop

bbs_above = as.data.frame(output)
write.csv(bbs_above, "intermed/bbs_above.csv", row.names = FALSE)
#updated 05/03 

####Merging occ data across scales####
bbs_above = read.csv("intermed/bbs_above.csv", header = TRUE)
bbs_below = read.csv("intermed/bbs_below_avgs.csv", header = TRUE)

#use area in km2 as appropriate scale metric
#area in km by # of routes * 50 stops in each rte * area of a stop 
bbs_below2 = bbs_below %>% 
  mutate(maxdist = c("NA")) %>%
  dplyr::select(focalrte, scale, everything()) %>%
  mutate(area = bbs_below$scale*(pi*(0.4^2)), 
         scale = paste("seg", scale, sep = ""))

bbs_above = bbs_above %>% 
  dplyr::mutate(area = scale*50*(pi*(0.4^2))) %>% 
  dplyr::select(focalrte, scale, meanOcc, pctCore, 
                pctTran, aveN, maxdist, area) 

bbs_above$scale = as.factor(bbs_above$scale)
bbs_allscales = rbind(bbs_below2, bbs_above)

bbs_allscales = bbs_allscales %>%
  mutate(logA = log10(area), 
         logN = log10(aveN),
         lnA = log(area), #log is the natural log 
         lnN = log(aveN))

write.csv(bbs_allscales, "intermed/bbs_allscales.csv", row.names = FALSE) 
#saved 05/03

####Extract coefficients from scale-occupancy relationships for analysis####
#normal cutoff of 67%

bbs_allscales = read.csv("intermed/bbs_allscales.csv", header = TRUE)
levels(bbs_allscales$scale)
unique(bbs_allscales$scale)
length(unique(bbs_allscales$scale))
length(unique(bbs_allscales$focalrte)) 
bbs_allscales2 = bbs_allscales %>% 
  dplyr::select(-scale, -maxdist) #irrelevant post-merge

PCA.df = data.frame(stateroute = numeric(), PCA.min = numeric(), PCA.max = numeric(), 
                    PCA.slope = numeric(), 
                    PCA.mid = numeric(), 
                    PCA.curvature = numeric())
PCN.df = data.frame(stateroute = numeric(), PCN.min = numeric(), PCN.max = numeric(), 
                    PCN.slope = numeric(), 
                    PCN.mid = numeric(), 
                    PCN.curvature = numeric())


stateroutes = unique(bbs_allscales2$focalrte)

for(s in stateroutes){
  logsub = subset(bbs_allscales2, bbs_allscales2$focalrte == s)  
  #PCA 
  PCA.min = logsub$pctCore[logsub$logA == min(logsub$logA)]
  PCA.max = logsub$pctCore[logsub$logA == max(logsub$logA)]
  PCA.mid = min(logsub$logA[logsub$pctCore >= 0.5]) 
  PCA.slope = ((PCA.max - PCA.min)/(max(logsub$logA) - min(logsub$logA)))
  #want the FIRST instance where it hits this range -> how? minimum scale at which it does that
  
    PCA.obline = logsub$pctCore #observed, the vector of x1
    b = PCA.min -(PCA.slope*min(logsub$logA)) # b = y1 - m*x1 solve for b
    PCA.pline = PCA.slope*logsub$logA+b #linear vector of occs that lie btw/min & max 
  PCA.curvature = sum(PCA.obline-PCA.pline) 
  #AUC proxy - taking diff between actual and predicted vals at EVERY scale
  
  PCAmodel = data.frame(stateroute = s, PCA.min, PCA.max, PCA.slope, 
                        PCA.mid, PCA.curvature)
  
  PCA.df = rbind(PCA.df, PCAmodel)
  
  
  #PCN 
  PCN.min = logsub$pctCore[logsub$logN == min(logsub$logN)]
  PCN.max = logsub$pctCore[logsub$logN == max(logsub$logN)]
  PCN.mid = min(logsub$logN[logsub$pctCore >= 0.5]) 
  PCN.slope = ((PCN.max - PCN.min)/(max(logsub$logN) - min(logsub$logN)))
  
    PCN.obline = logsub$pctCore 
    b2 = PCN.min -(PCN.slope*min(logsub$logN)) 
    PCN.pline = PCN.slope*logsub$logN+b2 
  PCN.curvature = sum(PCN.obline-PCN.pline) 
  
  PCNmodel = data.frame(stateroute = s, PCN.min, PCN.max, 
                        PCN.slope, PCN.mid, PCN.curvature)
  
  PCN.df = rbind(PCN.df, PCNmodel) 
  
}  

#join all together using inner_join by focal rte 
core_coefs = PCA.df %>% 
  inner_join(PCN.df, 
             PCA.df, 
             by = "stateroute") %>% 
  distinct()

write.csv(core_coefs, "intermed/core_coefs.csv", row.names = FALSE) 
#updated 05/08

