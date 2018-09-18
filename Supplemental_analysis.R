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
    groupedCols = paste("Stop", ((g-1)*scale + 1):(g*scale), sep = "")
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
far = dist.df %>% group_by(rte1) %>%
  arrange(dist, .by_group = TRUE) %>% 
  slice(1:66)

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
#right, so, I want to show species temporal occupancies and how those change with scale 
#we would expect vultures/scavengers and falconiformes to be core at much higher scales than passeriformes 
#given their reqs by life history for large foraging territories 
#while we would expect passeriformes to reach core occs at lower scales onward and plateau earlier
#given their smaller resource reqs, in spite of their capabilities for covering large dists 

#Preliminary decision to group by ORDER within AOU codes as grouping variable bc easier to work with 
#orders have common life history characterisitics and traits that may distinguish spp therein 
#and lead to some orders being more or less likely to plateau at earlier/lower scales than others


####Supplemental figures mirroring figs 4-6####
####Dummy data and real data pre-plotting tidy for figures 1 & 4####
fifty_allyears = read.csv("intermed/fifty_allyears.csv", header = TRUE) #using updated version, 50 stop data, 07/12


fifty_bestAous = fifty_allyears %>% 
  filter(aou > 2880 & !(aou >= 3650 & aou <= 3810) & !(aou >= 3900 & aou <= 3910) & 
           !(aou >= 4160 & aou <= 4210) & aou != 7010) #leaving out owls, waterbirds as less reliable data

#use occ_counts function for calculating occupancy at any scale to get raw occs for distribution plots
occ_counts = function(countData, countColumns, scale) {
  bbssub = countData[, c("stateroute", "year", "aou", countColumns)] #these are our grouping vars
  bbssub$groupCount = rowSums(bbssub[, countColumns]) 
  bbsu = unique(bbssub[bbssub[, "groupCount"]!= 0, c("stateroute", "year", "aou")]) 
  
  abun.summ = bbssub %>% #abundance
    group_by(stateroute, year) %>%  
    summarize(totalN = sum(groupCount))  #we want to go further and summarize across focal + secondary rtes tho
  
  occ.summ = bbsu %>% #occupancy
    count(stateroute, aou) %>%
    mutate(occ = n/15, scale = scale) %>% #, #may want to get rid of, this is at the column-counting scale
    #scale = scale) %>%
    left_join(abun.summ, by = 'stateroute')
  return(occ.summ)
}


# Generic calculation of occupancy for a specified scale
#fix to run all at once, so no sep run for above-scale, USE occ-counts for both 
b_scales = c(5, 10, 25, 50)
output = c()
for (s in b_scales) {
  numGroups = floor(50/s)
  for (g in 1:numGroups) {
    groupedCols = paste("stop", ((g-1)*s + 1):(g*s), sep = "")
    temp = occ_counts(fifty_bestAous, groupedCols, s) 
    output = rbind(output, temp) 
  } 
}

min_dist = output
min_dist = min_dist[, -3] #removes vestigial "n" count column 

#need to avg occs between unique stateroute-aou pairs since 5 for every 1 
min_dist3 = min_dist %>% 
  group_by(aou, stateroute, scale) %>% 
  summarise(occ = mean(occ)) %>% dplyr::select(everything()) 
min_out = as.data.frame(min_dist3)

write.csv(min_out, "intermed/min_out.csv", row.names = FALSE)

#check 
min_out2 = min_out %>% 
  filter(scale == "50")

fig1a = ggplot(min_out2, aes(occ))+
  geom_density(bw = "bcv", kernel = "gaussian", n = 2000, na.rm = TRUE)+
  labs(x = "Proportion of time present at site", y = "Probability Density", title = "Single Route Scale")+ 
  theme_classic() #coord_cartesian(xlim = c(0, 1), ylim = c(0, 2.5))+
fig1a #looks good!

#######################

#distribution at the 2-66 rte group scales#
dist.df = read.csv("intermed/dist_df.csv", header = TRUE)
bbs_above_guide = read.csv("intermed/bbs_above_guide.csv", header = TRUE) #generated by core_scale_analyses
#groupcounts for each aou for each year at scale of ONE stateroute 

uniqrtes = unique(bbs_above_guide$stateroute) #all routes present are unique, still 953 which is great
scales = c(2, 4, 8, 16, 32, 66) 
max_out = c()

for (nu in scales){
  #test example route 2010 and nu at 57 routes -> large scale, should have high occ 
  for (r in uniqrtes) { #for each focal route
    tmp_rte_group = dist.df %>% #changes with size of nu but caps at 66
      filter(rte1 == r) %>% 
      top_n(66, desc(dist)) %>% #fixed ordering by including arrange parm, 
      #remove/skip top row 
      arrange(dist) %>%
      slice(1:nu) %>% 
      dplyr::select(everything()) %>% data.frame()
    
    
    focal_clustr = bbs_above_guide %>% 
      filter(stateroute %in% tmp_rte_group$rte2) 
    
    occ.summ = focal_clustr %>% 
      dplyr::select(year, aou) %>% #duplicates remnant of distinct secondary routes - finally ID'd bug
      distinct() %>% #removing duplicates 09/20
      count(aou) %>% #how many times does that aou show up in that clustr that year 
      dplyr::mutate(occ = n/15, stateroute = r, scale = nu) 
    
    max_out = rbind(max_out, occ.summ)
    
  }
}

max_out = max_out[, -2] #rm vistigal "n" column 
max_out = as.data.frame(max_out)

fig1c = ggplot(max_out, aes(occ))+
  geom_density(bw = "bcv", kernel = "gaussian", n = 2000, na.rm = TRUE)+
  labs(x = "Proportion of time present at site", y = "Probability Density", title = "Maximum Scale")+theme_classic()
#so it was the limits giving me crap in the original 
fig1c

write.csv(max_out, "intermed/max_out.csv", row.names = FALSE)
##################

#read in occ data for merge for all_fig output creation 
min_out = read.csv("intermed/min_out.csv", header = TRUE)
max_out = read.csv("intermed/max_out.csv", header = TRUE)
dist.df = read.csv("intermed/dist_df.csv", header = TRUE)
#groupcounts for each aou for each year at scale of ONE stateroute 

#filter out to only routes that are up to 1000km radius away from each other before analyses 
far = dist.df %>% arrange(rte1, dist) %>% group_by(rte1) %>% slice(66)
hist(far$dist)
far2 = far %>% filter(dist < 1000)

min_out2 = min_out %>% filter(stateroute %in% far2$rte1)
max_out2 = max_out %>% filter(stateroute %in% far2$rte1)
#updated 05/07

#organize by scales; label and differentiate scales so that below-rtes are appropriately smaller
min_out = min_out2 %>% 
  dplyr::select(stateroute, aou, occ, scale) %>% 
  dplyr::mutate(area = scale*(pi*(0.4^2))) %>% #scale corresponds to the number of stops
  dplyr::select(stateroute, aou, occ, area)

max_out = max_out2 %>% 
  dplyr::select(stateroute, aou, occ, scale) %>% 
  dplyr::mutate(area = scale*50*(pi*(0.4^2))) %>% #scale corresponds to the number of agg routes; 50 stops per rte
  dplyr::select(stateroute, aou, occ, area)


all_fig = rbind(max_out, min_out)
length(unique(all_fig$stateroute)) #983, as it should be 
write.csv(all_fig, "intermed/all_figoutput.csv", row.names = FALSE)

####Plotting how CT distributions change across scale, using AOU INSTEAD OF AREA####
all_fig = read.csv("intermed/all_figoutput.csv", header = TRUE)
#all_fig$area = as.factor(all_fig$area)
all_fig$aou_f = factor(all_fig$aou),#want to bin by orders
                        levels = c(3, 5, 13, 25, 50, 101, 201, 402, 804, 1659), 
                        labels = c("2.5, 5 point count stops", "5", "13", "25, 1 BBS route", "50", "101", "201", "402", "804", "1659, 66 aggregate BBS routes")) 
all_fig = all_fig %>% 
  #first I have to take levels for area_f and lump everything together that isn't 50/25, 1 BBS route
  mutate(area_spec = as.numeric(area_f))

all_fig$area_spec[all_fig$area_spec != 4] <- 0

all_fig$area_spec = factor(all_fig$area_spec)


#make sure legend includes both color and dash 
all_figplot = ggplot(all_fig, aes(occ, group = area_f, color = area_f, linetype = area_spec))+
  stat_density(geom = "path", position = "identity", bw = "bcv", kernel = "gaussian", n = 4000, na.rm = TRUE, size = 1.3)+
  #stat_density(singlerte, aes(occ), geom = "path", position = "identity", bw = "bcv", kernel = "gaussian", n = 4000, na.rm = TRUE, size = 1.6)+
  labs(x = "Proportion of time present at site", y = "Probability Density")+theme_classic()+
  scale_color_viridis(discrete = TRUE, name = expression("Spatial Scale in km"^{2}))+
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 30, color = "black"))+
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30))+
  theme(legend.position = c(0.50, 0.50))+guides(linetype = FALSE) + theme(legend.key = element_rect(size = 4, color = 'white'), legend.key.size = unit(2, 'lines'))
all_figplot
ggsave(file = "output/Figure4.tiff", plot = all_figplot)


####Plotting NULL all routes with 3 highlighted "types####
bbs_allscales = read.csv("intermed/bbs_allscales.csv", header = TRUE)
bbs_allscales = bbs_allscales %>% 
  dplyr::filter(logN != "NA")

core_coefs = read.csv("intermed/core_coefs.csv", header = TRUE) #AUC etc.

coefs_ranked = core_coefs %>% 
  arrange(PCA.curvature) #middle teal line should be least curvy

coefs_avgs = core_coefs %>% 
  summarise_all(funs(mean, sd)) %>% 
  mutate(focalrte = "99999") 


#compare rtes in homogeneous vs heterogeneous regions, illustrate in color to prove point 
env_all = read.csv("intermed/env_all.csv", header = TRUE) #AUC etc.
ndvi_ranked = env_all %>% 
  group_by(stateroute) %>% 
  summarize(ndvi_m = mean(ndvi.var)) %>%
  arrange(desc(ndvi_m)) 
#lowest var in NDVI: rtes 72151, 72049, 72052, #mostly 72's and 2,000's 
#highest var in NDVI: rtes 14059, 14140, 69253, 69021, 85011

elev_ranked = env_all %>% 
  group_by(stateroute) %>% 
  summarize(elev_m = mean(elev.var)) %>%
  arrange(desc(elev_m))
#lowest var in elev: rtes 34027 (best, closest to normal avgs), mostly 34's, 35010, 
#highest var in elev: rtes 17221, 6012, 17044, 6071, 85169, 14059 mostly 14's, 17's, and 6,000's

central_alt = bbs_allscales %>%  
  dplyr::select(logA, pctCore) %>% 
  transmute(pctCore_m = rollapply(pctCore, width = 1, FUN = mean, na.rm = TRUE, fill = NULL),
            logA = logA) %>% 
  mutate(logA = round(logA, digits = 2)) %>%
  group_by(logA) %>%
  summarise(pctCore = mean(pctCore_m)) %>% 
  mutate(focalrte = "99999", logA = logA) %>% 
  dplyr::select(focalrte, logA, pctCore)

bbs_allsub = bbs_allscales %>% 
  filter(focalrte == 34054 | focalrte == 85169) %>%
  dplyr::select(focalrte, logA, pctCore)

bbs_allsub2 = rbind(bbs_allsub, central_alt)

bbs_allsub2$focalrte = factor(bbs_allsub2$focalrte,
                              levels=c( "99999","34054", "85169"),
                              labels=c("Mean",
                                       "Low Heterogeneity",
                                       "High Heterogeneity"))
#use this to assign diff colors for each factor level per what color scheme is ideal?
#72 is PA, 14 is Cali, 34 is Illinois, 17 is Colorado 

pred_plot = ggplot(bbs_allscales, aes(x = logA, y = pctCore))+geom_line(aes(group = focalrte), color = "grey")+
  theme_classic()+geom_abline(aes(intercept = 0.5, slope = 0), linetype = "dashed")+
  geom_line(data = bbs_allsub2, aes(x = logA, y = pctCore, group = as.factor(focalrte), color = as.factor(focalrte)), size = 2)+ #geom_smooth(model = lm, color = 'red')+
  labs(x = expression("Log"[10]*" Area"), y = "", title = "A")+
  scale_color_viridis(discrete = TRUE, name = "", option = "B", begin = 0.05, end = .75)+
  theme(plot.title = element_text(size = 34), axis.title = element_text(size = 30), 
        axis.text = element_text(size = 28, color = "black"), legend.text = element_text(size = 28), 
        legend.title = element_text(size = 30))+
  theme(legend.position = c(0.74, 0.18), legend.key = element_rect(size = 4, color = 'white'),
        legend.key.size = unit(2, 'lines')) 
pred_plot #yellow = high variation in habhet, purple = low variation, low habhet 

central2_alt = bbs_allscales %>%  
  dplyr::select(logN, pctCore) %>% 
  transmute(pctCore_m = rollapply(pctCore, width = 1, FUN = mean, na.rm = TRUE, fill = NULL),
            logN = logN) %>% 
  mutate(logN = round(logN, digits = 1)) %>%
  group_by(logN) %>%
  summarise(pctCore = mean(pctCore_m)) %>% 
  mutate(focalrte = "99999", logN = logN) %>% 
  dplyr::select(focalrte, logN, pctCore)

#calculate confidence intervals for central_alt and central2_alt vals


bbs_allsub = bbs_allscales %>% 
  filter(focalrte == 34054 | focalrte == 85169) %>%
  dplyr::select(focalrte, logN, pctCore)

bbs_allsub3 = rbind(bbs_allsub, central2_alt)

bbs_allsub3$focalrte = factor(bbs_allsub3$focalrte,
                              levels=c( "99999","34054", "85169"),
                              labels=c("Mean",
                                       "Low Heterogeneity",
                                       "High Heterogeneity"))
#use this to assign diff colors for each factor level per what color scheme is ideal?
#72 is PA, 14 is Cali, 34 is Illinois, 17 is Colorado 

pred_abuns = ggplot(bbs_allscales, aes(x = logN, y = pctCore))+geom_line(aes(group = focalrte), color = "grey")+
  theme_classic()+geom_abline(aes(intercept = 0.5, slope = 0), linetype = "dashed")+
  geom_line(data = bbs_allsub3, aes(x = logN, y = pctCore, group = as.factor(focalrte), color = as.factor(focalrte)), size = 2)+ #geom_smooth(model = lm, color = 'red')+
  labs(x = expression("Log"[10]*" Community Size"), y = "", title = "B")+
  scale_color_viridis(discrete = TRUE, name = "", option = "B", begin = 0.05, end = .75)+
  theme(plot.title = element_text(size = 30), axis.title = element_text(size = 30), axis.text = element_text(size = 28, color = "black"), legend.text = element_text(size = 30), legend.title = element_text(size = 30))+
  theme(legend.position = "none") 
pred_abuns #yellow = high variation in habhet, purple = low variation, low habhet 


p1 = gridExtra::grid.arrange(pred_plot, pred_abuns, ncol = 2, 
                             left = ggpubr::text_grob("Proportion Core Species in Community", 
                                                      rot = 90, vjust = 1, size = 30))
ggsave(file = "output/Figure5.tiff", plot = p1)

