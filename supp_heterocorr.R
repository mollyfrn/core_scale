#This script performs a supplemental analysis to answer the question: 
# How well does heterogeneity correlate with itself across scales? 
#start with the variance in heterogeneity at the scale of a single stateroute
#even tho analyses in coef_hetero_analyses.R starts with the scale of 3:66

#read in env_all (data of all env variances for each stateroute calc'd across each agg scale pool)
#rearrange df so that: 
# columns = scale 
# rows = stateroute/focal rtes 
# values = variance vals 

#then, corr(df[2:end]) 
#plot as a function of column # e.g. scale so that mirrors figure 7 in text 

#################################################################################
#read in data
env_all = read.csv("intermed/env_all.csv", header = TRUE) 
env_all = env_all %>% dplyr::select(1:4)


#rearrange df: 
ndvi_supp = env_all %>% 
  dplyr::select(stateroute:ndvi.var) %>% 
  #tidyr::unite(hetero = ndvi.var:elev.var) %>% 
  tidyr::spread(key = scale,
                value = ndvi.var)

elev_supp = env_all %>% 
  dplyr::select(-3) %>%
  tidyr::spread(key = scale, 
                value = elev.var)


#preview:
ndvi_corr = round(cor(ndvi_supp[, 2:ncol(ndvi_supp)]), 2) #round to 2 decimal digits
ndvi_corr = as.data.frame(ndvi_corr)
ndvi_corr = ndvi_corr %>% mutate(env_var = "ndvi")
write.csv(ndvi_corr, "intermed/ndvi_corr.csv", row.names = FALSE)

elev_corr = round(cor(elev_supp[, 2:ncol(elev_supp)]), 2) #round to 2 decimal digits
elev_corr = as.data.frame(elev_corr)
elev_corr = elev_corr %>% mutate(env_var = "elev")
write.csv(elev_corr, "intermed/elev_corr.csv", row.names = FALSE)


#rearrange for easy plotting 

ndvi_new = ndvi_corr %>% 
  tidyr::gather(key = scale, 
                value = corr_r, 
                1:66)

elev_new = elev_corr %>% 
  tidyr::gather(key = scale, 
                value = corr_r, 
                1:66)

vars_supp = full_join(elev_new, ndvi_new)
write.csv(vars_supp, "intermed/vars_supp.csv", row.names = FALSE)

####Plotting####

supp_heterocorr = read.csv("intermed/vars_supp.csv", header = TRUE) 

supp_heterocorr$dep = factor(supp_heterocorr$env_var, 
                             levels=c("elev", "ndvi"),
                             labels=c("Elevation", "NDVI"))

#scale on x and r on y, panel by coef of interest, line color by var measure
ggplot(supp_heterocorr, aes(x = scale, y = corr_r))+
  geom_line()+facet_wrap(~dep)+theme_classic()+
  geom_abline(intercept = 0, slope = 0)+
  theme_classic()+theme(text = element_text(size = 18))+
  labs(color = "Environmental Heterogeneity", x = "Number of aggregated BBS Routes", y = "Pearson's correlation coefficient")+theme(legend.position = c(0.84, 0.20))+
  scale_color_viridis(begin = 0, end = 0.7, discrete = TRUE, option = "D")+
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4))
ggsave(file = "output/Figure10.tiff", plot = last_plot())
