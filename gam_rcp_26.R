CR40s<-brick("cru_ts3.23.1941.1950.tmp.dat.nc")
CR50s<-brick("cru_ts3.23.1951.1960.tmp.dat.nc")
CR60s<-brick("cru_ts3.23.1961.1970.tmp.dat.nc")
CR70s<-brick("cru_ts3.23.1971.1980.tmp.dat.nc")
CR80s<-brick("cru_ts3.23.1981.1990.tmp.dat.nc")
CR90s<-brick("cru_ts3.23.1991.2000.tmp.dat.nc")
CR00s<-brick("cru_ts3.23.2001.2010.tmp.dat.nc")

obs<-stack(CR40s[[109:120]], CR50s, CR60s, CR70s, CR80s, CR90s,CR00s[[1:48]] )


###26 rcp

ts26_05_30<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_rcp26_r1i1p1_200512-203011.nc")
ts26_30_55<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_rcp26_r1i1p1_203012-205511.nc")
ts26_55_80<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_rcp26_r1i1p1_205512-208011.nc")
ts26_80_99<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_rcp26_r1i1p1_208012-209911.nc")
ts26_99_24<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_rcp26_r1i1p1_209912-212411.nc")

ts26_05_00<-stack(ts26_05_30[[2:300]],ts26_30_55,ts26_55_80,ts26_80_99,ts26_99_24[[1]])

modelcc_26<-rotate(ts26_05_00)

obs_r_26<-resample(obs, modelcc_26)

obs_rk_26<-obs_r_26+273.15

mask_mod_26<-mask(modelcc_26, obs_rk_26[[1]])

#obs_mod<-stack(obs_rk, modelcc)

obs_mod_26<-stack(obs_rk_26, mask_mod_26)

obs_mod_mon_mean_26<-matrix(cellStats(obs_mod_26, mean))

plot(obs_mod_mon_mean_26)
##################################################
year_mon<-rep(1:(nlayers(obs_mod_26)/12), each=12)

#changing from monthly to annual data
obs_mod_ann_26<-stackApply(obs_mod_26, indices=year_mon, fun=mean, na.rm=TRUE)


obs_mat_ann_26<-matrix(cellStats(obs_mod_ann_26, mean))
plot(obs_mat_ann_26)

x = c(1:149)
g = mgcv:::gam(obs_mat_ann_26~s(x, k=8), fx=TRUE)
plot(g)


g_pred = predict(g)
g_pred = g_pred[55:148]
d = diff(g_pred)
plot(d)

#bird values
centre_tempb<-0.04526127
scale_tempb<-0.07075669


#mammal values
centre_tempm<-0.01479149
scale_tempm<-0.0793229


gam_diff<-scale(d, center=centre_tempb, scale=scale_tempb)

ave_mod<-data.frame(gam_diff, rep(0, length(gam_diff)), rep(0, length(gam_diff)))
colnames(ave_mod)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred<-predict(m1c, ave_mod, re.form=NA, se.fit=T)

plot(pop_pred$fit, ylim=c(-0.04, 0.01))
lines(pop_pred$fit + 1.96*pop_pred$se.fit)
lines(pop_pred$fit - 1.96*pop_pred$se.fit)
10^sum(pop_pred$fit + 1.96*pop_pred$se.fit)
10^sum(pop_pred$fit - 1.96*pop_pred$se.fit)
10^sum(pop_pred$fit)

# pop_pred$fit<-c(0, pop_pred$fit)
# pop_pred$se.fit<-c(0, pop_pred$se.fit)

plot(y=10^cumsum(pop_pred$fit), x=2006:2098, ylim=c(0, 4.5), type="l", ylab="Bird Population Index", xlab="Year", main="Predicted Bird Population Trends Under RCP 26")
lines(y=10^cumsum(pop_pred$fit - 1.96*pop_pred$se.fit), x=2006:2098, lty=2, col="red")
lines(y=10^cumsum(pop_pred$fit + 1.96*pop_pred$se.fit), x=2006:2098, lty=2, col="red")



#install.packages("merTools")
library(merTools)
preds <- predictInterval(m1c, newdata = ave_mod, n.sims = 999)


x = c(1:149)

fun<-function(y){
  
  z<-matrix(y)
  if (sum(is.na(z))<1){
    g<-mgcv:::gam(matrix(z) ~ s(x, k=8), fx=TRUE)
    g_pred<-predict(g)
    #g_pred<-g_pred[55:148]
    g_pred<-g_pred[55:100] #2050
    d<-diff(g_pred)
    gam_diff<-scale(d, center = centre_tempb, scale = scale_tempb)
    ave_mod<-data.frame(gam_diff, rep(0, length(gam_diff)), rep(0, length(gam_diff)))
    colnames(ave_mod)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
    pop_pred<-predict(m1c, ave_mod, re.form=NA, se.fit=T)
  } else{
    
    pop_pred$fit<-rep(NA,45) #2050
    #pop_pred$fit<-rep(NA,94)
  }
  return(pop_pred$fit)
}

gam_pred_rast_2050<-calc(obs_mod_ann_26, fun)

pop_2100<-function(x){10^sum(x)}

pop_2050<-function(x){10^sum(x)}
# year<-2007:2099
# png(file="%04d_pop_prod_rcp26.png", width=1440, height=720)
# for (i in 2:nlayers(gam_pred_rast)){
#   
#   gam_stack<-10^prod(gam_pred_rast[[1:i]])
#   plot(gam_stack,axes=FALSE,legend=FALSE, main=year[i])
#   #writeRaster(gam_stack,paste(year[i], "pop_product_rcp26.tif", sep = "_") )
#   print(i)
# }
# dev.off()

pred_2100<-calc(gam_pred_rast, pop_2100)

pred_2050<-calc(gam_pred_rast_2050, pop_2050)

plot(pred_2100)

library(RColorBrewer)
library(rgdal)
library(ggplot2)
world<-readOGR(dsn="C:/Users/Fiona/Documents/GIS/ne_10m_land", layer="ne_10m_land")
e<-c(-180, 180, -63, 83.6341)
world_c<-crop(world, e)
plot(world_c)
world_df<-fortify(world_c)
colnames(world_df)[c(1:2)]<-c("Longitude", "Latitude")


pred_2100_df<-rasterToPoints(pred_2050)
pred_2100_df<-rasterToPoints(log10(pred_2050))
pred_2100_df2<-data.frame(pred_2100_df)
colnames(pred_2100_df2)<-c("Longitude", "Latitude", "Population_Decline")


theme_opts<-list(theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_rect(fill = 'white', colour = NA),
                       plot.background = element_rect(),
                       axis.line = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       plot.title = element_text()))

ggplot(data=pred_2100_df2, aes(Longitude, Latitude))+
  geom_raster(aes(fill=Population_Decline))+
  scale_fill_gradientn(colors=brewer.pal(7, "RdYlGn"), guide=guide_colourbar(title="Bird Population Change"))+#,
                       #labels=c("-99%", "-50%", "0%", "+50%", "+100%"), breaks=c(0.01,0.5, 1, 1.5, 2))+   
                      #labels=c("-99.9%", "-99%", "-90%","0%"), breaks=c(-3, -2, -1, 0))+
  geom_path(data=world_df, aes(Longitude, Latitude, group=group))+
  coord_fixed(ratio = 1)+
  theme_opts




