CR40s<-brick("cru_ts3.23.1941.1950.tmp.dat.nc")
CR50s<-brick("cru_ts3.23.1951.1960.tmp.dat.nc")
CR60s<-brick("cru_ts3.23.1961.1970.tmp.dat.nc")
CR70s<-brick("cru_ts3.23.1971.1980.tmp.dat.nc")
CR80s<-brick("cru_ts3.23.1981.1990.tmp.dat.nc")
CR90s<-brick("cru_ts3.23.1991.2000.tmp.dat.nc")
CR00s<-brick("cru_ts3.23.2001.2010.tmp.dat.nc")

obs<-stack(CR40s[[109:120]], CR50s, CR60s, CR70s, CR80s, CR90s,CR00s[[1:48]] )


###85 rcp

ts85_05_30<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_200512-203011.nc")
ts85_30_55<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_203012-205511.nc")
ts85_55_80<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_205512-208011.nc")
ts85_80_00<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_208012-210011.nc")

ts85_05_00<-stack(ts85_05_30[[2:300]],ts85_30_55,ts85_55_80,ts85_80_00[[1:229]])
#ts_05_00<-brick(ts_05_00)

modelcc<-rotate(ts85_05_00)

obs_r<-resample(obs, modelcc)

obs_rk<-obs_r+273.15

mask_mod<-mask(modelcc, obs_rk[[1]])

#obs_mod<-stack(obs_rk, modelcc)

obs_mod<-stack(obs_rk, mask_mod)

obs_mod_mon_mean<-matrix(cellStats(obs_mod, mean))

plot(obs_mod_mon_mean)
##################################################
year_mon<-rep(1:(nlayers(obs_mod)/12), each=12)

#changing from monthly to annual data
obs_mod_ann<-stackApply(obs_mod, indices=year_mon, fun=mean, na.rm=TRUE)

obs_mat_ann<-matrix(cellStats(obs_mod_ann, mean))
plot(obs_mat_ann)

#for lambda_mean
# x = c(1:149)
# g = mgcv:::gam(obs_mat_ann~s(x, k=8), fx=TRUE)
# plot(g)
# 
# g_pred = predict(g)
# g_pred = g_pred[55:148]
# d = diff(g_pred)
# plot(g_pred)
# plot(d)

#for lambda sum
d<-cumsum(diff(obs_mat_ann[55:148]))

#bird values
# centre_tempb<-0.04526127
# scale_tempb<-0.07075669
# 
# #mammal values
# centre_tempm<-0.01479149
# scale_tempm<-0.0793229

###with lambda sum

#bird values
centre_tempb<-
scale_tempb<-

#mammal values
centre_tempm<-
scale_tempm<-

#both values

centre_temp<-0.3605663
scale_temp<-0.6856779

gam_diff<-scale(d, center=centre_temp, scale=scale_temp)
plot(gam_diff)

ave_mod<-data.frame(gam_diff, rep(0, length(gam_diff)), rep(0, length(gam_diff)))
colnames(ave_mod)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred<-predict(m1c, ave_mod, re.form=NA, se.fit=T)



plot(pop_pred$fit*100, type="l", ylim=c(-100, 10))
#plot(pop_pred$fit, ylim=c(-0.045, 0.005))
lines(pop_pred$fit*100 + 1.96*pop_pred$se.fit*100, col="red", lty=2)
lines(pop_pred$fit*100 - 1.96*pop_pred$se.fit*100, col="red", lty=2)
#abline(h=0)


fun<-function(y){
  
  z<-matrix(y)
  if (sum(is.na(z))<1){

    d<-cumsum(diff(z[56:149]))
    gam_diff<-scale(d, center = centre_temp, scale = scale_temp)
    ave_mod<-data.frame(gam_diff, rep(0, length(gam_diff)), rep(0, length(gam_diff)))
    colnames(ave_mod)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
    pop_pred<-predict(m1c, ave_mod, re.form=NA, se.fit=T)
  } else{
    
    #pop_pred$fit<-rep(NA,45) #2050
    pop_pred$fit<-rep(NA,93)
  }
  return(pop_pred$fit)
}

#

pred_rast_2100<-calc(obs_mod_ann, fun)

pred_rast_2100[[93]][pred_rast_2100[[93]] < -1] <- -1

plot(pred_rast_2100[[93]])
pred_rast_2100_mean<-matrix(cellStats(pred_rast_2100, mean))

plot(pred_rast_2100_mean, type="l")








##########################################



10^sum(pop_pred$fit + 1.96*pop_pred$se.fit)
10^sum(pop_pred$fit - 1.96*pop_pred$se.fit)
10^sum(pop_pred$fit)

pop_pred$fit<-c(0, pop_pred$fit)
pop_pred$se.fit<-c(0, pop_pred$se.fit)

plot(y=10^cumsum(pop_pred$fit), x=2006:2099, ylim=c(0, 3), type="l", ylab="Mammal Population Index", xlab="Year", main="Predicted Mammal Population Trends Under RCP 85")
lines(y=10^cumsum(pop_pred$fit - 1.96*pop_pred$se.fit), x=2006:2099, lty=2, col="Red")
lines(y=10^cumsum(pop_pred$fit + 1.96*pop_pred$se.fit), x=2006:2099, lty=2, col="Red")

x = c(1:149)

fun<-function(y){
  
  z<-matrix(y)
  if (sum(is.na(z))<1){
    g<-mgcv:::gam(matrix(z) ~ s(x, k=8), fx=TRUE)
    g_pred<-predict(g)
    g_pred<-g_pred[55:148]
    #g_pred<-g_pred[55:100] #2050
    d<-diff(g_pred)
    gam_diff<-scale(d, center = centre_tempb, scale = scale_tempb)
    ave_mod<-data.frame(gam_diff, rep(0, length(gam_diff)), rep(0, length(gam_diff)))
    colnames(ave_mod)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
    pop_pred<-predict(m1c, ave_mod, re.form=NA, se.fit=T)
  } else{
    
    #pop_pred$fit<-rep(NA,45) #2050
    pop_pred$fit<-rep(NA,94)
  }
  return(pop_pred$fit)
}

gam_pred_rast_2050<-calc(obs_mod_ann, fun)

gam_pred_rast_2100<-calc(obs_mod_ann, fun)

pop_2100<-function(x){10^sum(x)}

pred_2100<-calc(gam_pred_rast_2100, pop_2100)

pred_2050<-calc(gam_pred_rast_2050, pop_2100)

library(RColorBrewer)
library(rgdal)
library(ggplot2)
world<-readOGR(dsn="C:/Users/Fiona/Documents/GIS/ne_10m_land", layer="ne_10m_land")
e<-c(-180, 180, -63, 83.6341)
world_c<-crop(world, e)
plot(world_c)
world_df<-fortify(world_c)
colnames(world_df)[c(1:2)]<-c("Longitude", "Latitude")


pred_2100_df<-rasterToPoints(pred_2100)
pred_2100_df<-rasterToPoints(log10(pred_2100))
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
  scale_fill_gradientn(colors=brewer.pal(7, "RdYlGn"), guide=guide_colourbar(title="Bird Population Change"), 
  #labels=c("-99%", "-50%", "0%", "+50%", "+100%"), breaks=c(0.01,0.5, 1, 1.5, 2))+   
  labels=c("-99.9999%", "-99.99%", "-99%","0%"))+
  geom_path(data=world_df, aes(Longitude, Latitude, group=group))+
  coord_fixed(ratio = 1)+
  theme_opts
  




  brk <- c(1, 0, -1, -2, -3, -4)
plot(log10(pred_2100), col=brewer.pal(7, "RdYlGn"))

#plot(world, lwd=0.1, add=T)

legend("right",inset=F, title="Predicted Population Decline",
       c("0%","90%","99%", "99.9%", "99.99%", "99.999%", "99.9999%"), fill=terrain.colors(7))









plot(log10(pred_2050))
hist(pred_2050)



plot((pred_2100 - 1)*100,col=terrain.colors(12), breaks=brk)



#install.packages("merTools")
library(merTools)
preds <- predictInterval(m1c, newdata = ave_mod, n.sims = 999)



