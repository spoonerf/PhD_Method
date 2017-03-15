library(raster)

mtc<-raster("Global_Rate_Mean_Temp_Change.tif")
luc<-raster("Global_Rate_Land_Use_Change.tif")


plot(mtc)
plot(luc)

centre_tempb<-0.04526127
scale_tempb<-0.07006285

centre_lucb<- -0.0001411668
scale_lucb<- 0.005582236

centre_bmassb<-2.349304967 
scale_bmassb<-0.916571082

####with no vultures and no realm missing

centre_tempb<-0.045383714
scale_tempb<-0.070236640

centre_lucb<- -0.0001416388
scale_lucb<- 0.005613524

centre_bmassb<-2.328359012 
scale_bmassb<-0.908937973



mtc_s<-scale(mtc, scale=scale_tempb, center=centre_tempb)
luc_s<-scale(luc, scale=scale_lucb, center=centre_lucb)
body_s<-scale(body, scale=scale_bmassb, center=centre_bmassb)

#luc_sr<-resample(luc_s, mtc_s, method="bilinear")

mtc_s<-resample(mtc_s,luc_s, method="bilinear")
luc_sr<-luc_s

body<-(mtc_s - mtc_s)

body_s<-body + log10(1)
body_s = body*0

luc_sr = luc_sr*0

rlm<-readOGR(dsn="C:/Users/Fiona/Documents/GIS/Ecoregions", layer ="ecoregions_dissolved_realm" )
rp <- rasterize(rlm, luc_s, 'WWF_REALM2')


plot(rp)
Afrotropic<-rp
Antarctic<-rp
Australasia<-rp
IndoMalay<-rp
Nearctic<-rp
Neotropic<-rp
Oceania<-rp
Palearctic<-rp

Afrotropic[values(Afrotropic) != 1] <- NA
Antarctic[values(Antarctic) != 2] <- NA
Australasia[values(Australasia) != 3] <- NA
IndoMalay[values(IndoMalay) != 4] <- NA
Nearctic[values(Nearctic) != 5] <- NA
Neotropic[values(Neotropic) != 6] <- NA
Oceania[values(Oceania) != 7] <- NA
Palearctic[values(Palearctic) != 8] <- NA

v = values(luc_s)

realm_mod<-function(rlm, name){
  v = as.numeric(values(rlm))
  realm = v*NA
  realm[!is.na(values(rlm))] = name
  realm = factor(realm, levels(dt$WWF_REALM2))
  stack_pred<-stack(mtc_s, luc_s, body_s)
  temp = as.data.frame(stack_pred)
  temp = cbind(temp, realm)
  colnames(temp)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale", "WWF_REALM2")
  pred_rast<-predict(newdata=temp,m1cr,re.form=NA)
  e<-extent(luc_s)
  pred_rast = raster(matrix(pred_rast, ncol=4320, nrow=2160, byrow=T))
  extent(pred_rast)<-e
  return(pred_rast)
}


Afr<-realm_mod(Afrotropic, "Afrotropic")
Aus<-realm_mod(Australasia, "Australasia")
IdM<-realm_mod(IndoMalay, "Indo-Malay")
Nea<-realm_mod(Nearctic, "Nearctic")
Neo<-realm_mod(Neotropic, "Neotropic")
Oce<-realm_mod(Oceania, "Oceania")
Pal<-realm_mod(Palearctic, "Palearctic")



wrld<-merge(Afr,Aus, IdM, Nea, Neo, Oce, Pal)
plot(wrld)


levelplot(10^(wrld*40), zscaleLog=10, pretty=T,par.settings=GrTheme())

## Change the color theme
levelplot(r, par.settings=GrTheme())
levelplot(r, par.settings=PuOrTheme())

myTheme=rasterTheme(region=brewer.pal('Blues', n=9))
levelplot(r, par.settings=myTheme)





############################################


plot(mtc_s)

luc_r<-resample(luc, mtc, method="bilinear")
luc_sr<-resample(luc_s, mtc_s, method="bilinear")

stack_predus<-stack(mtc, luc_r, body)
stack_pred<-stack(mtc_s, luc_sr, body)

names(stack_predus)<- c("Estimate", "both_change", "Boydmass")
names(stack_pred)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")


predus_rast<-predict(stack_predus, m1cus, re.form=NA)
pred_rast<-predict(stack_pred,m1c, re.form=NA)


predus_pcnt<-(((10^predus_rast) - 1))
pred_pcnt<-(((10^pred_rast) - 1))

plot(predus_pcnt)
plot(pred_pcnt)        
#plot(mtc)

#Year 2100 rates
#1.5
rate_1.5<-0.000895 #based on 2005 anomaly of 0.65, so annual rate of increase to reach 1.5 deg by 2100
rate_2.0<-0.01421
model_rate<-0.0701

rate_x<-rate_1.5/model_rate
rate_x2<-rate_2.0/model_rate

extrap<-function(x){
  
  percent<-rate_x2*x
  pop_perc_rem<-(1 + percent)^95
  return(pop_perc_rem)
  }

perc_change<-calc(predus_pcnt,extrap)
plot(perc_change)



check_stack<-stack(perc_change, pred_pcnt,pred_rast, stack_pred)

cellStats(perc_change, mean)
writeRaster(pred_pcnt, "Predicted_Population_Declines.tif")


library(RColorBrewer)
library(rgdal)
library(ggplot2)
library(rgeos)
world<-readOGR(dsn=getwd(), layer="ne_10m_land")
e<-c(-180, 180, -63, 83.6341)
world_c<-crop(world, e)
plot(world_c)
world_df<-fortify(world_c)
colnames(world_df)[c(1:2)]<-c("Longitude", "Latitude")


pred_2100_df<-rasterToPoints(pred_pcnt)
#pred_2100_df<-rasterToPoints(log10(pred_2100))
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
                       axis.title.y = element_blank()))

ggplot(data=pred_2100_df2, aes(Longitude, Latitude))+
  geom_raster(aes(fill=Population_Decline))+
  scale_fill_gradientn(colors=brewer.pal(7, "RdYlGn"), guide=guide_colourbar(title=""))+ 
  #labels=c("-99%", "-50%", "0%", "+50%", "+100%"), breaks=c(0.01,0.5, 1, 1.5, 2))+   
  #labels=c("-99.9999%", "-99.99%", "-99%","0%"))+
  geom_path(data=world_df, aes(Longitude, Latitude, group=group))+
  coord_fixed(ratio = 1)+
  theme_bw()+
  labs(title="Rate of Land Use Change (Degree Celsius per Year)")+
  theme(text = element_text(size=20))









v = values(luc_sr)
realm = v*NA
realm[!is.na(values(luc_sr))] = "Afrotropic"
realm = factor(realm, levels(dt$WWF_REALM2))

#test_realm = raster(matrix(as.numeric(realm), ncol=720, nrow=360, byrow=T))
#plot(test_realm)

