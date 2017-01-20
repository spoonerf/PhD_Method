mtc<-raster("Global_Rate_Mean_Temp_Change.tif")
luc<-raster("Global_Rate_Land_Use_Change.tif")


plot(mtc)
plot(luc)

centre_tempb<-0.04526127
scale_tempb<-0.07006285

centre_lucb<- -0.0001411668
scale_lucb<- 0.005582236

mtc_s<-scale(mtc, scale=scale_tempb, center=centre_tempb)
luc_s<-scale(luc, scale=scale_lucb, center=centre_lucb)

luc_sr<-resample(luc_s, mtc_s, method="bilinear")

body<-(mtc_s- mtc_s)

stack_pred<-stack(mtc_s, luc_sr, body)

names(stack_pred)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")

pred_rast<-predict(stack_pred,mav, re.form=NA)

plot(pred_rast)

pred_pcnt<-(((10^pred_rast) - 1))

plot(pred_pcnt)        
#plot(mtc)

#Year 2100 rates
#1.5
rate_1.5<-0.000895 #based on 2005 anomaly of 0.65
rate_2.0<-0.01421
model_rate<-0.0701
rate_x<-rate_1.5/model_rate
rate_x2<-rate_2.0/model_rate

extrap<-function(x){
  
  percent<-rate_x*x
  pop_perc_rem<-(1 + percent)^95
  return(pop_perc_rem)
  }

perc_change<-calc(pred_pcnt,extrap)


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


        