library(rgdal)
library(ggplot2)

wmap <- readOGR(dsn="ne_110m_land.shp", layer="ne_110m_land")
#countries <- readOGR("ne_50m_admin_0_countries.shp", layer="ne_50m_admin_0_countries")
bbox<-readOGR("ne_10m_wgs84_bounding_box.shp", layer="ne_10m_wgs84_bounding_box") 

temp<-read.csv("All_LPI_All_Years_Nobuff_1931_moreLPI_end2005.csv")
luc<-read.csv("Hyde_crop_pasture_annual_change.csv")
LPI<-read.csv("LPI_pops_20160523_edited.csv")

body<-read.csv("bird_and_mammal_traits2.csv")
body2<-read.csv("LPI_traits.csv")

body3<-rbind(body, body2)

body4<-unique(body3[,c(2:5)])

pop<-read.csv("Global_Population_Trends_Rsq_Lambda_07_10_16.csv")

temp<-temp[,c("ID", "Estimate")]

LPI<-LPI[,c("ID","Binomial","Common_name","Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat")]

df<-merge(merge(temp,luc, by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

df<-merge(df, body4, by="ID")

nrow(df)

df2<-subset(df, !is.na(Estimate)&r_sq >= 0.499999  & !is.na(both_change) & !is.na(Bodymass_g) &
            length_time >=5 & System!="Marine" &Specific_location == 1 &(Class=="Mammalia"|Class=="Aves"))

nrow(df2)

library(plyr)
#counting duplicates at each location
sp_dups<-data.frame(ddply(df2,.(Longitude,Latitude),nrow))
sp_dups$loc_id<-1:length(sp_dups$Longitude)
sp_dups_df<-merge(sp_dups, df2, by=c("Longitude","Latitude"))



loc<-data.frame(sp_dups_df$Longitude,sp_dups_df$Latitude,sp_dups_df$V1)
head(loc)
loc<-unique(loc)
colnames(loc)<-c("Longitude", "Latitude", "V1")

coordinates(loc)<-c("Longitude","Latitude")
proj4string(loc) <- CRS("+proj=longlat")

# loc_wmerc<-spTransform(loc, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
# wmap_wmerc<-spTransform(wmap, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
# countries_wmerc<-spTransform(countries, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))

wmap_df<-fortify(wmap)
#countries_df<-fortify(countries)
loc_df<-data.frame(loc)
bbox_df<-fortify(bbox)
# wmap_wmerc_df<-fortify(wmap_wmerc)
# countries_wmerc_df<-fortify(countries_wmerc)
# loc_wmerc_df<-data.frame(loc_df)

# theme_opts<-list(theme(panel.grid.minor = element_blank(),
#                        panel.grid.major = element_blank(),
#                        panel.background = element_rect(fill = 'light blue', colour = NA),
#                        plot.background = element_rect(fill="light grey",
#                                                       size=1,linetype="solid",color="black"),
#                        axis.line = element_blank(),
#                        axis.text.x = element_blank(),
#                        axis.text.y = element_blank(),
#                        axis.ticks = element_blank(),
#                        axis.title.x = element_blank(),
#                        axis.title.y = element_blank(),
#                        plot.title = element_text(size=22)))

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_rect(fill="#e6e8ed"),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_text(size=22)))



ggplot(data=wmap_df, aes(long,lat, group=group, fill=hole)) +       #bbox_df, aes(long,lat, group=group)
  geom_polygon()+
  geom_point(data=loc_df, aes(Longitude, Latitude, group=NULL,fill=NULL,size=V1),
             color="orange",alpha=I(7/10)) +
  scale_size(range=c(1,7), guide = "legend",labs(size="No. of Populations"))+
  #scale_fill_manual(values=c("#262626", "#e6e8ed"), guide="none")+
  scale_fill_manual(values=c("#262626", "#e6e8ed"), guide="none")+
  theme_opts

#data=wmap_df, aes(long,lat,group=group), fill="black"
###web mercator projection
# ggplot() + 
#   geom_polygon(data=wmap_wmerc_df, aes(long,lat,group=group), fill="white")+
#   geom_point(data=loc_wmerc_df, aes(Longitude, Latitude, group=NULL,fill=NULL,size=V1),
#              color="black",alpha=I(6/10)) +
#   scale_size(range=c(1,7), guide = "legend",labs(size="Count")) +
#   coord_equal() +
#   theme(aspect.ratio=1)+
#   theme_opts







