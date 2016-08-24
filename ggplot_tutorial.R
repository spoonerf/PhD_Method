library(ggplot2)
library(rgdal)

# read shapefile
wmap <- readOGR(dsn="ne_10m_land.shp", layer="ne_10m_land")

# convert to dataframe
wmap_df<-fortify(wmap)

# create a blank ggplot theme

theme_opts<-list(theme(panel.grid.minor = element_blank(),
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

# plot map

ggplot(wmap_df, aes(long,lat,group=group)) + 
  geom_polygon() + 
  labs(title = "World map (longlat)") +
  coord_equal() +
  theme_opts

ggsave("map1.png", width=12.5, height=8.25, dpi=72)

# reproject from longlat to robinson

wmap_robin <- spTransform(wmap, CRS("+proj=robin"))

wmap_df_robin <- fortify(wmap_robin)
ggplot(wmap_df_robin, aes(long,lat, group=group)) + 
  geom_polygon() + 
  labs(title="World map (robinson)") + 
  coord_equal() +
  theme_opts

ggsave("map2.png", width=12.5, height=8.25, dpi=72) 


# show hole

ggplot(wmap_df_robin, aes(long,lat, group=group, fill=hole)) +
  geom_polygon() + 
  labs(title="World map (robin)") +
  coord_equal() + 
  theme_opts
ggsave("map3.png", width=12.5, height=8.25, dpi=72) 


# change colors

ggplot(wmap_df_robin, aes(long,lat, group=group, fill=hole)) + 
  geom_polygon() + 
  labs(title="World map (Robinson)") + 
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("#262626", "#e6e8ed"), guide="none") # change colors & remove legend

ggsave("map4.png", width=12.5, height=8.25, dpi=72) 


# add graticule and bounding box (longlat)

grat <- readOGR("ne_10m_graticules_15.shp", layer="ne_10m_graticules_15") 
grat_df <- fortify(grat)

bbox <- readOGR("ne_10m_wgs84_bounding_box.shp", layer="ne_10m_wgs84_bounding_box") 
bbox_df<- fortify(bbox)

ggplot(bbox_df, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=wmap_df, aes(long,lat, group=group, fill=hole)) + 
  geom_path(data=grat_df, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  labs(title="World map + graticule (longlat)") + 
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("black", "white"), guide="none") # change colors & remove legend

ggsave("map5.png", width=12.5, height=8.25, dpi=72) 



# graticule (Robin)

grat_robin <- spTransform(grat, CRS("+proj=robin"))  # reproject graticule
grat_df_robin <- fortify(grat_robin)
bbox_robin <- spTransform(bbox, CRS("+proj=robin"))  # reproject bounding box
bbox_robin_df <- fortify(bbox_robin)

ggplot(bbox_robin_df, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=wmap_df_robin, aes(long,lat, group=group, fill=hole)) + 
  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  labs(title="World map (Robinson)") + 
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("black", "white"), guide="none") # change colors & remove legend

ggsave("map6.png", width=12.5, height=8.25, dpi=72) 


# add country borders

countries <- readOGR("ne_10m_admin_0_countries.shp", layer="ne_10m_admin_0_countries") 
countries_robin <- spTransform(countries, CRS("+init=ESRI:54030"))
countries_robin_df <- fortify(countries_robin)

ggplot(bbox_robin_df, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=countries_robin_df, aes(long,lat, group=group, fill=hole)) + 
  geom_path(data=countries_robin_df, aes(long,lat, group=group, fill=hole), color="white", size=0.3) +
  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  labs(title="World map (Robinson)") + 
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("black", "white"), guide="none") # change colors & remove legend

ggsave("map7.png", width=12.5, height=8.25, dpi=72)


# bubble plot
places <- readOGR("ne_10m_populated_places.shp", layer="ne_10m_populated_places") 
places_df <- as(places, "data.frame")
places_robin_df <- project(cbind(places_df$LONGITUDE, places_df$LATITUDE), proj="+init=ESRI:54030") 
places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("LONGITUDE", "LATITUDE")
places_robin_df$POP2000 <- places_df$POP2000 

places_robin_df_sub<-subset(places_robin_df,POP2000>=5000)

ggplot(bbox_robin_df, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=countries_robin_df, aes(long,lat, group=group, fill=hole)) + 
  geom_point(data=places_robin_df_sub, aes(LONGITUDE, LATITUDE, group=NULL, fill=NULL, size=POP2000), color="#32caf6", alpha=I(8/10)) +
  geom_path(data=countries_robin_df, aes(long,lat, group=group, fill=hole), color="white", size=0.3) +
  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  labs(title="World map (Robinson)") + 
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("black", "white"), guide="none")+
  scale_size_continuous(range=c(1,20), guide="none")# change colors & remove legend

ggsave("map8.png", width=12.5, height=8.25, dpi=72) 


# Winkel tripel projection
countries_wintri <- spTransform(countries, CRS("+proj=wintri"))
bbox_wintri <- spTransform(bbox, CRS("+proj=wintri"))
wmap_wintri <- spTransform(wmap, CRS("+proj=wintri"))
grat_wintri <- spTransform(grat, CRS("+proj=wintri"))

p<-ggplot(bbox_wintri, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=countries_wintri, aes(long,lat, group=group, fill=hole)) + 
  geom_path(data=countries_wintri, aes(long,lat, group=group, fill=hole), color="white", size=0.3) +
  geom_path(data=grat_wintri, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  labs(title="World map (Winkel Tripel)") + 
  coord_equal(ratio=1) + 
  theme_opts +
  scale_fill_manual(values=c("black", "white"), guide="none") # change colors & remove legend

ggsave(plot=p, "map9.png", width=12.5, height=8.25, dpi=72)




















