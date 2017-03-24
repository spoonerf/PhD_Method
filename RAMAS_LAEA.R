#Sourcing and formatting environmental data for RAMAS
install.packages("landsat")

library(raster)

LPI<-read.csv("LPI_pops_20160523_edited.csv")
alp<-subset(LPI, Common_name=="Alpine ibex" & Specific_location ==1)

xy<-cbind(alp$Longitude, alp$Latitude)
xy<-SpatialPoints(xy)
proj4string(xy)<-CRS("+init=epsg:4326")
xy<-spTransform(xy, CRS("+init=epsg:3035"))


lu<-raster("eu27ch2010.asc")
plot(lu)

# 111 - Settlement
# 222 - Cropland
# 333 - Forest
# 444 - Grassland
# 555 - Other Land
# 666 - Water

library(rgdal)
proj4string(lu) <- CRS("+init=epsg:3035")
plot(lu)
points(xy)

xmn<-xmin(xy) - 10000
ymn<-ymin(xy) - 10000

xmx<-xmax(xy) + 10000
ymx<-ymax(xy) + 10000

e<-extent(xmn,xmx,ymn,ymx)

luc<-crop(lu,e)

elev1<-getData("SRTM", lat=46, lon= 8)
elev2<-getData("SRTM", lat=46, lon= 10)

proj4string(elev1)<-CRS("+init=epsg:4326")
proj4string(elev2)<-CRS("+init=epsg:4326")

newproj<- "+init=epsg:3035"

elev1r<-projectRaster(elev1,luc, crs=CRS(newproj), method="bilinear")
elev2r<-projectRaster(elev2,luc, crs=CRS(newproj), method="bilinear")

elev<-merge(elev1r, elev2c)

proj4string(elev) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
plot(elev)


elevm<-as.matrix(elev, nrow=266, ncol=362)

library(landsat)

asp<-slopeasp(elevm, 1000, 1000, smoothing=1)
rasp<-raster(asp$aspect)
extent(rasp)<-e
plot(rasp)


bioc<-getData("worldclim", var="bio",lat=46, lon= 8, res=0.5)

# bio_5 http://www.tandfonline.com/doi/abs/10.1080/08927014.2004.9522636  #max temp of hottest month -  During the hottest summer females moved over larger ranges at higher altitudes. 
# bio_6 http://www.tandfonline.com/doi/abs/10.1080/08927014.2004.9522636  #min temp of coldest month - In the presence of thick snow cover, females significantly reduced winter home range sizes. 
# bio_10  same paper ^^ - mean temp of warmest quarter 
# bio_11  same paper ^^ - mean temp of coldest quarter 


bio<-bioc[[c(5,6,10,11)]]

biocr<-projectRaster(bio,luc, crs=CRS(newproj), method="bilinear")


biocr<-resample(bio, luc)

ds<-stack(elev, luc, biocr)

writeRaster(ds[[1]], "Elevation_Alps", format="RST", overwrite=T)
writeRaster(ds[[2]], "Land_Use_Alps", format="RST", overwrite=T)
writeRaster(ds[[3]], "Bio_5_Alps", format="RST", overwrite=T)
writeRaster(ds[[4]], "Bio_6_Alps", format="RST", overwrite=T)
writeRaster(ds[[5]], "Bio_10_Alps", format="RST", overwrite=T)
writeRaster(ds[[6]], "Bio_11_Alps", format="RST", overwrite=T)











