#Sourcing and formatting environmental data for RAMAS
install.packages("landsat")

library(raster)

LPI<-read.csv("LPI_pops_20160523_edited.csv")
alp<-subset(LPI, Common_name=="Alpine ibex" & Specific_location ==1)
#elevation
#elev<-getData('alt', country='ETH', mask=TRUE) #doesn't work HTTP 404 error

country<-getData("GADM", country="CHE", level=0)
elev1<-getData("SRTM", lat=46, lon= 8)
elev2<-getData("SRTM", lat=46, lon= 10)

e<-extent(min(alp$Longitude)-0.5, max(alp$Longitude)+0.5, min(alp$Latitude)-0.5, max(alp$Latitude)+0.5)
elev1c<-crop(elev1,e)
elev2c<-crop(elev2,e)

elev<-merge(elev1c, elev2c)

proj4string(elev) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
plot(elev)

plot(country, add=TRUE)
points(alp$Longitude, alp$Latitude)

elevm<-as.matrix(elev, nrow=6001, ncol=12001)

library(landsat)

asp<-slopeasp(elevm, 0.00083333333, 0.00083333333, smoothing=1)
rasp<-raster(asp$aspect)
extent(rasp)<-e
rerasp<-resample(rasp, elev)
plot(rerasp)
points(alp$Longitude, alp$Latitude)

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
newproj <- "+init=epsg:4326"

lu_tr<-projectRaster(lu,rerasp, crs=CRS(newproj), method="ngb")
plot(lu_tr)



plot(ds,1)
points(alp$Longitude, alp$Latitude)


bioc<-getData("worldclim", var="bio",lat=46, lon= 8, res=0.5)

# bio_5 http://www.tandfonline.com/doi/abs/10.1080/08927014.2004.9522636  #max temp of hottest month -  During the hottest summer females moved over larger ranges at higher altitudes. 
# bio_6 http://www.tandfonline.com/doi/abs/10.1080/08927014.2004.9522636  #min temp of coldest month - In the presence of thick snow cover, females significantly reduced winter home range sizes. 
# bio_10  same paper ^^ - mean temp of warmest quarter 
# bio_11  same paper ^^ - mean temp of coldest quarter 


bio<-bioc[[c(5,6,10,11)]]

biocr<-resample(bio, lu_tr)

ds<-stack(elev, rerasp, lu_tr, biocr)

writeRaster(ds[[1]], "Elevation_Alps", format="RST")
writeRaster(ds[[2]], "Aspect_Alps.asc")
writeRaster(ds[[3]], "Land_Use_Alps.asc")
writeRaster(ds[[4]], "Bio_5_Alps.asc")
writeRaster(ds[[5]], "Bio_6_Alps.asc")
writeRaster(ds[[6]], "Bio_10_Alps.asc")
writeRaster(ds[[5]], "Bio_11_Alps.asc")



