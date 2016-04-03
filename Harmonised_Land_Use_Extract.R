library(raster)

crop<-brick("D:/Fiona/PhD1/Land_Use/Harmonised/LUHa_u2t1.v1_gcrop.nc4")
prim<-brick("D:/Fiona/PhD1/Land_Use/Harmonised/LUHa_u2t1.v1_gothr.nc4")
past<-brick("D:/Fiona/PhD1/Land_Use/Harmonised/LUHa_u2t1.v1_gpast.nc4")
secd<-brick("D:/Fiona/PhD1/Land_Use/Harmonised/LUHa_u2t1.v1_gsecd.nc4")
urban<-brick("D:/Fiona/PhD1/Land_Use/Harmonised/LUHa_u2t1.v1_gurbn.nc4")


cropn<-crop[[251:306]]
primn<-prim[[251:306]]
pastn<-past[[251:306]]
secdn<-secd[[251:306]]
urbann<-urban[[251:306]]

writeRaster(cropn, "Crop_1950.tif")
writeRaster(primn, "Prim_1950.tif")
writeRaster(pastn, "Past_1950.tif")
writeRaster(secdn, "Secd_1950.tif")
writeRaster(urbann, "Urbn_1950.tif")

cropn<-brick("Crop_1950.tif")
primn<-brick("Prim_1950.tif")
pastn<-brick("Past_1950.tif")
secdn<-brick("Secd_1950.tif")
urbn<-brick("Urbn_1950.tif")


LPI<-read.csv("D:/Fiona/Git_Method/Git_Method/LPI_populations_IP_fishedit_20140310_nonconf.csv")

LPIsp<-subset(LPI, Specific_location==1 & System !="Marine" & Class != "Actinopterygii"& Class != "Cephalaspidomorphi" )

xy<-data.frame(LPIsp$Longitude, LPIsp$Latitude)

xy<-unique(xy)     #identifying unique locations to extract climate data from 

xy_df<-data.frame(xy)
colnames(xy_df)<-c("lon", "lat")
coordinates(xy_df) <- c("lon", "lat")

length(xy_df)


n<-6  #number of cores to use - not sure how many I can go up to
cl<-makeCluster(n)
registerDoParallel(cl)  

days<-nlayers(urbn)    #splitting the data evenly between the cores
step<-floor(days/n)

ptime <- system.time({   
  df<- foreach(days, .combine=cbind) %dopar%{
    rasterex <- raster:::extract(urbn[[1:days]], xy_df)
  }
}) 
ptime 
beep(1)
stopCluster(cl)

dates<-1950:2005

datesr<-rep(dates, each=1078)

dfm<-melt(df)

lon<-xy[,1]
lat<-xy[,2]

dfm2<-cbind(lon,lat,datesr, dfm)

#crop2<-dfm2
#prim2<-dfm2
#past2<-dfm2
#secdn2<-dfm2
#urbn2<-dfm2

land_use<-data.frame(crop2$lon, crop2$lat, crop2$datesr, crop2$value, prim2$value, past2$value, secdn2$value, urbn2$value)

colnames(land_use)<-c("Longitude", "Latitude","Year", "Crop", "Primary", "Pasture", "Secondary", "Urban")


LPI_ID<-LPIsp[,c("ID", "Longitude", "Latitude")]


LPILU<-merge(LPI_ID,land_use, by=c("Longitude", "Latitude"))



test<-subset(LPILU, ID==5485)

plot(test$Year, test$Crop, ylim=c(0,1))
lines(test$Year,test$Primary, type="l", col="green")
points(test$Year,test$Pasture, col="brown")
points(test$Year,test$Secondary, col="hot pink")
points(test$Year,test$Urban, col="grey")


sums<-rowSums(LPILU[,c(5:9)])
hist(sums)

LPILU$Other<- 1 - rowSums(LPILU[,c(5:9)])
 
write.csv(LPILU, "Global_Land_Use_All_LPI.csv")
  
