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

days<-nlayers(cropn)    #splitting the data evenly between the cores
step<-floor(days/n)

ptime <- system.time({   
  df<- foreach(1, .combine=cbind) %dopar%{
    rasterex <- raster:::extract(cropn[[days]], xy_df)
  }
}) 
ptime 
beep(3)
stopCluster(cl)

