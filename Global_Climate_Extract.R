library(raster)
library(doParallel)
library(beepr)
library(lubridate)
library(reshape2)

HC<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/HadCRUT.4.4.0.0.median.nc")
CT<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/CRUTEM.4.4.0.0.anomalies.nc")

CR40s<-brick("D:/Fiona/PhD1/Climate/Global/cru_ts3.23.1941.1950.tmp.dat.nc/cru_ts3.23.1941.1950.tmp.dat.nc")
CR50s<-brick("D:/Fiona/PhD1/Climate/Global/cru_ts3.23.1951.1960.tmp.dat.nc/cru_ts3.23.1951.1960.tmp.dat.nc")
CR60s<-brick("D:/Fiona/PhD1/Climate/Global/cru_ts3.23.1961.1970.tmp.dat.nc/cru_ts3.23.1961.1970.tmp.dat.nc")
CR70s<-brick("D:/Fiona/PhD1/Climate/Global/cru_ts3.23.1971.1980.tmp.dat.nc/cru_ts3.23.1971.1980.tmp.dat.nc")
CR80s<-brick("D:/Fiona/PhD1/Climate/Global/cru_ts3.23.1981.1990.tmp.dat.nc/cru_ts3.23.1981.1990.tmp.dat.nc")
CR90s<-brick("D:/Fiona/PhD1/Climate/Global/cru_ts3.23.1991.2000.tmp.dat.nc/cru_ts3.23.1991.2000.tmp.dat.nc")
CR00s<-brick("D:/Fiona/PhD1/Climate/Global/cru_ts3.23.2001.2010.tmp.dat.nc/cru_ts3.23.2001.2010.tmp.dat.nc")
CR10s<-brick("D:/Fiona/PhD1/Climate/Global/cru_ts3.23.2011.2014.tmp.dat.nc/cru_ts3.23.2011.2014.tmp.dat.nc")

plot(CR40s[[11]])

LPI<-read.csv("D:/Fiona/Git_Method/Git_Method/LPI_populations_IP_fishedit_20140310_nonconf.csv")

LPIsp<-subset(LPI, Specific_location==1 & System !="Marine" & Class != "Actinopterygii"& Class != "Cephalaspidomorphi" )

CR<-stack(CR40s,CR50s,CR60s,CR70s,CR80s,CR90s,CR00s,CR10s)

plot(CR[[1]])
points(LPIsp$Longitude, LPIsp$Latitude)

xy<-data.frame(LPIsp$Longitude, LPIsp$Latitude)

test<-extract(CR[[1]], xy, buffer=5000, na.rm=T)


ex<-matrix(unlist(test), ncol=1, byrow = TRUE)

cbind(ex,xy)

xy<-unique(xy)     #identifying unique locations to extract climate data from 

xy_df<-data.frame(xy)
colnames(xy_df)<-c("lon", "lat")
coordinates(xy_df) <- c("lon", "lat")

length(xy_df)


###################################

CR2<-CR[[1:12]]

n<-6#number of cores to use - not sure how many I can go up to
cl<-makeCluster(n)
registerDoParallel(cl)  
getDoParWorkers()

buff<- function (year,b)  {
  
  rasterex<-raster:::extract(CR[[year]], xy_df, buff=b, fun=mean, na.rm=TRUE)
  return(rasterex)
  
}

# yr<-rep(1950:1952,each=2)
# diam<-c(0,10000,25000)

#month<-seq(1941, 2014.917, by= 1/12)
diam<-c(50000)
lyr<-1:nlayers(CR)


stime <- system.time({
  sr <- foreach(1, .combine = cbind) %dopar% mapply(buff,lyr,diam)
})
stime
beep(3)

stopCluster(cl)


dates<-seq(ymd('1941-01-16'),ymd('2014-12-16'), by = 'months')


datesr<-rep(dates, each=1078)

srm<-melt(sr)

lon<-xy[,1]
lat<-xy[,2]

srm2<-cbind(lon,lat,datesr, srm)

srm3<-srm2[,c(1,2,3,6)]

colnames(srm3)<-c("Longitude", "Latitude", "Date", "Mean_T")

LPI_ID<-LPIsp[,c("ID", "Longitude", "Latitude")]

LPIclim<-merge(srm3, LPI_ID, by=c("Longitude", "Latitude"))

write.csv(LPIclim, "Global_Mean_Temp_All_LPI.csv")















