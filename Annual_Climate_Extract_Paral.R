library(raster)
library(ncdf4)

LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv", stringsAsFactors=FALSE)

LPI_EU<-subset(LPI, ClassX != "Fishes"& Specific_location =="1")

ID<-LPI_EU$ID
pop_data<- LPI_EU[,c(1,63:125)]

pop_datab <- (pop_data [,2:64] !="NULL")
points_per_pop1950_2012 = rowSums(pop_datab)
length_id <- data.frame(ID,points_per_pop1950_2012)

LPI_EU2<-merge(length_id, LPI_EU, by = "ID")
LPI_EU3<-subset(LPI_EU2, points_per_pop1950_2012 >=2)    #disregarding any populations with less than two time points

LPI_EU4<-subset(LPI_EU3, Latitude >= 25.5 & Latitude <= 75.5 & Longitude >= -40.5 & Longitude <= 75.5) 

id<-LPI_EU4$ID
xy<-cbind(LPI_EU4$Longitude,LPI_EU4$Latitude)

xy_df<-unique(xy)

rmean<-brick("tg_0.25deg_reg_v11.0.nc", varname = "tg")

tm <- seq(as.Date('1950-01-01'), as.Date('2014-12-31'), 'day')
s <- setZ(rmean, tm, 'days')

tmy <- seq(as.Date('1950-01-01'), as.Date('2014-12-31'), 'year')

rmn<-zApply(s, tmy, fun=mean, name='year')

#rmn<-brick("Annual_Meant.tif")

n<-8 #number of cores to use - not sure how many I can go up to
cl<-makeCluster(n)
registerDoParallel(cl)  
getDoParWorkers()

buff<- function (year,b)  {
  
  rasterex<-raster:::extract(rmn[[year-1949]], xy, buff=b, fun=mean, na.rm=TRUE)
  return(rasterex)
  
}

# yr<-rep(1950:1952,each=2)
# diam<-c(0,10000,25000)

yr<-rep(1950:2013, each=5)
diam<-c(0,10000,25000,50000,100000)

stime <- system.time({
  sr <- foreach(1, .combine = cbind) %dopar% mapply(buff,yr,diam)
})
stime

stopCluster(cl)







