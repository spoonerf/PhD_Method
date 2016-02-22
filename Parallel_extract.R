library(doParallel)

##adding in large raster brick - about 4GB
rmax<-brick("D:/Fiona/Git_Method/Git_Method/tx_0.25deg_reg_v11.0.nc", varname = "tx")

##adding in dataframe we want to extract from and subsetting it so that only relevant locations are left
LPI<-read.csv("D:/Fiona/Git_Method/Git_Method/LPI_populations_IP_fishedit_20140310_nonconf.csv", stringsAsFactors=FALSE)
LPI_EU<-subset(LPI, ClassX != "Fishes"& Specific_location =="1")
ID<-LPI_EU$ID
pop_data<- LPI_EU[,c(1,63:125)]
pop_datab <- (pop_data [,2:64] !="NULL")
points_per_pop1950_2012 = rowSums(pop_datab)
length_id <- data.frame(ID,points_per_pop1950_2012)
LPI_EU2<-merge(length_id, LPI_EU, by = "ID")
LPI_EU3<-subset(LPI_EU2, points_per_pop1950_2012 >=2)    #disregarding any populations with less than two time points
LPI_EU4<-subset(LPI_EU3, Latitude >= 25.5 & Latitude <= 75.5 & Longitude >= -40.5 & Longitude <= 75.5) 

xy<-cbind(LPI_EU4$Longitude,LPI_EU4$Latitude)     #removing points which are in the sea/don't have any data on the most recent day
test<-extract(rmax[[23741]], xy, buffer=50000, fun=mean, na.rm=TRUE)
df<-data.frame(id, test)
df2<-subset(df, !is.na(test))

LPI_EU5<-LPI_EU4[LPI_EU4$ID %in% df2$id,]

xy<-cbind(LPI_EU5$Longitude,LPI_EU5$Latitude)
xy<-unique(xy)     #identifying unique locations to extract climate data from 

###parallellise

n<-6  #number of cores to use - not sure how many I can go up to
cl<-makeCluster(n)
registerDoParallel(cl)  

days<-nlayers(rmax)    #splitting the data evenly between the cores
step<-floor(days/n)

ptime <- system.time({   
  df<- foreach(lyr=seq(1,days,step)[1:6], .combine=cbind) %dopar%{
    rasterex <- raster:::extract(rmax[[lyr:(lyr+step-1)]], xy, buffer=50000, fun=mean, na.rm=TRUE)
  }
}) 
ptime  

stopCluster(cl)



