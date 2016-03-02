library(raster)
library(plyr)
library(broom)
library(taRifx)
library(zoo)
library(doParallel)
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

colnames(sr)<-yr

srm<-melt(sr)

lon<-rep(xy_df[,1])
lat<-rep(xy_df[,2])

db<-rep(diam, each=479)

srm2<-cbind(lon,lat,db, srm)

srm3<-srm2[,c(1,2,3,5,6)]

colnames(srm3)<-c("Longitude", "Latitude", "Buffer", "Year", "Mean_T")


ID_xy<-data.frame(LPI_EU4$ID, LPI_EU4$Longitude, LPI_EU4$Latitude)

colnames(ID_xy)<-c("ID","Longitude", "Latitude")

srm4<-merge(ID_xy, srm3, by=c("Longitude", "Latitude"))


write.csv(srm4, "Annual_Meant_Buff.csv")

all_EU_mean_data<-read.csv("Annual_Meant_Buff.csv")


######

doMean = function(sp_name) {
  spid2 = subset(LPI_EU4, ID == sp_name)   #subsetting the population data by each population
  spid = spid2[,64:126]                     #subsetting only the dates
  colnames(spid)<-1950:2012              #renaming the date column names as R doesn't like numbered column names
  Year<-1950:2012   
  climid=subset(all_EU_mean_data, ID == sp_name & Buffer == 100000)  #subsetting the climate data by each population
  
  name<-spid2$Binomial
  id<-spid2$ID
  Date<-as.numeric(colnames(spid))
  spidt<-destring(t(spid))
  Yrs<-Year[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  
  Mean_yr<-climid[climid$Year %in% Yrs, ]
  Mean<-Mean_yr$Mean_T
  Yr<-Mean_yr$Year
  
  if (sum(is.na(Mean))!=length(Mean)){
    
    lm_mean<-lm(Mean~Yr)
    lm_mean_df<-tidy(lm_mean)[2,]  
    mean_df<-cbind(id,lm_mean_df)
    
  } else{
    
    mean_df<-matrix(c(id,NA,NA,NA,NA,NA), nrow=1, ncol=6)
    colnames(mean_df)<-c("id", "term", "estimate", "std.error", "statistic", "p.value")
    mean_df<-data.frame(mean_df)
  }
  
  
  print(mean_df)  
  return(mean_df)
}

all_df_list <- lapply(unique(LPI_EU4$ID), doMean)

all_matrix <- matrix(unlist(all_df_list), ncol=6, byrow=TRUE)
mean_df <- data.frame(all_matrix)
colnames(mean_df) <- c("ID", "Term","Estimate","SE","Statistic","p.val")

#clim_sens<-mean_df[,c(1,3,4)]
#colnames(clim_sens)[c(2,3)]<-c("Estimate","SE")

sens100k<-mean_df[,c(1,3)]
colnames(sens100k)[c(2)]<-c("Estimate100k")

sens<-cbind(sens0k, sens10k[,2], sens25k[,2], sens50k[,2], sens100k[,2])
colnames(sens)[2:6]<-c("MnSlope0k", "MnSlope10k", "MnSlope25k", "MnSlope50k", "MnSlope100k")

write.csv(sens, "Climate_Sensitivity_2016_03_02.csv")

