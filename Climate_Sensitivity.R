library(raster)
library(plyr)
library(broom)
library(taRifx)
library(zoo)
library(doParallel)

LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")
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
Binomial<-LPI_EU4$Binomial
xy<-cbind(LPI_EU4$Longitude,LPI_EU4$Latitude)

xy<-unique(xy)

xmean<- brick("Annual_Meant.tif")
names(xmean)<-seq(1950,2014,1)
xmean<-xmean[[1:64]]

all_EU_mean_data<-data.frame(lon=numeric(0), lat=character(0), Year=numeric(0), ex=numeric(0),ex10k=numeric(0),ex25k=numeric(0),ex50k=numeric(0),ex75k=numeric(0),ex100k=numeric(0))

for (i in 1:nlayers(xmean)) {
  
  ex<- extract(xmean[[i]], xy)
  ex10k<- extract(xmean[[i]], xy, buffer=10000, fun=mean, na.rm=TRUE)
  ex25k<- extract(xmean[[i]], xy, buffer=25000, fun=mean, na.rm=TRUE)
  ex50k<- extract(xmean[[i]], xy, buffer=50000, fun=mean, na.rm=TRUE)
  ex75k<- extract(xmean[[i]], xy, buffer=75000, fun=mean, na.rm=TRUE)
  ex100k<- extract(xmean[[i]], xy, buffer=100000, fun=mean, na.rm=TRUE)
  
  date<-names(xmean[[i]])
  date2 <- as.numeric(format(as.Date(date, "X%Y"), "%Y"))
  dataex<-data.frame(lon = xy[,1], lat = xy[,2], Year= date2,ex,ex10k,ex25k,ex50k,ex75k,ex100k)   
  print(date2)
  all_EU_mean_data = rbind(all_EU_mean_data, dataex)
}


colnames(all_EU_mean_data)[c(1:11)]<-c("Longitude","Latitude", "Year", "Month", "Mean_Temp","Mean_Temp1k","Mean_Temp10k","Mean_Temp25k","Mean_Temp50k","Mean_Temp75k","Mean_Temp100k")

write.csv(all_EU_mean_data, "Sensitivity_Climate.csv")

head(all_EU_mean_data)


all_EU_mean_data<-read.csv("Annual_Meant_Buff.csv")


######

doMean = function(sp_name) {
  spid2 = subset(LPI_EU4, ID == sp_name)   #subsetting the population data by each population
  spid = spid2[,64:126]                     #subsetting only the dates
  colnames(spid)<-1950:2012              #renaming the date column names as R doesn't like numbered column names
  climid=subset(all_EU_mean_data, ID == sp_name & Buffer == 0)  #subsetting the climate data by each population
  
  name<-spid2$Binomial
  id<-spid2$ID
  Date<-as.numeric(colnames(spid))
  spidt<-destring(t(spid))
  Mean_yr<-climid[climid$Year %in% Year, ]
  Mean<-Mean_yr$Mean_T
  Year<-Mean_yr$Year
  
  if (sum(is.nan(Mean))!=length(Mean)){
    
    lm_mean<-lm(Mean~Year)
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

all_df_list <- lapply(unique(LPI_clim$ID), doMean)

all_matrix <- matrix(unlist(all_df_list), ncol=6, byrow=TRUE)
mean_df <- data.frame(all_matrix)
colnames(mean_df) <- c("ID", "Term","Estimate","SE","Statistic","p.val")

#clim_sens<-mean_df[,c(1,3,4)]
#colnames(clim_sens)[c(2,3)]<-c("Estimate","SE")

sens100k<-mean_df[,c(1,3,4)]
colnames(sens100k)[c(2,3)]<-c("Estimate100k","SE100k")

clim_sens<-merge(clim_sens,sens100k,by="ID" )
head(clim_sens)

write.csv(clim_sens, "Climate_Sensitivity.csv")

