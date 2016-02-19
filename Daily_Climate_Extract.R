install.packages("foreach")
install.packages("doSNOW")
library(foreach)
library(doSNOW)
library(raster)

rmax<-brick("tx_0.25deg_reg_v11.0.nc", varname = "tx")
rmin<-brick("tn_0.25deg_reg_v11.0.nc", varname = "tn")
rpcp<-brick("rr_0.25deg_reg_v11.0.nc", varname = "rr")


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

test<-extract(rmax[[23741]], xy, buffer=50000, fun=mean, na.rm=TRUE)
df<-data.frame(id, test)

df2<-subset(df, !is.na(test))

LPI_EU5<-LPI_EU4[LPI_EU4$ID %in% df2$id,]

id<-LPI_EU5$ID
xy<-cbind(LPI_EU5$Longitude,LPI_EU5$Latitude)
binomial<-LPI_EU5$Binomial

all_EU_max_data<-data.frame(id=numeric(0), Binomial=character(0),xy=numeric(0), year=numeric(0), month=numeric(0), day=numeric(0), rasterex=numeric(0))

for (i in 1:nlayers(rmax)){
  
  rasterex <- raster:::extract(rmax[[i]], xy, buffer=50000, fun=mean, na.rm=TRUE)
  date<-names(rmax[[i]])
  dataex<-data.frame(id, binomial,xy, substr(date, 2,5),substr(date, 7,8), substr(date, 10,11), rasterex)   
  print(date)
  all_EU_max_data = rbind(all_EU_max_data, dataex)
  
}

write.csv(all_EU_max_data, "all_EU_daily_max_temp.csv")


id<-Inc5MammalsTemp$ID
xy<-cbind(Inc5MammalsTemp$Longitude,Inc5MammalsTemp$Latitude)

binomial<-Inc5MammalsTemp$Binomial





