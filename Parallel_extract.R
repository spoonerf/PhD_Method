library(doParallel)

##adding in large raster brick - about 4GB
rmax<-brick("D:/Fiona/Git_Method/Git_Method/tx_0.25deg_reg_v11.0.nc", varname = "tx")
rmin<-brick("D:/Fiona/Git_Method/Git_Method/tn_0.25deg_reg_v11.0.nc", varname = "tn")
rmean<-brick("D:/Fiona/Git_Method/Git_Method/tg_0.25deg_reg_v11.0.nc", varname = "tg")
rpcp<-brick("D:/Fiona/Git_Method/Git_Method/rr_0.25deg_reg_v11.0.nc", varname = "rr")


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
id<-LPI_EU4$ID
test<-extract(rmin[[23741]], xy, buffer=50000, fun=mean, na.rm=TRUE)
df<-data.frame(id, test)
df2<-subset(df, !is.na(test))

LPI_EU5<-LPI_EU4[LPI_EU4$ID %in% df2$id,]

xy<-cbind(LPI_EU5$Longitude,LPI_EU5$Latitude)
loc<-data.frame(LPI_EU5$ID, LPI_EU5$Longitude, LPI_EU5$Latitude)
colnames(loc)<-c("ID", "Longitude", "Latitude")

xy<-unique(xy)     #identifying unique locations to extract climate data from 

xy_df<-data.frame(xy)
colnames(xy_df)<-c("lon", "lat")
coordinates(xy_df) <- c("lon", "lat")

length(xy_df)
###parallellise

n<-6  #number of cores to use - not sure how many I can go up to
cl<-makeCluster(n)
registerDoParallel(cl)  

days<-nlayers(rpcp)    #splitting the data evenly between the cores
step<-floor(days/n)

ptime <- system.time({   
  df<- foreach(lyr=seq(1,days,step)[1:6], .combine=cbind) %dopar%{
    rasterex <- raster:::extract(rpcp[[lyr:(lyr+step-1)]], xy_df, buffer=50000, fun=mean, na.rm=TRUE)
  }
}) 
ptime  

stopCluster(cl)

df2<-data.frame(xy_df, df)
colnames(df2)[1:2]<-c("Longitude", "Latitude")

df3<-merge(loc, df2, by=c("Longitude", "Latitude"))

crex<-subset(df3, ID=327)

write.csv(df3, "Daily_Prec_All_EU.csv")

#########################
###formatting the data###
#########################
df3<-read.csv("Daily_Mean_Temp_All_EU.csv")

pcps<-stack(df3)


long<-pcps[c(1:1047),1]
lat<-pcps[c(1048:2094),1]
id<-pcps[c(2095:3141),1]
   #check theseeee

vals<-pcps[c(3142:24854733),]

#23736 days

longr<-rep(long, 23736)
latr<-rep(lat, 23736)
idr<-rep(id, 23736)

pcp_df<-data.frame(idr,latr,longr,vals)
pcp_df$ind<-as.Date(pcp_df$ind, "X%Y.%m.%d")

head(pcp_df)

write.csv(pcp_df, "Daily_Prec_All_EU_form.csv")




df3<-read.csv("Daily_Mean_Temp_All_EU_form.csv")



