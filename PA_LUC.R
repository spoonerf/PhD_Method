PA<-read.csv("PA_ID.csv")

crop_s<-stack(list.files(path = getwd(), pattern = "^.*crop*.*.tif$"))
gras_s<-stack(list.files(path = getwd(), pattern = "^.*gras*.*.tif$"))

LPI_LUC<-read.csv("LPI_pops_20160523_edited.csv")

LPI_PA<-merge(LPI_LUC, PA, by="ID")

head(LPI_PA)
names(LPI_PA)[names(LPI_PA) == 'Protected_status.y'] <- 'PA_status'

ID<-LPI_PA$ID
pop_data<- LPI_pop[,c(1,65:130)]

pop_datab <- (pop_data [,2:67] !="NULL")
points_per_pop1950_2012 = rowSums(pop_datab)
length_id <- data.frame(ID,points_per_pop1950_2012)

LPI_LUC<-merge(length_id, LPI_PA, by = "ID")
LPI_LUC<-subset(LPI_PA, points_per_pop1950_2012 >=2)


id<-LPI_LUC$ID
latlong<-cbind(LPI_LUC$Longitude,LPI_LUC$Latitude)

x<-latlong[,1]
y<-latlong[,2]

xy<-cbind(id,x,y) 
xy<-na.omit(xy)
id<-xy[,1]

points<-SpatialPoints(cbind(xy[,2], xy[,3]))


xmin<- colFromX(crop_s[[1]], points[1:length(points),]) - 1
xmax<- colFromX(crop_s[[1]], points[1:length(points),]) + 1
ymin<- rowFromY(crop_s[[1]], points[1:length(points),]) - 1
ymax<- rowFromY(crop_s[[1]], points[1:length(points),]) + 1 


grid_crop<-data.frame(ID=id,xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)   #setting extents of grid to extract

grid_crop2<-merge(grid_crop,LPI_LUC, by="ID")

grid_crop2<-subset(grid_crop2, Specific_location ==1) 

nrow(grid_crop2)

library(taRifx)
library(zoo)

result <- data.frame() #empty result dataframe

for (i in 1:length(grid_crop2$ID)){
  
  ID<-grid_crop2[i,1]
  Binomial<-as.character(grid_crop2[i,6])
  spid = grid_crop2[i,69:134]                     #subsetting only the dates
  colnames(spid)<-1950:2015
  
  Date<-as.numeric(colnames(spid))
  spidt<-destring(t(spid))
  
  Year<-Date[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  Population<-spidt[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  
  crop_check<-crop(crop_s, extent(crop_s, grid_crop2[i,4],grid_crop2[i,5],grid_crop2[i,2],grid_crop2[i,3]))
  gras_check<-crop(gras_s, extent(gras_s, grid_crop2[i,4],grid_crop2[i,5],grid_crop2[i,2],grid_crop2[i,3]))
  
  crop_df<-data.frame(as.matrix(crop_check))
  gras_df<-data.frame(as.matrix(gras_check))
  
  decs<-c("1940","1950","1960","1970","1980","1990","2000","2005")
  colnames(crop_df)<-decs
  colnames(gras_df)<-decs
  
  mean_crop_df<-colMeans(crop_df)
  mean_crop_df2<-data.frame(mean_crop_df, decs)
  colnames(mean_crop_df2)<-c("mean_crop", "Year")
  
  mean_gras_df<-colMeans(gras_df)
  mean_gras_df2<-data.frame(mean_gras_df, decs)
  colnames(mean_gras_df2)<-c("mean_grass", "Year")
  
  Year_na<-rep(NA, length(1940:2005))
  Year_all<-1940:2005
  Yr<-cbind(Year_all, Year_na)
  colnames(Yr)<-c("Year", "Year_NA")
  
  Yr_intrp_crop<-merge(Yr, mean_crop_df2, by="Year", all=T)
  Yr_intrp_gras<-merge(Yr, mean_gras_df2, by="Year", all=T)
  
  min_yr<-min(Year)
  max_yr<-max(Year)
  
  
  if (min_yr != max_yr & min_yr < 2005 & sum(is.na(Yr_intrp_crop$mean_crop))<60) {
    
    Yr_intrp_crop$mean_crop_int<-na.approx(Yr_intrp_crop$mean_crop)
    Yr_intrp_gras$mean_gras_int<-na.approx(Yr_intrp_gras$mean_grass)
    Yr_intrp_crop$both<-Yr_intrp_crop$mean_crop_int + Yr_intrp_gras$mean_gras_int
    
    crop_change<-mean(diff(subset(Yr_intrp_crop, Year >= min_yr & Year <= max_yr)$mean_crop_int))
    gras_change<-mean(diff(subset(Yr_intrp_gras, Year >= min_yr & Year <= max_yr)$mean_gras_int))
    both_change<-mean(diff(subset(Yr_intrp_crop, Year >= min_yr & Year <= max_yr)$both))
    
  }else{
    crop_change<-NA
    gras_change<-NA
    both_change<-NA
  }
  
  final<-cbind(ID,Binomial,crop_change, gras_change, both_change) 
  result<-rbind(final,result)
  print(final)
  
}
