#install.packages("devtools")
library(devtools)
install_github("timnewbold/GISOperations")
library(GISOperations)

#cellareas <- DegreeCellAreaKM(lat=coordinates(crop1940)[,2],height=res(crop1940)[2],width=res(crop1940)[1])


#crop_pcnt_1940<-values(crop1940)/cellareas

library(raster)
crop1940<-raster("crop1940AD.asc")
gras1940<-raster("gras1940AD.asc")

crop1950<-raster("crop1950AD.asc")
gras1950<-raster("gras1950AD.asc")

crop1960<-raster("crop1960AD.asc")
gras1960<-raster("gras1960AD.asc")

crop1970<-raster("crop1970AD.asc")
gras1970<-raster("gras1970AD.asc")

crop1980<-raster("crop1980AD.asc")
gras1980<-raster("gras1980AD.asc")

crop1990<-raster("crop1990AD.asc")
gras1990<-raster("gras1990AD.asc")

crop2000<-raster("crop2000AD.asc")
gras2000<-raster("gras2000AD.asc")

crop2005<-raster("crop2005AD.asc")
gras2005<-raster("gras2005AD.asc")


rlist<-list(crop1940,crop1950,crop1960,crop1970,crop1980,crop1990,crop2000,crop2005,gras1940,gras1950,gras1960,gras1970,gras1980,gras1990,gras2000,gras2005)
cellareas <- DegreeCellAreaKM(lat=coordinates(crop1940)[,2],height=res(crop1940)[2],width=res(crop1940)[1])



for (map in rlist){
  
  map_pcnt<-signif(values(map),digits=4)/cellareas
  map_pcnt_m<-matrix(map_pcnt, nrow=nrow(map), ncol=ncol(map), byrow=TRUE)
  map_pcnt_r<-raster(map_pcnt_m, xmn=map@extent[1], xmx=map@extent[2], ymn=map@extent[3], ymx=map@extent[4])
  writeRaster(map_pcnt_r, paste(names(map),"raster.tif", sep="_"),overwrite=TRUE)
  print(names(map))
}


crop_s<-stack(list.files(path = getwd(), pattern = "^.*crop*.*.tif$"))
gras_s<-stack(list.files(path = getwd(), pattern = "^.*gras*.*.tif$"))


nat_1940<-(1 - (crop_s[[1]] + gras_s[[1]]))
nat_2005<-(1 - (crop_s[[7]] + gras_s[[7]]))
plot(nat_2005 - nat_1940)

###########LPI data

LPI_LUC<-read.csv("LPI_pops_20160523_edited.csv")

LPI_pop<-subset(LPI_LUC, Specific_location==1 & System !="Marine" & Class != "Actinopterygii"& Class != "Cephalaspidomorphi")

ID<-LPI_pop$ID
pop_data<- LPI_pop[,c(1,65:130)]

pop_datab <- (pop_data [,2:67] !="NULL")
points_per_pop1950_2012 = rowSums(pop_datab)
length_id <- data.frame(ID,points_per_pop1950_2012)

LPI_LUC<-merge(length_id, LPI_pop, by = "ID")
LPI_LUC<-subset(LPI_pop, points_per_pop1950_2012 >=2)


id<-LPI_LUC$ID
latlong<-cbind(LPI_LUC$Longitude,LPI_LUC$Latitude)

x<-latlong[,1]
y<-latlong[,2]

xy<-cbind(id,x,y) #creating the grid around the population - 5km either side gives a grid of 121km^2
xy<-na.omit(xy)
id<-xy[,1]

points<-SpatialPoints(cbind(xy[,2], xy[,3]))


xmin<- colFromX(crop_s[[1]], points[1:length(points),]) - 1
xmax<- colFromX(crop_s[[1]], points[1:length(points),]) + 1
ymin<- rowFromY(crop_s[[1]], points[1:length(points),]) - 1
ymax<- rowFromY(crop_s[[1]], points[1:length(points),]) + 1    #changed from 5 to 12


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
    both_sum<-sum(diff(subset(Yr_intrp_crop, Year >= min_yr & Year <= max_yr)$both))
  }else{
    crop_change<-NA
    gras_change<-NA
    both_change<-NA
    both_sum<-NA
  }
  
  years<-max(Year)-min(Year)
  final<-cbind(ID,Binomial,crop_change, gras_change, both_change, both_sum, years) 
  result<-rbind(final,result)
  print(final)
  
}

#write.csv(result, "Hyde_crop_pasture_annual_change.csv")
write.csv(result, "Hyde_crop_pasture_annual_change_sum.csv")


head(result)

hist(result$both_change)

hist(as.numeric(as.character(result$crop_change)))


