library(devtools)
install_github("timnewbold/GISOperations")
library(GISOperations)

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

crop_na<-crop_s[[1]]
values(crop_na)<-NA

crop_s<-crop_s[[2:8]]
gras_s<-gras_s[[2:8]]

agri_s<-crop_s+gras_s

s <- stack(replicate(length(1950:2005), crop_na))

hyde_years<-c(1950,1960,1970,1980,1990,2000,2005)
years<-1950:2005

hyde_layers<-which(years %in% hyde_years)

for (i in 1:length(hyde_layers)){
  s[[hyde_layers[i]]]<-agri_s[[i]]
  print(i)
}


s[[1]][is.na(s[[1]])]<-0
s[[56]][is.na(s[[56]])]<-0


hyde_interp<-calc(s, fun = na.approx)

names(hyde_interp)<-paste("Year_", 1950:2005, sep="")

writeRaster(hyde_interp, "Hyde_Interpolated.tif")









