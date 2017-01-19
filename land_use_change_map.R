library(zoo)
library(raster)

crop_s<-stack(list.files(path = getwd(), pattern = "^.*crop*.*.tif$"))
gras_s<-stack(list.files(path = getwd(), pattern = "^.*gras*.*.tif$"))

crop_s<-crop_s[[2:8]]
gras_s<-gras_s[[2:8]]

extent(crop_s) == extent(gras_s)

both_1950<-crop_s[[1]]+gras_s[[1]]
both_1960<-crop_s[[2]]+gras_s[[2]]
both_1970<-crop_s[[3]]+gras_s[[3]]
both_1980<-crop_s[[4]]+gras_s[[4]]
both_1990<-crop_s[[5]]+gras_s[[5]]
both_2000<-crop_s[[6]]+gras_s[[6]]
both_2005<-crop_s[[7]]+gras_s[[7]]

both_s<-stack(both_1950,both_1960,both_1970,both_1980,both_1990,both_2000, both_2005)

both_s[both_s>1] <- 1
plot(both_s[[1]])
plot(both_s[[7]])
plot(both_s[[7]] - both_s[[1]])



names(both_s)<-c(seq(1950,2000,by=10),2005)

Year_na<-rep(NA, length(1950:2005))
Year_all<-1950:2005
Yr<-cbind(Year_all, Year_na)
colnames(Yr)<-c("Year", "Year_NA")

Year<-c(seq(1950,2000,by=10),2005)
fun_interp<-function(x){
  
  x<-matrix(x)
  year_x<-cbind(x,Year)
  colnames(year_x)<-c("LU", "Year")
  year_merge<-merge(year_x,Yr, by="Year", all=TRUE)
  LUv<-na.approx(year_merge$LU)
  LUC<-mean(diff(LUv))
  return(LUC)
}

LUC_cell<-calc(both_s, fun_interp)





