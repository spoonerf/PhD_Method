library(zoo)
library(raster)

crop_s<-stack(list.files(path = getwd(), pattern = "^.*crop*.*.tif$"))
gras_s<-stack(list.files(path = getwd(), pattern = "^.*gras*.*.tif$"))

crop_s<-crop_s[[2:8]]
gras_s<-gras_s[[2:8]]


both_s<-crop_s+gras_s
plot(both_s[[1]])

both_s[both_s>1] <- 1
plot(both_s[[1]])

names(both_s)<-c(seq(1950,2000,by=10),2005)

Year_na<-rep(NA, length(1950:2005))
Year_all<-1950:2005
Yr<-cbind(Year_all, Year_na)
colnames(Yr)<-c("Year", "Year_NA")

Yr_intrp<-merge(Yr, mean_crop_df2, by="Year", all=T)

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





