install.packages("ncdf4")
library(raster)
library(ncdf4)

rmean<-brick("tg_0.25deg_reg_v11.0.nc", varname = "tg")

tm <- seq(as.Date('1950-01-01'), as.Date('2014-12-31'), 'day')
s <- setZ(rmean, tm, 'days')

tmy <- seq(as.Date('1950-01-01'), as.Date('2014-12-31'), 'year')


rmn<-zApply(s, tmy, fun=mean, name='year')

rasterex <- raster:::extract(rmn, xy_df, fun=mean, na.rm=TRUE)
raster10 <- raster:::extract(rmn, xy_df, buffer=10000, fun=mean, na.rm=TRUE)
raster25 <- raster:::extract(rmn, xy_df, buffer=25000, fun=mean, na.rm=TRUE)
raster50 <- raster:::extract(rmn, xy_df, buffer=50000, fun=mean, na.rm=TRUE)
raster100 <- raster:::extract(rmn, xy_df, buffer=100000, fun=mean, na.rm=TRUE)



n<-8  #number of cores to use - not sure how many I can go up to
cl<-makeCluster(n)
registerDoParallel(cl) 

yrs<-nlayers(rmn)
step<-floor(yrs/n)

ptime <- system.time({df<- foreach(lyr=seq(1,yrs,step)[1:n], .combine=cbind) %dopar%{
  
}
}) 
ptime



stopCluster(cl)










