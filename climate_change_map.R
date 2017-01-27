CR40s<-brick("cru_ts3.23.1941.1950.tmp.dat.nc")
CR50s<-brick("cru_ts3.23.1951.1960.tmp.dat.nc")
CR60s<-brick("cru_ts3.23.1961.1970.tmp.dat.nc")
CR70s<-brick("cru_ts3.23.1971.1980.tmp.dat.nc")
CR80s<-brick("cru_ts3.23.1981.1990.tmp.dat.nc")
CR90s<-brick("cru_ts3.23.1991.2000.tmp.dat.nc")
CR00s<-brick("cru_ts3.23.2001.2010.tmp.dat.nc")

CR<-stack(CR40s[[109:120]],CR50s,CR60s,CR70s,CR80s,CR90s,CR00s[[c(1:60)]])

year_mon<-rep(1:(nlayers(CR)/12), each=12)
CR_year<-stackApply(CR, indices=year_mon, fun=mean, na.rm=TRUE)

time<-1:nlayers(CR_year)
X <- cbind(1, time)

## pre-computing constant part of least squares
invXtX <- solve(t(X) %*% X) %*% t(X)

## much reduced regression model; [2] is to get the slope
quickfun <- function(y) (invXtX %*% y)[2]
x4 <- calc(CR_year, quickfun) 

writeRaster(x4, "Global_Rate_Mean_Temp_Change.tif", overwrite=TRUE)
#
cellStats(x4, mean)

plot(x4)
















