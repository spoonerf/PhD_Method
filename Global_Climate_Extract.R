library(raster)
HC<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/HadCRUT.4.4.0.0.median.nc")
CT<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/CRUTEM.4.4.0.0.anomalies.nc")

CR40s<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/cru_ts3.23.1941.1950.tmp.dat.nc/cru_ts3.23.1941.1950.tmp.dat.nc")
CR50s<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/cru_ts3.23.1951.1960.tmp.dat.nc/cru_ts3.23.1951.1960.tmp.dat.nc")
CR60s<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/cru_ts3.23.1961.1970.tmp.dat.nc/cru_ts3.23.1961.1970.tmp.dat.nc")
CR70s<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/cru_ts3.23.1971.1980.tmp.dat.nc/cru_ts3.23.1971.1980.tmp.dat.nc")
CR80s<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/cru_ts3.23.1981.1990.tmp.dat.nc/cru_ts3.23.1981.1990.tmp.dat.nc")
CR90s<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/cru_ts3.23.1991.2000.tmp.dat.nc/cru_ts3.23.1991.2000.tmp.dat.nc")
CR00s<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/cru_ts3.23.2001.2010.tmp.dat.nc/cru_ts3.23.2001.2010.tmp.dat.nc")
CR10s<-brick("C:/Users/Fiona/Desktop/PhD/Climate/Global_HadCRUT/cru_ts3.23.2011.2014.tmp.dat.nc/cru_ts3.23.2011.2014.tmp.dat.nc")

plot(CR40s[[11]])

LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")

LPIsp<-subset(LPI, Specific_location==1 & System !="Marine" & Class != "Actinopterygii"& Class != "Cephalaspidomorphi" )

CR<-stack(CR40s,CR50s,CR60s,CR70s,CR80s,CR90s,CR00s,CR10s)

plot(CR[[1]])
points(LPIsp$Longitude, LPIsp$Latitude)

xy<-data.frame(LPIsp$Longitude, LPIsp$Latitude)

test<-extract(CR[[1]], xy, buffer=5000, na.rm=T)


unlist(test)





