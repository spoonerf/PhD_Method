library(raster)

CR40s<-brick("cru_ts3.23.1941.1950.tmp.dat.nc")
CR50s<-brick("cru_ts3.23.1951.1960.tmp.dat.nc")
CR60s<-brick("cru_ts3.23.1961.1970.tmp.dat.nc")
CR70s<-brick("cru_ts3.23.1971.1980.tmp.dat.nc")
CR80s<-brick("cru_ts3.23.1981.1990.tmp.dat.nc")
CR90s<-brick("cru_ts3.23.1991.2000.tmp.dat.nc")
CR00s<-brick("cru_ts3.23.2001.2010.tmp.dat.nc")

CR<-stack(CR40s[[109:120]],CR50s,CR60s,CR70s,CR80s,CR90s,CR00s[[c(1:60)]])

Jan<-seq(1, nlayers(CR), 12)
Dec<-Jan+11

year_s<-stack()

for (i in 1:length(Jan)){
  
  year<-mean(CR[[Jan[i]:Dec[i]]])
  year_s<-stack(year, year_s)
   
}

years<-1950:2005
names(year_s)<-years

##plotting each year
# for (i in 1:nlayers(year_s)){
#   
#   plot(year_s[[i]], main=i+1949)
# }
# 

year_test<-year_s[10000:11000]

year_row<-matrix()

for (i in 1:ncell(year_test)){
  
  cell<-as.vector(year_s[i])
  
  if ((sum(is.na(cell))) != length(cell)){
    
  slope<-summary(lm(cell~years))$coefficients[2]
  
  } else {
    
    slope<-NA
  }
  
  year_row<-rbind(slope, year_row)
  
}

cell_val<-numeric(ncell(year_s))

for (i in 1:ncell(year_s)){
  
  cell_na<-!is.na(year_s[[1]][i])
  cell_val<-rbind(cell_na,cell_val)
  
}








