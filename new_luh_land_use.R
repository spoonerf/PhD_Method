library(raster)

yr<-as.character(850:2015)
date<-as.Date(yr, format="%Y" )

primf<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="primf")
primf<-setZ(primf, date, name="year")

primn<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="primn")
primn<-setZ(primn, date, name="year")

secdf<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="secdf")
secdf<-setZ(secdf, date, name="year")

secdn<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="secdn")
secdn<-setZ(secdn, date, name="year")

urban<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="urban")
urban<-setZ(urban, date, name="year")

pastr<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="pastr")
pastr<-setZ(pastr, date, name="year")

rnge<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="range")
rnge<-setZ(rnge, date, name="year")

c3ann<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="c3ann")
c3ann<-setZ(c3ann, date, name="year")

c4ann<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="c4ann")
c4ann<-setZ(c4ann, date, name="year")

c3per<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="c3per")
c3per<-setZ(c3per, date, name="year")

c4per<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="c4per")
c4per<-setZ(c4per, date, name="year")

c3nfx<-brick("D:/Fiona/PhD1/Land_Use/LUH2/states.nc", varname="c3nfx")
c3nfx<-setZ(c3nfx, date, name="year")


primf_cr<-primf[[1091:1166]]
primn_cr<-primn[[1091:1166]]
secdf_cr<-secdf[[1091:1166]]
secdn_cr<-secdn[[1091:1166]]
urban_cr<-urban[[1091:1166]]
pastr_cr<-pastr[[1091:1166]]
rnge_cr<-rnge[[1091:1166]]
c3ann_cr<-c3ann[[1091:1166]]
c4ann_cr<-c4ann[[1091:1166]]
c3per_cr<-c3per[[1091:1166]]
c4per_cr<-c4per[[1091:1166]]
c3nfx_cr<-c3nfx[[1091:1166]]


writeRaster(primf_cr, "primf_1940.tif")
writeRaster(primn_cr, "primn_1940.tif")
writeRaster(secdf_cr, "secdf_1940.tif")
writeRaster(secdn_cr, "secdn_1940.tif")
writeRaster(urban_cr, "urban_1940.tif")
writeRaster(pastr_cr, "pastr_1940.tif")
writeRaster(rnge_cr, "rnge_1940.tif")
writeRaster(c3ann_cr, "c3ann_1940.tif")
writeRaster(c4ann_cr, "c4ann_1940.tif")
writeRaster(c3per_cr, "c3per_1940.tif")
writeRaster(c4per_cr, "c4per_1940.tif")
writeRaster(c3nfx_cr, "cnfx_1940.tif")


primf<-brick("primf_1940.tif")
primn<-brick("primn_1940.tif")
secdf<-brick("secdf_1940.tif")
secdn<-brick("secdn_1940.tif")
urban<-brick("urban_1940.tif")
pastr<-brick("pastr_1940.tif")
rnge<-brick("rnge_1940.tif")
c3ann<-brick("c3ann_1940.tif")
c4ann<-brick("c4ann_1940.tif")
c3per<-brick("c3per_1940.tif")
c4per<-brick("c4per_1940.tif")
c3nfx<-brick("cnfx_1940.tif")

LPI<-read.csv("D:/Fiona/Git_Method/Git_Method/LPI_populations_IP_fishedit_20140310_nonconf.csv")

LPIsp<-subset(LPI, Specific_location==1 & System !="Marine" & Class != "Actinopterygii"& Class != "Cephalaspidomorphi" )

xy<-data.frame(LPIsp$Longitude, LPIsp$Latitude)

xy<-unique(xy)     #identifying unique locations to extract climate data from 

xy_df<-data.frame(xy)
colnames(xy_df)<-c("lon", "lat")
coordinates(xy_df) <- c("lon", "lat")

head(xy_df)


library(doParallel)

layer<-c3nfx

n<-6  #number of cores to use - not sure how many I can go up to
cl<-makeCluster(n)
registerDoParallel(cl)  

days<-nlayers(layer)    #splitting the data evenly between the cores
step<-floor(days/n)

ptime <- system.time({   
  df<- foreach(days, .combine=cbind) %dopar%{
    rasterex <- raster:::extract(layer[[1:days]], xy_df)
  }
}) 
ptime 
#beep(3)
stopCluster(cl)

dates<-1940:2015

datesr<-rep(dates, each=1078)

dfm<-melt(df)

lon<-xy[,1]
lat<-xy[,2]

dfm2<-cbind(lon,lat,datesr, dfm)
#dfm2<-cbind(lon,lat,dates, dfm)


#primf2<-dfm2   #do each individually
#primn2<-dfm2
#secdf2<-dfm2
#secdn2<-dfm2
#urban2<-dfm2
#pastr2<-dfm2
#rnge2<-dfm2
#c3ann2<-dfm2
#c4ann2<-dfm2
#c3per2<-dfm2
##c4per2<-dfm2
#c3nfx2<-dfm2

land_use_nat<-data.frame(primf2$lon, primf2$lat, primf2$datesr, primf2$value, primn2$value,secdf2$value, secdn2$value)

colnames(land_use_nat)<-c("Longitude", "Latitude", "Year", "Primary_Forest", "Primary_Non_Forest", 
                          "Secondary_Forest", "Secondary_Non_Forest") 

land_use_nat$total<-land_use$Primary_Forest+land_use$Primary_Non_Forest+land_use$Secondary_Forest+
  land_use$Secondary_Non_Forest

land_use<-data.frame(primf2$lon, primf2$lat, primf2$datesr, primf2$value, primn2$value, 
                     secdf2$value, secdn2$value, urban2$value, pastr2$value, rnge2$value,
                     c3ann2$value, c4ann2$value, c3per2$value, c4per2$value, c3nfx2$value)

colnames(land_use)<-c("Longitude", "Latitude","Year", "Primary_Forest", "Primary_Non_Forest", 
                      "Secondary_Forest", "Secondary_Non_Forest", "Urban", "Pasture", "Rangelands",
                     "C3_Annual", "C4_Annual", "C3_Perennial", "C4_Perennial",  "C3_Nitrogen_Fix" )

land_use$natural<-land_use$Primary_Forest+land_use$Primary_Non_Forest+land_use$Secondary_Forest+
  land_use$Secondary_Non_Forest

#land_use$total<-land_use$Primary_Forest+land_use$Primary_Non_Forest+land_use$Secondary_Forest+
 # land_use$Secondary_Non_Forest+land_use$Urban+land_use$Pasture+land_use$Rangelands+land_use$C3_Annual+
  #land_use$C4_Annual+land_use$C3_Perennial+land_use$C4_Perennial+land_use$C3_Nitrogen_Fix

LPI_ID<-LPIsp[,c("ID", "Longitude", "Latitude")]

LPILU<-merge(LPI_ID,land_use, by=c("Longitude", "Latitude"))


write.csv(LPILU, "LUH2_Land_Use_All_LPI.csv")


library(taRifx)

doDist = function(sp_name) {
  spid2 = subset(LPIsp, ID == sp_name)   #subsetting the population data by each population
  spid = spid2[,63:118]                     #subsetting only the dates
  colnames(spid)<-1950:2005              #renaming the date column names as R doesn't like numbered column names
  lu_id=subset(LPILU, ID == sp_name)  #subsetting the climate data by each population
  
  id<-spid2$ID
  Date<-as.numeric(colnames(spid))
  spidt<-destring(t(spid))
  
  if (sum(!is.na(spidt)) > 1) {
    Yr<-Date[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
    luyr<-as.matrix(subset(lu_id, Year>=min(Yr) & Year<= max(Yr))[,c(4,17)])   #6:11 for all categories 12:13 for nat vs anth
    luyr<-luyr[order(luyr[,1]),]
    natdiff<-mean(diff(luyr[,2]))
    
    # lum<-matrix(luyr, ncol=ncol(luyr))
    # m<-dist(lum, method="euclidean")
    # n<-as.matrix(m)
    # lu<-diag(n[c(2:nrow(n)),c(1:ncol(n))])
    # euc<-mean(lu)
    
  } else{
   # euc<-NA
    natdiff<-NA
  }
  #euc_df<-data.frame(id,euc)
 # euc_df<-data.frame(id,euc,natdiff)
  euc_df<-data.frame(id, natdiff)
  print(euc_df)  
  return(euc_df)
}

all_df_list <- lapply(unique(LPIsp$ID), doDist)

all_matrix <- matrix(unlist(all_df_list), ncol=2, byrow=TRUE)
mean_df <- data.frame(all_matrix)
colnames(mean_df) <- c("ID", "Nat_change")

write.csv(mean_df, "LUC_average_annual_change_nat_025.csv")
