body<-read.csv("bird_and_mammal_traits2.csv")
luc<-read.csv("LUC_average_annual_change_nat_025.csv")
LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")
Realm<-read.csv("selected_pops_Ecoregion.csv")
Realm<-Realm[,c("ID", "WWF_REALM2")]

pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_18.csv")
EurHil<-read.csv("Europe_HILDA_5_year_pops.csv")  # data from Euro-centric analysis

temp<-temp[,c("ID", "Estimate")]

LPI<-LPI[,c("ID","Binomial","Common_name", "Order", "Protected_status", "Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat", "Migratory")]

df<-merge(merge(temp,luc, by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

df<-merge(df, Realm, by="ID")

df<-merge(df, body[,c(3:5)], by="ID", all=TRUE)     #41 pops bodysizes missing for birds

df2<-subset(df, !is.na(Estimate)&r_sq >= 0.5  & !is.na(Nat_change)&length_time >=10 & System!="Marine" &Specific_location == 1 & (Class=="Mammalia") & !is.na(Bodymass) )

nrow(df2)

df2$Nat_loss<-0-df2$Nat_change 

df2[is.na(df2$lambda_mean),]$lambda_mean<-0

pyr<-subset(df2, Common_name =="Pyrenean chamois")    #record 11470 had wrong longitude - in Russia!

plot(pyr$Longitude, pyr$Latitude)

pyrs<-pyr[,c("ID","Longitude","Latitude","lambda_mean")]

id<-pyrs$ID
lam<-as.numeric(pyrs$lambda_mean)

library(rgdal)

pyrxy<-SpatialPoints(pyr[,c("Longitude","Latitude")])

library(raster)

e<-extent(pyrxy)

r<-raster(e, resolution=0.25)

rz<-rasterize(pyrxy,r,lam )
rid<-rasterize(pyrxy,r,id)
plot(rz)
plot(rid)

rz_spdf<-xyFromCell(rz, 1:108)

rzm<-as.vector(rz)
ridm<-as.vector(rid)
patchID<-1:length(rid)

df<-cbind(patchID,rz_spdf,rzm,ridm)
colnames(df)[c(4,5)]<-c("area", "ID")

write.csv(df, "Rupicapra_pyrenaica_populations.csv")


lnd<-read.csv("LUH2_Land_Use_All_LPI.csv")

lnd_rp<-merge(df,lnd, by="ID")

lnd_rp2<-unique(lnd_rp[,c("ID","patchID", "x", "y")])

library(reshape2)

lnd_rp2<-lnd_rp[,c("ID","Year","natural")]

lnd_rp_cast<-dcast(lnd_rp2, ID ~ Year )

####geographic data

library(raster)

yr<-as.character(850:2015)
date<-as.Date(yr, format="%Y" )

primf<-brick("states.nc", varname="primf")
primf<-setZ(primf, date, name="year")

primn<-brick("states.nc", varname="primn")
primn<-setZ(primn, date, name="year")

secdf<-brick("states.nc", varname="secdf")
secdf<-setZ(secdf, date, name="year")

secdn<-brick("states.nc", varname="secdn")
secdn<-setZ(secdn, date, name="year")

primf_cr<-primf[[1091:1166]]
primn_cr<-primn[[1091:1166]]
secdf_cr<-secdf[[1091:1166]]
secdn_cr<-secdn[[1091:1166]]


# writeRaster(primf_cr, "primf_1940.tif")
# writeRaster(primn_cr, "primn_1940.tif")
# writeRaster(secdf_cr, "secdf_1940.tif")
# writeRaster(secdn_cr, "secdn_1940.tif")

primf<-brick("primf_1940.tif")
primn<-brick("primn_1940.tif")
secdf<-brick("secdf_1940.tif")
secdn<-brick("secdn_1940.tif")

n40<-primf[[1]]+ primn[[1]]+secdf[[1]]+secdn[[1]]
n40c<-crop(n40, e)

n50<-primf[[11]]+ primn[[11]]+secdf[[11]]+secdn[[11]]
n50c<-crop(n50,e)

n60<-primf[[21]]+ primn[[21]]+secdf[[21]]+secdn[[21]]
n60c<-crop(n60,e)

n70<-primf[[31]]+ primn[[31]]+secdf[[31]]+secdn[[31]]
n70c<-crop(n70,e)

n80<-primf[[41]]+ primn[[41]]+secdf[[41]]+secdn[[41]]
n80c<-crop(n80,e)

n90<-primf[[51]]+ primn[[51]]+secdf[[51]]+secdn[[51]]
n90c<-crop(n90,e)

n2000<-primf[[61]]+ primn[[61]]+secdf[[61]]+secdn[[61]]
n2000c<-crop(n2000,e)

n2010<-primf[[71]]+ primn[[71]]+secdf[[71]]+secdn[[71]]
n2010c<-crop(n2010,e)


n40m<-as.vector(n40c)
n50m<-as.vector(n50c)
n60m<-as.vector(n60c)
n70m<-as.vector(n70c)
n80m<-as.vector(n80c)
n90m<-as.vector(n90c)
n2000m<-as.vector(n2000c)
n2010m<-as.vector(n2010c)

nall<-cbind(n40m,n50m,n60m,n70m,n80m,n90m,n2000m,n2010m)
colnames(nall)<-c("period_1940", "period_1950", "period_1960", "period_1970", "period_1980", "period_1990", "period_2000", "period_2010")

land_use_map<-cbind(df[,c(1:3)], nall) #percentage cover of natural land use (primary+secondary) for each decade 1940-2010

no_yrs_mine<-10

###### demographic information

library(popbio)

















