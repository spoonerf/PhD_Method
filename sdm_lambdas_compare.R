
wd<-getwd()

lf<-list.files(paste(wd, "Legion/snow_cervus_bias_new/", sep=""))
files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]

library(raster)
cervus<-stack(paste(wd, "/Legion/snow_cervus_bias_new/", files, sep = ""))
capra<-stack(paste( wd, "/Legion/snow_capra_bias/", files, sep = ""))
ursus<-stack(paste(wd,"/Legion/snow_bear_bias_new/", files, sep = ""))
gulo_n<-stack(paste(wd, "/Gulo_gulo_bias_NAm/SDM_folder_bias/", files, sep =""))
gulo_e<-stack(paste(wd, "/Gulo_gulo_bias_EU/SDM_folder_bias_EU/", files, sep =""))
lepus<-stack(paste(wd, "/Lepus_americanus_bias/SDM_folder_bias/", files, sep=""))
wtd<-stack(paste(wd, "/Odocoileus_virginianus_bias/SDM_folder_bias/",files, sep = ""))
rangifer_e<-stack(paste(wd, "/Rangifer_tarandus_bias/SDM_folder_bias/", files, sep = ""))
rangifer_n<-stack(paste(wd, "/Rangifer_tarandus_bias_NAm/SDM_folder_bias_NAm/", files, sep = ""))
polar_n<-stack(paste(wd,"/Ursus_maritimus_bias/SDM_folder_bias/", files, sep = ""))
harte<-stack(paste(wd, "/Alcelaphus_buselaphus_bias/SDM_folder_bias/", files, sep=""))
phaco<-stack(paste(wd, "/Phacochoerus_africanus_bias/SDM_folder_bias/", files, sep=""))
giraffe<-stack(paste(wd,"/Giraffa_camelopardalis_bias/SDM_folder_bias/", files, sep=""))  
kobus<-stack(paste(wd, "/Kobus_ellipsiprymnus_bias/SDM_folder_bias/", files, sep=""))
wilde<-stack(paste(wd,"/Connochaetes_taurinus_bias/SDM_folder_bias/", files, sep=""))
py_chamois<-stack(paste(wd,"/Rupicapra_pyrenaica_bias/SDM_folder_bias/", files, sep=""))
roe_deer<-stack(paste(wd,"/Capreolus_capreolus_bias/SDM_folder_bias/", files, sep=""))
equus<-stack(paste(wd,"/Equus_quagga_bias/SDM_folder_bias/", files, sep=""))


lpi<-read.csv("LPI_pops_20160523_edited.csv")

cervus_lpi<-lpi[lpi$Binomial == "Cervus_elaphus" & lpi$Specific_location == 1 & lpi$Region == "Europe", ]
capra_lpi<-lpi[lpi$Binomial == "Capra_ibex" & lpi$Specific_location == 1 & lpi$Region == "Europe", ]
ursus_lpi<-lpi[lpi$Binomial == "Ursus_arctos" & lpi$Specific_location == 1 & lpi$Region == "North America", ]
gulo_n_lpi<-lpi[lpi$Binomial == "Gulo_gulo" & lpi$Specific_location == 1 & lpi$Region == "North America", ]
gulo_e_lpi<-lpi[lpi$Binomial == "Gulo_gulo" & lpi$Specific_location == 1 & lpi$Region != "North America", ]
lepus_lpi<-lpi[lpi$Binomial == "Lepus_americanus" & lpi$Specific_location == 1, ]
wtd_lpi<-lpi[lpi$Binomial == "Odocoileus_virginianus" & lpi$Specific_location == 1 & lpi$Region == "North America", ]
rang_e_lpi<-lpi[lpi$Binomial == "Rangifer_tarandus" & lpi$Specific_location == 1 & lpi$Region != "North America", ]
rang_n_lpi<-lpi[lpi$Binomial == "Rangifer_tarandus" & lpi$Specific_location == 1 & lpi$Region == "North America", ]
polar_n_lpi<-lpi[lpi$Binomial == "Ursus_maritimus" & lpi$Specific_location == 1,]
harte_lpi<-lpi[lpi$Binomial == "Alcelaphus_buselaphus" & lpi$Specific_location == 1,]
phaco_lpi<-lpi[lpi$Binomial == "Phacochoerus_africanus" & lpi$Specific_location == 1,]
giraffe_lpi<-lpi[lpi$Binomial == "Giraffa_camelopardalis" & lpi$Specific_location == 1,]
kobus_lpi<-lpi[lpi$Binomial == "Kobus_ellipsiprymnus" & lpi$Specific_location == 1, ]
wilde_lpi<-lpi[lpi$Binomial == "Connochaetes_taurinus" & lpi$Specific_location == 1, ]
py_chamois_lpi<-lpi[lpi$Binomial == "Rupicapra_pyrenaica" & lpi$Specific_location == 1, ]
roe_deer_lpi<-lpi[lpi$Binomial == "Capreolus_capreolus" & lpi$Specific_location == 1, ]
equus_lpi<-lpi[lpi$Binomial == "Equus_burchellii" & lpi$Specific_location == 1, ]

sdm_cervus<-raster:::extract(cervus, cbind(cervus_lpi$Longitude, cervus_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_capra<-raster:::extract(capra, cbind(capra_lpi$Longitude, capra_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_ursus<-raster:::extract(ursus, cbind(ursus_lpi$Longitude, ursus_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_gulon<-raster:::extract(gulo_n, cbind(gulo_n_lpi$Longitude, gulo_n_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_guloe<-raster:::extract(gulo_e, cbind(gulo_e_lpi$Longitude, gulo_e_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_lepus<-raster:::extract(lepus, cbind(lepus_lpi$Longitude, lepus_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_wtd<-raster:::extract(wtd, cbind(wtd_lpi$Longitude, wtd_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_range<-raster:::extract(rangifer_e, cbind(rang_e_lpi$Longitude, rang_e_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_rangn<-raster:::extract(rangifer_n, cbind(rang_n_lpi$Longitude, rang_n_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_polar<-raster:::extract(polar_n, cbind(polar_n_lpi$Longitude, polar_n_lpi$Latitude), buffer= 50000, fun  = mean, na.rm = TRUE)
sdm_phaco<-raster:::extract(phaco, cbind(phaco_lpi$Longitude, phaco_lpi$Latitude), buffer= 50000, fun = mean, na.rm = TRUE)
sdm_harte<-raster:::extract(harte, cbind(harte_lpi$Longitude, harte_lpi$Latitude), buffer = 50000)#, fun = mean, na.rm = TRUE)
sdm_giraffe<-raster:::extract(giraffe, cbind(giraffe_lpi$Longitude, giraffe_lpi$Latitude), buffer= 50000)#, fun = mean, na.rm = TRUE)
sdm_kobus<-raster:::extract(kobus, cbind(kobus_lpi$Longitude, kobus_lpi$Latitude), buffer= 50000)#, fun = mean, na.rm = TRUE)
sdm_wilde<-raster:::extract(wilde, cbind(wilde_lpi$Longitude, wilde_lpi$Latitude), buffer= 50000, fun = mean, na.rm =TRUE)
sdm_py_chamois<-raster:::extract(py_chamois, cbind(py_chamois_lpi$Longitude, py_chamois_lpi$Latitude), buffer= 50000, fun = mean, na.rm =TRUE)
sdm_roe_deer<-raster:::extract(roe_deer, cbind(roe_deer_lpi$Longitude, roe_deer_lpi$Latitude), buffer= 50000, fun = mean, na.rm =TRUE)
sdm_equus<-raster:::extract(equus, cbind(equus_lpi$Longitude, equus_lpi$Latitude), buffer= 50000)#, fun = mean, na.rm =TRUE)


sdm_cervus2<-data.frame(cervus_lpi$ID,as.character(cervus_lpi$Binomial), sdm_cervus)
colnames(sdm_cervus2)<-c( "ID","Binomial", 1950:2005)
sdm_capra2<-data.frame(capra_lpi$ID,as.character(capra_lpi$Binomial), sdm_capra)
colnames(sdm_capra2)<-c( "ID","Binomial", 1950:2005)
sdm_ursus2<-data.frame(ursus_lpi$ID,as.character(ursus_lpi$Binomial) ,sdm_ursus)
colnames(sdm_ursus2)<-c( "ID","Binomial", 1950:2005)
sdm_gulon2<-data.frame(gulo_n_lpi$ID,as.character(gulo_n_lpi$Binomial) ,sdm_gulon)
names(sdm_gulon2)<-c( "ID","Binomial",1950:2005)
sdm_guloe2<-data.frame(gulo_e_lpi$ID,as.character(gulo_e_lpi$Binomial) ,sdm_guloe)
colnames(sdm_guloe2)<-c( "ID","Binomial", 1950:2005)
sdm_lepus2<-data.frame(lepus_lpi$ID,as.character(lepus_lpi$Binomial) ,sdm_lepus)
colnames(sdm_lepus2)<-c( "ID","Binomial",1950:2005)
sdm_wtd2<-data.frame(wtd_lpi$ID, as.character(wtd_lpi$Binomial), sdm_wtd)
colnames(sdm_wtd2)<-c( "ID","Binomial",1950:2005)
sdm_rangifere2<-data.frame(rang_e_lpi$ID, as.character(rang_e_lpi$Binomial), sdm_range)
colnames(sdm_rangifere2)<-c( "ID","Binomial",1950:2005)
sdm_rangifern2<-data.frame(rang_n_lpi$ID, as.character(rang_n_lpi$Binomial), sdm_rangn)
colnames(sdm_rangifern2)<-c( "ID","Binomial",1950:2005)
sdm_polar2<-data.frame(polar_n_lpi$ID, as.character(polar_n_lpi$Binomial), sdm_polar)
colnames(sdm_polar2)<-c( "ID","Binomial",1950:2005)
sdm_phaco2<-data.frame(phaco_lpi$ID, as.character(phaco_lpi$Binomial), sdm_phaco)
colnames(sdm_phaco2)<-c( "ID","Binomial",1950:2005)
sdm_wilde2<-data.frame(wilde_lpi$ID, as.character(wilde_lpi$Binomial), sdm_wilde)
colnames(sdm_wilde2)<-c( "ID","Binomial",1950:2005)
sdm_py_chamois2<-data.frame(py_chamois_lpi$ID, as.character(py_chamois_lpi$Binomial), sdm_py_chamois)
colnames(sdm_py_chamois2)<-c( "ID","Binomial",1950:2005)
sdm_roe_deer2<-data.frame(roe_deer_lpi$ID, as.character(roe_deer_lpi$Binomial), sdm_roe_deer)
colnames(sdm_roe_deer2)<-c( "ID","Binomial",1950:2005)




d_out<-matrix(ncol = 56)

for (i in 1:length(sdm_harte)){
  
  df<-sdm_harte[[i]]
  if(is.matrix(df)){
    y<-colMeans(df)
  } else {
      y<-df
    }
  d_out<-rbind(y,d_out)
  print(i)  
  }

d_out<-d_out[-nrow(d_out),]

sdm_harte2<-data.frame(harte_lpi$ID, as.character(harte_lpi$Binomial), d_out)
colnames(sdm_harte2)<-c( "ID","Binomial",1950:2005)


d_out<-matrix(ncol = 56)

for (i in 1:length(sdm_giraffe)){
  
  df<-sdm_giraffe[[i]]
  if(is.matrix(df)){
    y<-colMeans(df)
  } else {
    y<-df
  }
  d_out<-rbind(y,d_out)
  print(i)  
}

d_out<-d_out[-nrow(d_out),]


sdm_giraffe2<-data.frame(giraffe_lpi$ID, as.character(giraffe_lpi$Binomial), d_out)
colnames(sdm_giraffe2)<-c( "ID","Binomial",1950:2005)

d_out<-matrix(ncol = 56)

for (i in 1:length(sdm_kobus)){
  
  df<-sdm_kobus[[i]]
  if(is.matrix(df)){
    y<-colMeans(df)
  } else {
    y<-df
  }
  d_out<-rbind(y,d_out)
  print(i)  
}

d_out<-d_out[-nrow(d_out),]


sdm_kobus2<-data.frame(kobus_lpi$ID, as.character(kobus_lpi$Binomial), d_out)
colnames(sdm_kobus2)<-c( "ID","Binomial",1950:2005)


d_out<-matrix(ncol = 56)

for (i in 1:length(sdm_equus)){
  
  df<-sdm_equus[[i]]
  if(is.matrix(df)){
    y<-colMeans(df)
  } else {
    y<-df
  }
  d_out<-rbind(y,d_out)
  print(i)  
}

d_out<-d_out[-nrow(d_out),]


sdm_equus2<-data.frame(equus_lpi$ID, as.character(equus_lpi$Binomial), d_out)
colnames(sdm_equus2)<-c( "ID","Binomial",1950:2005)



library(reshape2)

sdm_cervus_melt<-melt(sdm_cervus2, id.vars = c("ID", "Binomial"))
colnames(sdm_cervus_melt)<-c("ID","Binomial", "Year", "HSI")

sdm_capra_melt<-melt(sdm_capra2, id.vars = c("ID", "Binomial"))
colnames(sdm_capra_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_ursus_melt<-melt(sdm_ursus2, id.vars = c("ID", "Binomial"))
colnames(sdm_ursus_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_gulon_melt<-melt(sdm_gulon2, id.vars = c("ID", "Binomial"))
colnames(sdm_gulon_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_guloe_melt<-melt(sdm_guloe2, id.vars = c("ID", "Binomial"))
colnames(sdm_guloe_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_lepus_melt<-melt(sdm_lepus2, id.vars = c("ID", "Binomial"))
colnames(sdm_lepus_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_wtd_melt<-melt(sdm_wtd2, id.vars = c("ID", "Binomial"))
colnames(sdm_wtd_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_range_melt<-melt(sdm_rangifere2, id.vars = c("ID", "Binomial"))
colnames(sdm_range_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_rangn_melt<-melt(sdm_rangifern2, id.vars = c("ID", "Binomial"))
colnames(sdm_rangn_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_polar_melt<-melt(sdm_polar2, id.vars = c("ID", "Binomial"))
colnames(sdm_polar_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_harte_melt<-melt(sdm_harte2, id.vars = c("ID", "Binomial"))
colnames(sdm_harte_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_phaco_melt<-melt(sdm_phaco2, id.vars = c("ID", "Binomial"))
colnames(sdm_phaco_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_giraffe_melt<-melt(sdm_giraffe2, id.vars = c("ID", "Binomial"))
colnames(sdm_giraffe_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_kobus_melt<-melt(sdm_kobus2, id.vars = c("ID", "Binomial"))
colnames(sdm_kobus_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_wilde_melt<-melt(sdm_wilde2, id.vars = c("ID", "Binomial"))
colnames(sdm_wilde_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_py_chamois_melt<-melt(sdm_py_chamois2, id.vars = c("ID", "Binomial"))
colnames(sdm_py_chamois_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_roe_deer_melt<-melt(sdm_roe_deer2, id.vars = c("ID", "Binomial"))
colnames(sdm_roe_deer_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_roe_deer_melt<-melt(sdm_roe_deer2, id.vars = c("ID", "Binomial"))
colnames(sdm_roe_deer_melt)<-c("ID","Binomial" , "Year", "HSI")

sdm_equus_melt<-melt(sdm_equus2, id.vars = c("ID", "Binomial"))
colnames(sdm_equus_melt)<-c("ID","Binomial" , "Year", "HSI")


#i love fiona

#sdm_values
all_sp<-rbind(sdm_cervus_melt, sdm_capra_melt, sdm_ursus_melt, sdm_gulon_melt, 
              sdm_guloe_melt, sdm_lepus_melt, sdm_wtd_melt, sdm_range_melt, sdm_rangn_melt, 
              sdm_polar_melt, sdm_harte_melt,sdm_phaco_melt, sdm_giraffe_melt, sdm_kobus_melt,
              sdm_wilde_melt, sdm_py_chamois_melt, sdm_roe_deer_melt, sdm_equus_melt)


#write.csv(all_sp, "all_species_sdm.csv")


lpi<-read.csv("LPI_pops_20160523_edited.csv")
all_sp<-read.csv("all_species_sdm.csv")

all_sp<-all_sp[,-1]

#lpi values 


library(taRifx)
library(plyr)
library(mgcv)
library(zoo)

sp_lpi<-lpi[lpi$ID %in% all_sp$ID,] 

pops<-sp_lpi[,c(1,65:120)]

colnames(pops)<-c("ID", 1950:2005)

colnames(pops)[2:ncol(pops)]<-paste("Year", 1950:2005, sep="_")
pops[pops=="NULL"]<-NA
pops$rep_id<-"Observed"
pops$md_id<-"Observed"

popsm<-as.matrix(pops)

gam_lpi<-function(x){
  #subsetting the population data by each population
  
  df<-pops[pops$ID == x,]
  spid = df[2:length(df)]                     #subsetting only the dates
  names(spid)<-1950:2005              #renaming the date column names as R doesn't like numbered column names
  spid<-as.numeric(spid)
 
  if(sum(!is.na(spid)) >2){
  pop_datab <- (!is.na(spid) )
  points = sum(pop_datab)
  id<-x[1]
  Date<-1950:2005
  spidt<-destring(t(spid))
  time<-length(min(which(!is.na(spidt))):max(which(!is.na(spidt))))
  missing<-time-points
  
  Year<-Date[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  Population<-spidt[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  Population[Population == 0] <- mean(Population, na.rm=TRUE)*0.01 #if a population is zero one year thhis is replaced with 1% of the average population estimate - because you can log zeros
  
  df<-data.frame(Year,Population)
  
  if (sum(na.omit(df$Population<1))>0) {
    df$Population<-df$Population+1
  }
  
  PopN = log10(df$Population)
  if (points >=6) {
    
    if (length(na.omit(PopN)) >=6) {
      SmoothParm = round(length(na.omit(PopN))/2)
    } else {
      SmoothParm=3
    }
    
    mg2<-mgcv:::gam(PopN ~ s(Year, k=SmoothParm), fx=TRUE)
    pv2 <- predict(mg2,df,type="response",se=TRUE)
    R_sq2<-summary(mg2)$r.sq
    model<-1
    pv2$fit[pv2$fit <= 0] <- NA
    
    lambda2<-diff(pv2$fit)
    
    ial<-data.frame(id, Year[-length(Year)], lambda2)
    colnames(ial)<-c("ID", "Year", "Lambdas")
  } else {
    
    lint<-na.approx(PopN)
    lint[lint<=0] <- NA
    lambda2<-diff(lint)
    ial<-data.frame(id, Year[-length(Year)], lambda2)
    colnames(ial)<-c("ID", "Year", "Lambdas")
  }
  } else {
    ial<-data.frame(x[1], NA, NA)
    colnames(ial)<-c("ID", "Year", "Lambdas")
  }
  print(x[1])
  return(ial)
}

gam_lpi_r<-lapply(pops$ID,  gam_lpi)
gam_r<-do.call( "rbind", gam_lpi_r)

gam_r<-gam_r[gam_r$Year <=2005,]

fill<-data.frame(rep(pops$ID, each=length(1950:2005)), 1950:2005)
colnames(fill)<-c("ID", "Year")

all_year_ab<-join(fill, gam_r, type="right")

all_year_ab$ID<-as.numeric(as.character(all_year_ab$ID))



#joining sdm and lambdas
all_sdm_lambdas<-join(all_year_ab,all_sp, type="right")
all_sdm_lambdas$ID<-as.factor(all_sdm_lambdas$ID)
plot(all_sdm_lambdas$HSI,all_sdm_lambdas$Lambdas)

library(lme4)
m1<-lmer(Lambdas ~ HSI + (1|ID) +(1|Binomial), data = all_sdm_lambdas)
summary(m1)

#all_sdm_lambdas<-all_sdm_lambdas[complete.cases(all_sdm_lambdas),]

# library(ggplot2)
# ggplot(all_sdm_lambdas, aes(x = HSI, y = Lambdas, group = ID, colour = Binomial))+
#   geom_point()+
#   geom_smooth(se =FALSE)+
#   theme(legend.position="none")+
#   ylim(-0.5, 0.5)+
#   facet_wrap(~ID)
# 


####Abundance
pops<-pops[,1:(ncol(pops)-2)]
colnames(pops)[2:ncol(pops)]<-1950:2005
melt_lpi<-data.table:::melt(pops, id.vars = "ID")
colnames(melt_lpi)<-c("ID", "Year","Abundance")
melt_lpi$Year<-as.numeric(as.character(melt_lpi$Year))
melt_lpi$Abundance<-as.numeric(melt_lpi$Abundance)


all_ab_sdm_lambdas<-join(all_sdm_lambdas, melt_lpi)



plot(all_ab_sdm_lambdas$HSI,log10(all_ab_sdm_lambdas$Abundance+1))

lmer(log10(Abundance+1)~  HSI+(1|ID)+ (1|Binomial), data = all_ab_sdm_lambdas)

lmer(log10(Abundance+1)~  HSI+(1|ID), data = all_ab_sdm_lambdas)


#####plotting trends
library(ggplot2)
library(dplyr)

all_sdm_lambdas<-all_sdm_lambdas %>%
  group_by(ID)%>%
  mutate(HSI_Lambdas = c(diff(log10(HSI)),NA))

##Deer
ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Cervus_elaphus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Cervus_elaphus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Capra_ibex",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Capra_ibex",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Gulo_gulo",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Gulo_gulo",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Lepus_americanus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Lepus_americanus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Odocoileus_virginianus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Odocoileus_virginianus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Rangifer_tarandus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Rangifer_tarandus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Ursus_arctos",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Ursus_arctos",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Ursus_maritimus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Ursus_maritimus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Alcelaphus_buselaphus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Alcelaphus_buselaphus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Capreolus_capreolus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Capreolus_capreolus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Connochaetes_taurinus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Connochaetes_taurinus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Equus_burchellii",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Equus_burchellii",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Giraffa_camelopardalis",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Giraffa_camelopardalis",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Kobus_ellipsiprymnus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Kobus_ellipsiprymnus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Phacochoerus_africanus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Phacochoerus_africanus",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Rupicapra_pyrenaica",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Rupicapra_pyrenaica",], aes(x = Year, y = HSI_Lambdas, group = ID), colour = "blue")+
  facet_grid(.~ID)


#removing a zebra population which is outside of the HSM - perhaps in the sea?
all_sdm_lambdas<-all_sdm_lambdas%>%
  group_by(ID)%>%
  mutate(non_na_count = sum(!is.na(HSI)))%>%
  filter(non_na_count != 0)%>%
  ungroup()

#lost one population to the sea

smooth_gam<-function(x){
  
  a<-all_sdm_lambdas[all_sdm_lambdas$ID == x,]
  a$k = floor(0.5*nrow(a)) - 1
  
  if (a$k[1] >= 3){
    a$smooth_HSI = fitted.values(gam(HSI ~ s(Year, k = a$k[1]), data = a))
  } else {
    a$smooth_HSI =  fitted.values(lm(HSI ~ Year, data = a))
  }
  print(x)
  return(a)
}


smooth_sdm_lambdas<-lapply(unique(all_sdm_lambdas$ID), smooth_gam)

smooth_sdm_lambdas<-do.call("rbind", smooth_sdm_lambdas)



all_sdm_lambdas_new<-smooth_sdm_lambdas %>%
  group_by(ID)%>%
  mutate(smooth_HSI_Lambdas = c(diff(log10(smooth_HSI)),NA))


all_sdm_lambdas<-all_sdm_lambdas_new[complete.cases(all_sdm_lambdas_new),]

#195 populations left - lost all the ones with data from 2005 onwards only


#removing populations with only one record (would have had two pre-lambda-ing)
keep_id<-all_sdm_lambdas%>%
  group_by(ID)%>%
  #summarise(non_na_count = sum(!is.na(col_2)))%>%
  count(.) %>%
  filter(n>=5)


all_sdm_lambdas<-all_sdm_lambdas[all_sdm_lambdas$ID %in% keep_id$ID,]

populations<-unique(all_sdm_lambdas$ID)


all_sdm_lambdas%>%
  group_by(ID)%>%
  count()

ccf_get<-function(x){

    sdm_x<-all_sdm_lambdas[all_sdm_lambdas$ID == populations[x],]
    sdm_ccf<-ccf(sdm_x$HSI_Lambdas,sdm_x$Lambdas, type="correlation",lag.max=5, plot = F)
    sdm_ccf_smooth<-ccf(sdm_x$smooth_HSI_Lambdas,sdm_x$Lambdas,type="correlation")
    

    lag_n5_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -5)]
    lag_n4_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -4)]
    lag_n3_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -3)]
    lag_n2_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -2)]
    lag_n1_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -1)]
    lag_0_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 0)]
    nused<-sdm_ccf$n.used
    # lag_1_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 1)]
    # lag_2_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 2)]
    # lag_3_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 3)]
    # lag_4_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 4)]
    # lag_5_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 5)]
    
    sm_n5_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -5)]
    sm_n4_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -4)]
    sm_n3_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -3)]
    sm_n2_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -2)]
    sm_n1_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -1)]
    sm_0_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == 0)]
    
    # sm_1_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == 1)]
    # sm_2_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == 2)]
    # sm_3_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == 3)]
    # sm_4_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == 4)]
    # sm_5_sdm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == 5)]
    
    
    empty_check<-function(x){
      if(length(x) ==0){
       x<-NA 
      }
      return(x)
    }
    
    # lags<-list(lag_n5_sdm,lag_n4_sdm,lag_n3_sdm,lag_n2_sdm, lag_n1_sdm,lag_0_sdm, lag_1_sdm, lag_2_sdm, lag_3_sdm,
    #            lag_4_sdm,lag_5_sdm,sm_n5_sdm,sm_n4_sdm,sm_n3_sdm,sm_n2_sdm,sm_n1_sdm,sm_0_sdm,sm_1_sdm,sm_2_sdm,
    #            sm_3_sdm,sm_4_sdm,sm_5_sdm)
    
    lags<-list(lag_0_sdm, lag_n1_sdm, lag_n2_sdm, lag_n3_sdm,lag_n4_sdm,lag_n5_sdm,sm_0_sdm,sm_n1_sdm,sm_n2_sdm,
               sm_n3_sdm,sm_n4_sdm,sm_n5_sdm)
    
    lags_out<-t(data.frame(unlist(lapply(lags, empty_check))))
   
 
    ccf_out<-data.frame(unique(sdm_x$ID), unique(sdm_x$Binomial),nused,lags_out)
    # colnames(ccf_out)<-c("ID","Binomial","sdm_n5","sdm_n4","sdm_n3", "sdm_n2","sdm_n1", "sdm_0", "sdm_1","sdm_2",
    #                      "sdm_3","sdm_4","sdm_5","sm_sdm_n5","sm_sdm_n4","sm_sdm_n3", "sm_sdm_n2","sm_sdm_n1", "sm_sdm_0", "sm_sdm_1",
    #                      "sm_sdm_2","sm_sdm_3","sm_sdm_4","sm_sdm_5")
    
    colnames(ccf_out)<-c("ID","Binomial","N_used" ,"sdm_0", "sdm_n1","sdm_n2","sdm_n3","sdm_n4","sdm_n5","sm_sdm_0", "sm_sdm_n1",
                       "sm_sdm_n2","sm_sdm_n3","sm_sdm_n4","sm_sdm_n5")
    rownames(ccf_out) <- NULL
    print(x)
    return(ccf_out)
    
}
#p values from max cor
#https://stackoverflow.com/questions/38173544/how-to-calculate-p-values-from-cross-correlation-function-in-r

ccf_scores<-lapply(1:length(populations), ccf_get)

#ccf_all<-ccf_scores[lapply(ccf_scores,length)>0] 

ccf_df<-do.call( "rbind",ccf_scores)
#colnames(ccf_df)<-c("ID","Binomial" ,"sdm lag -5", "sdm lag -4", "sdm lag -3", "sdm lag -2","sdm lag -1",
#                    "sdm lag 0", "sdm lag 1", "sdm lag 2", "sdm lag 3", "sdm lag 4", "sdm lag 5","sm lag -5",
#                    "sm lag -4", "sm lag -3", "sm lag -2","sm lag -1", "sm lag 0", "sm lag 1", "sm lag 2",
#                    "sm lag 3", "sm lag 4", "sm lag 5")

colnames(ccf_df)<-c("ID","Binomial", "N_used","sdm lag 0", "sdm lag -1", "sdm lag -2", "sdm lag -3", "sdm lag -4", "sdm lag -5",
                    "sm lag 0", "sm lag -1", "sm lag -2","sm lag -3", "sm lag -4", "sm lag -5")


#colnames(ccf_df)<-c("ID","Binomial","sdm_n3", "sdm_n2","sdm_n1", "sdm_0", "sdm_1","sdm_2","sdm_3")

library(reshape2)
ccf_melt<-melt(ccf_df, id.vars = c("ID", "Binomial", "N_used"))
ccf_melt<-unique(ccf_melt)



ccf_melt<-ccf_melt[complete.cases(ccf_melt),]

#

ccf_melt$Binomial <- factor(ccf_melt$Binomial, levels = c("Capra_ibex","Connochaetes_taurinus","Ursus_arctos", 
                                                        "Phacochoerus_africanus", "Giraffa_camelopardalis", 
                                                        "Alcelaphus_buselaphus", "Equus_burchellii",
                                                        "Ursus_maritimus", "Rupicapra_pyrenaica",
                                                        "Cervus_elaphus", "Rangifer_tarandus",
                                                        "Capreolus_capreolus", "Lepus_americanus", 
                                                        "Kobus_ellipsiprymnus", "Odocoileus_virginianus", "Gulo_gulo"))


ccf_melt0<-ccf_melt[ccf_melt$variable == "sdm lag 0",]

#cross correlation function between habitat suitability and lambdas
ggplot(ccf_melt0,aes(x=Binomial, y=value, group = Binomial, fill =Binomial)) + 
  geom_boxplot()+
  geom_point(size = 3)+
  geom_hline(yintercept=0, linetype = "dashed")+
  theme_bw()+
  ylim(-1,1)+
  xlab("Species")+
  ylab("Correlation Coefficient")+
  theme(legend.position="none",axis.text=element_text(size=16),
                                      axis.title=element_text(size=20))+
  scale_x_discrete(labels = c("Alpine\nibex","Blue\nwildebeest", "Brown\nbear","Common\nwarthog", "Giraffe",
                              "Hartebeest","Plain's\nzebra", 
                              "Polar\nbear","Pyrenean\nchamois","Red\ndeer",  "Reindeer", "Roe\ndeer", 
                              "Snowshoe\nhare","Waterbuck", "White-tailed\ndeer", "Wolverine"))



ggplot(data =ccf_melt0,aes(x=value)) + 
  geom_histogram(binwidth = 0.05,fill="grey50", col="white")+
  #scale_x_continuous(breaks = seq(-0.7,0.7,0.05), lim = c(-0.7,0.7))+
  xlab("Growth Rate - Habitat Suitability Correlation Coefficient")+
  ylab("Frequency")+
  theme_bw()



 ccf_melt_sm0<-ccf_melt[ccf_melt$variable == "sm lag 0",]

#cross correlation function between habitat suitability and lambdas
ggplot(ccf_melt_sm0,aes(x=Binomial, y=value, group = Binomial, fill = Binomial)) + 
  geom_boxplot()+
  geom_point(size = 3)+
  geom_hline(yintercept=0, linetype = "dashed")+
  theme_bw()+
  ylim(-1,1)+
  xlab("Species")+
  ylab("Correlation Coefficient")+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=20))+
  scale_x_discrete(labels = c("Alpine\nibex","Blue\nwildebeest", "Brown\nbear","Common\nwarthog", "Giraffe",
                              "Hartebeest","Plain's\nzebra", 
                              "Polar\nbear","Pyrenean\nchamois","Red\ndeer",  "Reindeer", "Roe\ndeer", 
                              "Snowshoe\nhare","Waterbuck", "White\ntailed\ndeer", "Wolverine"))



#plot - hsi lambdas against pop lambdas?


ccf_melt1<-ccf_melt[ccf_melt$variable == "sdm lag -1",]

#cross correlation function between habitat suitability and lambdas
ggplot(ccf_melt1,aes(x=Binomial, y=value, group = Binomial, colour = Binomial)) + 
  geom_violin()+
  geom_jitter(width = 0.15)+
  geom_hline(yintercept=0)


ggplot(data =ccf_melt1,aes(x=value)) + 
  geom_histogram(binwidth = 0.05,fill="grey50", col="white")+
  #scale_x_continuous(breaks = seq(-0.7,0.7,0.05), lim = c(-0.7,0.7))+
  xlab("Growth Rate - Habitat Suitability Correlation Coefficient")+
  ylab("Frequency")+
  theme_bw()




#######unsmoothed hab suit


ccf_max<-ccf_melt%>%
  group_by(ID)%>%
  filter(grepl("sdm",variable))%>%
  mutate(max_value = max(value))%>%
  dplyr:::select(ID, Binomial, N_used,value,lag,max_value)%>%
  filter(lag >= -5 & value == max_value)

table(ccf_max$lag)

pvals<-(2 * (1 - pnorm(abs(ccf_max$max_value), mean = 0, sd = 1/sqrt(ccf_max$N_used))))


ccf_max$Binomial <- factor(ccf_max$Binomial, levels = c("Capra_ibex","Connochaetes_taurinus","Ursus_arctos", 
                                                   "Phacochoerus_africanus", "Giraffa_camelopardalis", 
                                                   "Alcelaphus_buselaphus", "Equus_burchellii",
                                                   "Ursus_maritimus", "Rupicapra_pyrenaica",
                                                   "Cervus_elaphus", "Rangifer_tarandus",
                                                   "Capreolus_capreolus", "Lepus_americanus", 
                                                   "Kobus_ellipsiprymnus", "Odocoileus_virginianus", "Gulo_gulo"))

ccf_max$lag<-as.factor(ccf_max$lag)

ggplot(ccf_max,aes(x=Binomial, y=max_value, group = Binomial, fill =Binomial)) + 
  geom_boxplot()+
  geom_point(size = 3)+
  #geom_point(data = ccf_max,aes(colour = lag), size = 3)+
  #geom_jitter(width = 0.15)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()+
  ylim(-1,1)+
  ylab("Maximum Coefficient Value\n (0-5 year lag)")+
  xlab("")+
  theme(legend.position="none",axis.text=element_text(size=17),
        axis.title=element_text(size=20))+
  scale_x_discrete(labels = c("Alpine\nibex","Blue\nwildebeest", "Brown\nbear","Common\nwarthog", "Giraffe",
                              "Hartebeest","Plain's\nzebra", 
                              "Polar\nbear","Pyrenean\nchamois","Red\ndeer",  "Reindeer", "Roe\ndeer", 
                              "Snowshoe\nhare","Waterbuck", "White\ntailed\ndeer", "Wolverine"))


out<-strsplit(as.character(ccf_melt$variable), " ")
lag_id<-do.call("rbind", out)

ccf_melt$lag<-as.numeric(lag_id[,3])

#smoothed hab suit


ccf_max_sm<-ccf_melt%>%
  group_by(ID)%>%
  filter(grepl("sdm",variable))%>%
  mutate(max_value = max(value))%>%
  select(ID, Binomial, N_used,value,lag,max_value)%>%
  filter(lag >= -5 & value == max_value)

ccf_max_sm$Binomial <- factor(ccf_max_sm$Binomial, levels = c("Capra_ibex","Connochaetes_taurinus","Ursus_arctos", 
                                                        "Phacochoerus_africanus", "Giraffa_camelopardalis", 
                                                        "Alcelaphus_buselaphus", "Equus_burchellii",
                                                        "Ursus_maritimus", "Rupicapra_pyrenaica",
                                                        "Cervus_elaphus", "Rangifer_tarandus",
                                                        "Capreolus_capreolus", "Lepus_americanus", 
                                                        "Kobus_ellipsiprymnus", "Odocoileus_virginianus", "Gulo_gulo"))



table(ccf_max_sm$lag)

pvals_sm<-(2 * (1 - pnorm(abs(ccf_max_sm$max_value), mean = 0, sd = 1/sqrt(ccf_max_sm$N_used))))

ggplot(ccf_max_sm,aes(x=Binomial, y=max_value, group = Binomial, fill = Binomial)) + 
  geom_boxplot()+
  geom_point(size = 3)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()+
  ylim(-1,1)+
  ylab("Maximum Coefficient Value\n (0-5 year lag)")+
  xlab("Species")+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=20))+
  scale_x_discrete(labels = c("Alpine\nibex","Blue\nwildebeest", "Brown\nbear","Common\nwarthog", "Giraffe",
                              "Hartebeest","Plain's\nzebra", 
                              "Polar\nbear","Pyrenean\nchamois","Red\ndeer",  "Reindeer", "Roe\ndeer", 
                              "Snowshoe\nhare","Waterbuck", "White-tailed\ndeer", "Wolverine"))


######rmse


rmse_get<-function(x){
  
  sdm_x<-all_sdm_lambdas[all_sdm_lambdas$ID == populations[x],]
  sdm_rmse<-rmse(sdm_x$Lambdas,sdm_x$HSI_Lambdas)
  sdm_rmse_smooth<-rmse(sdm_x$Lambdas,sdm_x$smooth_HSI_Lambdas)

    empty_check<-function(x){
    if(length(x) ==0){
      x<-NA 
    }
    return(x)
  }
  
  rmses<-list(sdm_rmse, sdm_rmse_smooth)
  rmses_out<-t(data.frame(unlist(lapply(rmses, empty_check))))
  rmse_out<-data.frame(unique(sdm_x$ID), unique(sdm_x$Binomial),rmses_out)
  # colnames(ccf_out)<-c("ID","Binomial","sdm_n5","sdm_n4","sdm_n3", "sdm_n2","sdm_n1", "sdm_0", "sdm_1","sdm_2",
  #                      "sdm_3","sdm_4","sdm_5","sm_sdm_n5","sm_sdm_n4","sm_sdm_n3", "sm_sdm_n2","sm_sdm_n1", "sm_sdm_0", "sm_sdm_1",
  #                      "sm_sdm_2","sm_sdm_3","sm_sdm_4","sm_sdm_5")
  colnames(rmse_out)<-c("ID","Binomial", "rmse", "rmse_smooth")
  rownames(rmse_out) <- NULL
  print(x)
  return(rmse_out)
}

rmse_scores<-lapply(1:length(populations), rmse_get)

#ccf_all<-ccf_scores[lapply(ccf_scores,length)>0] 

rmse_df<-do.call( "rbind",rmse_scores)








###sum lambdas






sum_lambdas<-all_sdm_lambdas%>%
  group_by(ID) %>%
  mutate(sum_Lambdas = sum(Lambdas), sum_HSI_Lambdas = sum(HSI_Lambdas), sum_smooth_HSI_Lambdas = sum(smooth_HSI_Lambdas))%>%
  select(ID, Binomial, sum_Lambdas, sum_HSI_Lambdas, sum_smooth_HSI_Lambdas)%>%
  distinct()


plot(sum_lambdas$sum_Lambdas, sum_lambdas$sum_HSI_Lambdas)

plot(sum_lambdas$sum_Lambdas, sum_lambdas$sum_smooth_HSI_Lambdas)

region_lpi<-data.frame(lpi$ID, lpi$Region)
colnames(region_lpi)<-c("ID", "Region")

region_lpi<-region_lpi[region_lpi$ID %in% sum_lambdas$ID,]

sum_lambdas$ID<-as.integer(as.character(sum_lambdas$ID))

sum_lambdas_reg<-join(region_lpi, sum_lambdas)

ggplot(sum_lambdas_reg, aes(x = sum_HSI_Lambdas, y = sum_Lambdas, group =Region, colour = Region))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlim(-2,2)

# +
#   facet_wrap(.~Region)

ggplot(sum_lambdas, aes(x = sum_smooth_HSI_Lambdas, y = sum_Lambdas, group = Binomial, colour = Binomial))+
  geom_point()+
  facet_wrap(.~Binomial)


sum_lambdas$Lam_HSI<-ifelse(sum_lambdas$sum_Lambdas * sum_lambdas$sum_HSI_Lambdas >0, 1,0)

sum_lambdas$Lam_sm_HSI<-ifelse(sum_lambdas$sum_Lambdas * sum_lambdas$sum_smooth_HSI_Lambdas >0, 1,0)


#rmse 



auc<-read.csv("AUC_all_species.csv")

mean_auc_df<-auc %>%
  group_by(Binomial, Model) %>%
  mutate(mean_auc = mean(AUC), sd_auc = sd(AUC)) %>%
  ungroup()%>%
  dplyr::select(Binomial, Model, mean_auc, sd_auc) %>%
  distinct()

write.csv(mean_auc_df, "Mean_AUC_all_species.csv")


auc %>%
 group_by(Binomial) %>%
  mutate(mean_auc = mean(AUC), sd_auc = sd(AUC)) %>%
  ungroup()%>%
  dplyr::select(Binomial, mean_auc, sd_auc) %>%
  distinct()










