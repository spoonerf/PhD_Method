
lf<-list.files("C:/Users/Fiona/Documents/PhD/PhD_Method/Legion/snow_cervus_bias_new/")
files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]

library(raster)
cervus<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Legion/snow_cervus_bias_new/", files, sep = ""))
capra<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Legion/snow_capra_bias/", files, sep = ""))
ursus<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Legion/snow_bear_bias_new/", files, sep = ""))
gulo_n<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Gulo_gulo_bias_NAm/SDM_folder_bias/", files, sep =""))
gulo_e<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Gulo_gulo_bias_EU/SDM_folder_bias_EU/", files, sep =""))
lepus<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Lepus_americanus_bias/SDM_folder_bias/", files, sep=""))
wtd<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Odocoileus_virginianus_bias/SDM_folder_bias/",files, sep = ""))
rangifer_e<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Rangifer_tarandus_bias/SDM_folder_bias/", files, sep = ""))
rangifer_n<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Rangifer_tarandus_bias_NAm/SDM_folder_bias_NAm/", files, sep = ""))
polar_n<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Ursus_maritimus_bias/SDM_folder_bias/", files, sep = ""))
harte<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Alcelaphus_buselaphus_bias/SDM_folder_bias/", files, sep=""))
phaco<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Phacochoerus_africanus_bias/SDM_folder_bias/", files, sep=""))
giraffe<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Giraffa_camelopardalis_bias/SDM_folder_bias/", files, sep=""))  
kobus<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Kobus_ellipsiprymnus_bias/SDM_folder_bias/", files, sep=""))
wilde<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Connochaetes_taurinus_bias/SDM_folder_bias/", files, sep=""))
py_chamois<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Rupicapra_pyrenaica_bias/SDM_folder_bias/", files, sep=""))

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

sdm_harte<-raster:::extract(harte, cbind(harte_lpi$Longitude, harte_lpi$Latitude), buffer = 50000, fun = mean, na.rm = TRUE)
sdm_giraffe<-raster:::extract(giraffe, cbind(giraffe_lpi$Longitude, giraffe_lpi$Latitude), buffer= 50000, fun = mean, na.rm = TRUE)
sdm_kobus<-raster:::extract(kobus, cbind(kobus_lpi$Longitude, kobus_lpi$Latitude), buffer= 50000, fun = mean, na.rm = TRUE)

sdm_wilde<-raster:::extract(wilde, cbind(wilde_lpi$Longitude, wilde_lpi$Latitude), buffer= 50000, fun = mean, na.rm =TRUE)
sdm_py_chamois<-raster:::extract(py_chamois, cbind(py_chamois_lpi$Longitude, py_chamois_lpi$Latitude), buffer= 50000, fun = mean, na.rm =TRUE)
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
sdm_harte2<-data.frame(harte_lpi$ID, as.character(harte_lpi$Binomial), sdm_harte)
colnames(sdm_harte2)<-c( "ID","Binomial",1950:2005)
sdm_phaco2<-data.frame(phaco_lpi$ID, as.character(phaco_lpi$Binomial), sdm_phaco)
colnames(sdm_phaco2)<-c( "ID","Binomial",1950:2005)
sdm_giraffe2<-data.frame(giraffe_lpi$ID, as.character(giraffe_lpi$Binomial), sdm_giraffe)
colnames(sdm_giraffe2)<-c( "ID","Binomial",1950:2005)
sdm_kobus2<-data.frame(kobus_lpi$ID, as.character(kobus_lpi$Binomial), sdm_kobus)
colnames(sdm_kobus2)<-c( "ID","Binomial",1950:2005)
sdm_wilde2<-data.frame(wilde_lpi$ID, as.character(wilde_lpi$Binomial), sdm_wilde)
colnames(sdm_wilde2)<-c( "ID","Binomial",1950:2005)
sdm_py_chamois2<-data.frame(py_chamois_lpi$ID, as.character(py_chamois_lpi$Binomial), sdm_py_chamois)
colnames(sdm_py_chamois2)<-c( "ID","Binomial",1950:2005)



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
colnames(sdm_rangn_melt)<-c("ID","Binomial" , "Year", "HSI")

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

#i love fiona

#sdm_values
all_sp<-rbind(sdm_cervus_melt, sdm_capra_melt, sdm_ursus_melt, sdm_gulon_melt, 
              sdm_guloe_melt, sdm_lepus_melt, sdm_wtd_melt, sdm_range_melt, sdm_rangn_melt, 
              sdm_polar_melt, sdm_harte_melt,sdm_phaco_melt, sdm_giraffe_melt, sdm_kobus_melt,
              sdm_wilde_melt, sdm_py_chamois_melt)


#write.csv(all_sp, "all_species_sdm.csv")

all_sp<-read.csv("all_species_sdm.csv")

all_sp<-all_sp[,-1]

#lpi values 


library(taRifx)
library(plyr)
library(mgcv)
library(zoo)

sp_lpi<-lpi[lpi$ID%in% all_sp$ID,] 

pops<-sp_lpi[,c(1,65:120)]

colnames(pops)<-c("ID", 1950:2005)

colnames(pops)[2:ncol(pops)]<-paste("Year", 1950:2005, sep="_")
pops[pops=="NULL"]<-NA
pops$rep_id<-"Observed"
pops$md_id<-"Observed"

popsm<-as.matrix(pops)

gam_lpi<-function(x){
  #subsetting the population data by each population
  spid = x[2:(length(x)-2)]                     #subsetting only the dates
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
  
  #not sure what this does - adding a constant of 1 so that logging doesn't go weird?
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

gam_lpi_r<-apply(popsm,  1, gam_lpi)
gam_r<-do.call( "rbind", gam_lpi_r)

gam_r<-gam_r[gam_r$Year <=2005,]

fill<-data.frame(rep(pops$ID, each=length(1950:2005)), 1950:2005)
colnames(fill)<-c("ID", "Year")

all_year_ab<-join(fill, gam_r, type="right")

all_year_ab$ID<-as.numeric(as.character(all_year_ab$ID))



#joining sdm and lambdas
all_sdm_lambdas<-join(all_year_ab,all_sp)
all_sdm_lambdas$ID<-as.factor(all_sdm_lambdas$ID)
plot(all_sdm_lambdas$HSI,all_sdm_lambdas$Lambdas)

library(lme4)
m1<-lmer(Lambdas ~ HSI + (1|ID) +(1|Binomial), data = all_sdm_lambdas)
summary(m1)

all_sdm_lambdas<-all_sdm_lambdas[complete.cases(all_sdm_lambdas),]

library(ggplot2)
ggplot(all_sdm_lambdas, aes(x = HSI, y = Lambdas, group = ID, colour = Binomial))+
  geom_point()+
  geom_smooth(se =FALSE)+
  theme(legend.position="none")+
  ylim(-0.5, 0.5)+
  facet_wrap(~ID)



####Abundance
pops<-pops[,1:(ncol(pops)-2)]
colnames(pops)[2:ncol(pops)]<-1950:2005
melt_lpi<-data.table:::melt(pops, id.vars = "ID")
colnames(melt_lpi)<-c("ID", "Year","Abundance")
melt_lpi$Year<-as.numeric(as.character(melt_lpi$Year))
melt_lpi$Abundance<-as.numeric(melt_lpi$Abundance)


all_ab_sdm_lambdas<-join(all_sdm_lambdas, melt_lpi)



plot(log10(all_ab_sdm_lambdas$Abundance+1), all_ab_sdm_lambdas$Lambdas)

lmer(log10(Abundance+1)~  HSI+(1|ID)+ (1|Binomial), data = all_ab_sdm_lambdas)



#####plotting trends

all_sdm_lambdas<-all_sdm_lambdas[complete.cases(all_sdm_lambdas),]




##Deer
ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Cervus_elaphus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Cervus_elaphus",], aes(x = Year, y = HSI, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Capra_ibex",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Capra_ibex",], aes(x = Year, y = HSI, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Gulo_gulo",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Gulo_gulo",], aes(x = Year, y = HSI, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Lepus_americanus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Lepus_americanus",], aes(x = Year, y = HSI, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Odocoileus_virginianus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Odocoileus_virginianus",], aes(x = Year, y = HSI, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Rangifer_tarandus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Rangifer_tarandus",], aes(x = Year, y = HSI, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Ursus_arctos",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Ursus_arctos",], aes(x = Year, y = HSI, group = ID), colour = "blue")+
  facet_grid(.~ID)

ggplot()+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Ursus_maritimus",], aes(x = Year, y = Lambdas, group = ID), colour = "black")+
  geom_smooth(data = all_sdm_lambdas[all_sdm_lambdas$Binomial == "Ursus_maritimus",], aes(x = Year, y = HSI, group = ID), colour = "blue")+
  facet_grid(.~ID)



populations<-unique(all_sdm_lambdas$ID)

ccf_get<-function(x){

    sdm_x<-all_sdm_lambdas[all_sdm_lambdas$ID == populations[x],]
    sdm_ccf<-ccf(sdm_x$HSI, sdm_x$Lambdas, type="correlation")
    
    lag_n3_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -3)]
    lag_n2_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -2)]
    lag_n1_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -1)]
    lag_0_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 0)]
    lag_1_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 1)]
    lag_2_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 2)]
    lag_3_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 3)]
    
    empty_check<-function(x){
      if(length(x) ==0){
       x<-NA 
      }
      return(x)
      }
    lags<-list(lag_n3_sdm,lag_n2_sdm, lag_n1_sdm,lag_0_sdm, lag_1_sdm, lag_2_sdm, lag_3_sdm)
    lags_out<-t(data.frame(unlist(lapply(lags, empty_check))))
   
    
    # ccf_out<-data.frame(sdm_x$ID, sdm_x$Binomial, lag_n3_sdm,lag_n2_sdm, lag_n1_sdm, lag_0_sdm, lag_1_sdm,lag_2_sdm, lag_3_sdm)
   # colnames(ccf_out)<-c("ID","Binomial","sdm_n3", "sdm_n2","sdm_n1", "sdm_0", "sdm_1","sdm_2","sdm_3")
    
    ccf_out<-data.frame(unique(sdm_x$ID), unique(sdm_x$Binomial),lags_out)
    colnames(ccf_out)<-c("ID","Binomial","sdm_n3", "sdm_n2","sdm_n1", "sdm_0", "sdm_1","sdm_2","sdm_3")
    rownames(ccf_out) <- NULL
    print(x)
    return(ccf_out)
    
}

ccf_scores<-lapply(1:length(populations), ccf_get)

#ccf_all<-ccf_scores[lapply(ccf_scores,length)>0] 

ccf_df<-do.call( "rbind",ccf_scores)
#colnames(ccf_df)<-c("ID","ldd","SD","cnd lag -5", "sdm lag -5","cnd lag -4", "sdm lag -4","cnd lag -3", "sdm lag -3","cnd lag -2", "sdm lag -2","cnd lag -1", "sdm lag -1", "cnd lag 0", "sdm lag 0","cnd lag 1", "sdm lag 1","cnd lag 2", "sdm lag 2","cnd lag 3", "sdm lag 3","cnd lag 4", "sdm lag 4","cnd lag 5", "sdm lag 5")
colnames(ccf_df)<-c("ID","Binomial","sdm_n3", "sdm_n2","sdm_n1", "sdm_0", "sdm_1","sdm_2","sdm_3")

ccf_melt<-melt(ccf_df, id.vars = c("ID", "Binomial"))
ccf_melt<-unique(ccf_melt)

#

ccf_melt0<-ccf_melt[ccf_melt$variable == "sdm_0",]

#cross correlation function between habitat suitability and lambdas
ggplot(ccf_melt0,aes(x=Binomial, y=value, group = Binomial, colour = Binomial)) + 
  geom_violin()+
  geom_jitter(width = 0.2)


#plot - hsi lambdas against pop lambdas?











