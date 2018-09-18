#getting trends from correlative model in chapter 2

temp<-read.csv("All_LPI_All_Years_Nobuff_1931_moreLPI_end2005.csv")


body4<-read.csv("LPI_BodyMass_Amniote_Database_edit.csv")
body4<-body4[,-3]

hyde<-read.csv("Hyde_crop_pasture_annual_change.csv")

hyde2<-read.csv("Hyde_crop_pasture_annual_change_sum.csv")

temp<-merge(temp, hyde[,c(2,3)], by="ID")

LPI<-read.csv("LPI_pops_20160523_edited.csv")
LPI_bin<-LPI[,c("ID", "Binomial")]

body4<-merge(body4[,c(2:4)], LPI_bin, by="Binomial")

body4<-unique(body4)

pop<-read.csv("Global_Population_Trends_Rsq_Lambda_07_10_16.csv")

temp<-temp[,c("ID", "Estimate")]

LPI<-LPI[,c("ID","Binomial","Confidential","Common_name", "Order","Family", "Protected_status", "Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat",  "Confidential", "Forest","Savanna", "Shrubland", "Grassland", "Wetland_Inland", "Rocky_areas", "Caves", "Desert")]

df<-merge(merge(temp,body4[,c(2:4)], by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

dfd<-merge(df, hyde[,c(-1,-3)], by="ID")

head(dfd)

df2<-subset(dfd, !is.na(Estimate)  & r_sq >= 0.4999999 &length_time >=5& System!="Marine"
            &Specific_location == 1 &!is.na(both_change) & !is.na(Log_Body_Mass_g)
            & (Class=="Mammalia") & Protected_status != "Unknown"  & Protected_status != "Both")
nrow(df2)


df2[is.na(df2$lambda_mean),]$lambda_mean<-0


nrow(df2)

library(plyr)
#counting duplicates at each location
sp_dups<-data.frame(ddply(df2,.(Longitude,Latitude),nrow))
sp_dups$loc_id<-1:length(sp_dups$Longitude)
sp_dups_df<-merge(sp_dups, df2, by=c("Longitude","Latitude"))

library(data.table)
dt = as.data.table(sp_dups_df)

parm_df<-sp_dups_df[,c("ID","Estimate", "both_change", "Log_Body_Mass_g")]  ##ID, land use, and climate  use "LUC_dist" or "Nat_change" for purely annual change in summed primary, secondary and other

parm_mat<-as.matrix(parm_df)
parm_scale<-scale(parm_mat[,c("Estimate", "both_change", "Log_Body_Mass_g")])       #use the scaling factors at the bottom of these to scale the rasters

parm_id<-parm_mat[,"ID"]

parm_df_scale<-data.frame(parm_id,parm_scale)

colnames(parm_df_scale)<-c("ID","mean_slope_scale", "change_rate_scale", "Bodymass_scale")

sp_df_scale<-merge(sp_dups_df, parm_df_scale, by="ID")

dt<-data.table(sp_df_scale)

length(unique(dt$loc_id))

nrow(dt)

source("rsquaredglmm.R")

library(lme4) 

m0<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0f<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0a<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0af<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0b<-lmer(lambda_mean ~ change_rate_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0bf<-lmer(lambda_mean ~ change_rate_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0c<-lmer(lambda_mean ~ mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0cf<-lmer(lambda_mean ~ mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0d<-lmer(lambda_mean ~ Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0df<-lmer(lambda_mean ~ Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1f<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1a<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1af<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1b<-lmer(lambda_mean ~ change_rate_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1bf<-lmer(lambda_mean ~ change_rate_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)  

m1c<-lmer(lambda_mean ~ mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1cf<-lmer(lambda_mean ~ mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

mnull<-lmer(lambda_mean ~ 1+(1|Binomial)+(1|loc_id),data=dt, REML=F)

AIC(m0)



# #Weights
library(MuMIn)

msAICc <- model.sel(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull)
msAICc <- model.sel(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull,m0f,m0af,m0bf,m0cf,m0df,m1f,m1af,m1bf,m1cf)

msAIC<-msAICc
msAIC$model<-rownames(msAICc)
msAIC<-data.frame(msAIC)
msAIC

((10^msAIC[,c(1:5)]) - 1)*100


mAIC<-AIC(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull,m0f,m0af,m0bf,m0cf,m0df,m1f,m1af,m1bf,m1cf)

#Rsq
models_list<-list(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull)
models_list<-list(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull,m0f,m0af,m0bf,m0cf,m0df,m1f,m1af,m1bf,m1cf)

modelsR<-lapply(models_list,rsquared.glmm)
modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)
#rownames(modelsRsq)<-c("m0","m0a","m0b","m0c","m0d","m1","m1a","m1b","m1c","mnull")
rownames(modelsRsq)<-c("m0","m0a","m0b","m0c","m0d","m1","m1a","m1b","m1c","mnull","m0f","m0af","m0bf","m0cf","m0df","m1f","m1af","m1bf","m1cf")
modelsRsq

library(MuMIn)
var_imp<-summary(model.avg(msAICc))


mav<-model.avg(msAICc, subset =  cumsum(weight) <= .95)


new_data2<-new_data[rep(seq_len(nrow(new_data)), each=2),]

predict(m0, newdata=new_data2)

####extracting climate data
library(ncdf4)
library(raster)

#CR30s<-brick("cru_ts3.23.1931.1940.tmp.dat.nc")
CR40s<-brick("cru_ts3.23.1941.1950.tmp.dat.nc")
CR50s<-brick("cru_ts3.23.1951.1960.tmp.dat.nc")
CR60s<-brick("cru_ts3.23.1961.1970.tmp.dat.nc")
CR70s<-brick("cru_ts3.23.1971.1980.tmp.dat.nc")
CR80s<-brick("cru_ts3.23.1981.1990.tmp.dat.nc")
CR90s<-brick("cru_ts3.23.1991.2000.tmp.dat.nc")
CR00s<-brick("cru_ts3.23.2001.2010.tmp.dat.nc")
CR10s<-brick("cru_ts3.23.2011.2014.tmp.dat.nc")

plot(CR50s[[11]], main="Mean Temperature, Novermber 1951")

#LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")
LPI<-read.csv("LPI_pops_20160523_edited.csv")

LPIsp<-subset(LPI, Specific_location==1 & Binomial == "Cervus_elaphus" & Region =="Europe" )

CR<-stack(CR30s,CR40s,CR50s,CR60s,CR70s,CR80s,CR90s,CR00s,CR10s)

plot(CR[[1]])
points(LPIsp$Longitude, LPIsp$Latitude)

LPIsp<-LPIsp[LPIsp$ID != 4382,]

xy<-data.frame(LPIsp$Longitude, LPIsp$Latitude)

xy<-unique(xy)     #identifying unique locations to extract climate data from 

xy_df<-data.frame(xy)
colnames(xy_df)<-c("lon", "lat")
coordinates(xy_df) <- c("lon", "lat")

length(xy_df)

library(doParallel)

indices <- format(as.Date(names(CR), format = "X%Y.%m.%d"), format = "%Y")
indices <- as.numeric(indices)

#sum layers
#yearly_CR<- stackApply(CR[[min(which(indices == 1950)):max(which(indices == 2005))]], indices[min(which(indices == 1950)):max(which(indices == 2005))], fun = mean)


temp_out<-extract(CR, xy_df)

colnames(temp_out)<-indices
temp_melt<-melt(temp_out)
colnames(temp_melt)<-c("ID", "Year", "Temp")


annual_mean_temp<-temp_melt%>%
  group_by(ID, Year)%>%
  summarise(mean_temp = mean(Temp))


annual_mean_temp$ID<-rep(LPIsp$ID, each = (nrow(annual_mean_temp)/length(LPIsp$ID)))

annual_mean_temp<-annual_mean_temp%>%
  group_by(ID)%>%
  mutate(temp_lambdas =c(0,diff(log10(mean_temp))))



ggplot(annual_mean_temp, aes(x = Year, y = mean_temp, group=ID))+
  geom_line()+
  facet_wrap(.~ID)


ggplot(annual_mean_temp, aes(x = Year, y = temp_lambdas, group=ID))+
  geom_line()+
  facet_wrap(.~ID)
  


scale_centre<-0.016147418

scale_scale<-0.078936102 


ids<-dt[dt$ID %in% LPIsp$ID,]

annual_mean_temp$mean_slope_scale<-(annual_mean_temp$temp_lambdas - scale_centre)/scale_scale

annual_mean_temp<-annual_mean_temp[annual_mean_temp$ID %in% ids$ID & annual_mean_temp$Year >=1950 & annual_mean_temp$Year<=2005,]

annual_mean_temp$loc_id<-rep(ids$loc_id, each = (nrow(annual_mean_temp)/length(ids$ID)))

annual_mean_temp$Binomial<-"Cervus_elaphus"


annual_mean_temp$predicted_lambdas<-predict(m1c, annual_mean_temp)


all_year_ab$ID<-as.numeric(as.character(all_year_ab$ID))

all_year_ab_sub<-all_year_ab[all_year_ab$ID %in% ids$ID,]

ggplot()+
  geom_line(data = annual_mean_temp, aes(x = Year, y = predicted_lambdas, group=ID))+
  geom_line(data = all_year_ab_sub, aes(x = Year, y= Lambdas, group=ID), colour="red", fill="lightcoral")+    #real
  facet_wrap(.~ID)










