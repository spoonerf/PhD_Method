
temp<-read.csv("All_LPI_Mean_Temp_Slope.csv")
luc<-read.csv("LUC_distance_all.csv")
LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")
Realm<-read.csv("selected_pops_Ecoregion.csv")
Realm<-Realm[,c("ID", "WWF_REALM2")]

pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_18.csv")
EurHil<-read.csv("Europe_HILDA_5_year_pops.csv")  # data from Euro-centric analysis

temp<-temp[,c("ID", "Estimate")]
LPI<-LPI[,c("ID","Binomial","Common_name","Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat", "T_realm", "FW_realm")]

df<-merge(merge(temp,luc, by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

df<-merge(df, Realm, by="ID")

nrow(df)

df2<-subset(df, !is.na(Estimate)&r_sq >= 0.5  & !is.na(LUC_dist)&length_time 
            >=5 & System!="Marine" &Specific_location == 1 & !is.na(lambda_sum) 
            & !is.na(lambda_mean) & Class=="Aves")


nrow(df2)
#write.csv(df2, "selected_pops_LPI.csv")

library(plyr)
#counting duplicates at each location
sp_dups<-data.frame(ddply(df2,.(Longitude,Latitude),nrow))
sp_dups$loc_id<-1:length(sp_dups$Longitude)
sp_dups_df<-merge(sp_dups, df2, by=c("Longitude","Latitude"))

library(data.table)
#dt = as.data.table(sp_dups_df)

parm_df<-sp_dups_df[,c("ID","Estimate", "LUC_dist")]  ##ID, land use, and climate
#for hilda data
parm_mat<-as.matrix(parm_df)
parm_scale<-scale(parm_mat[,c("Estimate", "LUC_dist")])       #use the scaling factors at the bottom of these to scale the rasters

centre_temp<-attr(parm_scale, 'scaled:center')[1]
centre_luc<-attr(parm_scale, 'scaled:center')[2]
scale_temp<-attr(parm_scale, 'scaled:scale')[1] 
scale_luc<-attr(parm_scale, 'scaled:scale')[2] 

parm_id<-parm_mat[,"ID"]

parm_df_scale<-data.frame(parm_id,parm_scale)

colnames(parm_df_scale)<-c("ID","mean_slope_scale", "change_rate_scale")

sp_df_scale<-merge(sp_dups_df, parm_df_scale, by="ID")

dt<-data.table(sp_df_scale)


#dt<-data.table(Euro)

land<-seq(min(dt$change_rate_scale),max(dt$change_rate_scale),length.out=100)
temp<-seq(min(dt$mean_slope_scale),max(dt$mean_slope_scale),length.out=100)

# df<-expand.grid(land,temp)
# colnames(df)<-c("change_rate_scale", "mean_slope_scale")

#unscaling the data
landus<-(land  * scale_luc) + centre_luc
tempus<-(temp * scale_temp) + centre_temp

dfus<-expand.grid(landus,tempus)
colnames(dfus)<-c("change_rate_scale", "mean_slope_scale")


library(lme4)
library(MuMIn)
library(raster)
source("rsquaredglmm.R")

R=299
AIC_m1= numeric(R)
AIC_m1a= numeric(R)
AIC_m1b= numeric(R)
AIC_m1c= numeric(R)
AIC_mnull= numeric(R)

marg_Rsq_m1 = numeric(R)
marg_Rsq_m1a = numeric(R)
marg_Rsq_m1b = numeric(R)
marg_Rsq_m1c = numeric(R)

cond_Rsq_m1 = numeric(R)
cond_Rsq_m1a = numeric(R)
cond_Rsq_m1b = numeric(R)
cond_Rsq_m1c = numeric(R)
cond_Rsq_mnull = numeric(R)

m1_w = numeric(R)
m1a_w = numeric(R)
m1b_w = numeric(R)
m1c_w = numeric(R)
mnull_w = numeric(R)

MTC_i = numeric(R)
LUC_i = numeric(R)
LUC_MTC_i = numeric(R)

int_av = numeric(R)
MTC_av = numeric(R)
LUC_av = numeric(R)
LUC_MTC_av = numeric(R)

MTC_avus = numeric(R)
LUC_avus = numeric(R)


for (i in 1:R) {
  dt2<-data.frame(dt[, ID[sample.int(.N, 1, TRUE)], by = loc_id])     #.N     signifies the number of rows when using data.table
  colnames(dt2)[2]<-"ID"
  sp_dups_df2<-sp_df_scale[sp_df_scale$ID %in% dt2$ID,]

  m1<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1T<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_dups_df2)
  
  m1a<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1aT<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=sp_dups_df2)
  
  m1b<-lmer(lambda_mean ~ change_rate_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1bT<-lmer(lambda_mean ~ change_rate_scale+(1|Binomial),data=sp_dups_df2)
  
  m1c<-lmer(lambda_mean ~ mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1cT<-lmer(lambda_mean ~ mean_slope_scale+(1|Binomial),data=sp_dups_df2)
  
  mnull<-lmer(lambda_mean ~ 1+(1|Binomial),data=sp_dups_df2, REML=F)
  mnullT<-lmer(lambda_mean ~ 1+(1|Binomial),data=sp_dups_df2)

  #AIC
  AIC_m1[i]<-AIC(m1)
  AIC_m1a[i]<-AIC(m1a)
  AIC_m1b[i]<-AIC(m1b)
  AIC_m1c[i]<-AIC(m1c)
  AIC_mnull[i]<-AIC(mnull)

  # #Weights

  msAICc <- model.sel(m1,m1a,m1b,m1c,mnull)
  msAICc$model<-rownames(msAICc)
  msAICc<-data.frame(msAICc)

  m1_w[i]<-subset(msAICc, model=="m1")$weight
  m1a_w[i]<-subset(msAICc, model=="m1a")$weight
  m1b_w[i]<-subset(msAICc, model=="m1b")$weight
  m1c_w[i]<-subset(msAICc, model=="m1c")$weight
  mnull_w[i]<-subset(msAICc, model=="mnull")$weight

  #Rsq
  models_list<-list(m1,m1a,m1b,m1c,mnull)
  modelsR<-lapply(models_list,rsquared.glmm)
  modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)

  marg_Rsq_m1[i]<-modelsRsq[1,4]
  marg_Rsq_m1a[i]<-modelsRsq[2,4]
  marg_Rsq_m1b[i]<-modelsRsq[3,4]
  marg_Rsq_m1c[i]<-modelsRsq[4,4]
  cond_Rsq_m1[i]<-modelsRsq[1,5]
  cond_Rsq_m1a[i]<-modelsRsq[2,5]
  cond_Rsq_m1b[i]<-modelsRsq[3,5]
  cond_Rsq_m1c[i]<-modelsRsq[4,5]
  cond_Rsq_mnull[i]<-modelsRsq[5,5]

  var_imp<-summary(model.avg(models_list))
  # MTC_i[i]<-var_imp$importance["mean_slope_scale"]
  # LUC_i[i]<-var_imp$importance["change_rate_scale"]
  # LUC_MTC_i[i]<-var_imp$importance["change_rate_scale:mean_slope_scale"]
  # #BM_i[i]<-var_imp$importance["Bodymass"]

  int_av[i]<-var_imp$coefmat.subset["(Intercept)","Estimate"]
  MTC_av[i]<-var_imp$coefmat.subset["mean_slope_scale","Estimate"]
  LUC_av[i]<-var_imp$coefmat.subset["change_rate_scale","Estimate"]
  LUC_MTC_av[i]<-var_imp$coefmat.subset["change_rate_scale:mean_slope_scale","Estimate"]
  
  MTC_avus[i]<-(MTC_av[i] * scale_temp) + centre_temp
  LUC_avus[i]<-(LUC_av[i] * scale_luc) + centre_luc
  
  mav<-model.avg(models_list)
  pred<-predict(mav, dfus, re.form=NA)
  pdf<-data.frame(dfus, pred)
  pred2<-matrix(pred, ncol=length(landus), nrow=length(tempus), byrow=T)
  #image(land,temp,pred2)
  # head(pred2)
  # pras<-raster(pred2, xmn=min(landus), xmx=max(landus), ymn=min(tempus), ymx=max(tempus))
  # 
  pred_pcnt<-(10^pred2) - 1    #for lambda sum this is total population change, for lambda mean this is average annual rate of change
  pcntras<-raster(pred_pcnt, xmn=min(landus), xmx=max(landus), ymn=min(tempus), ymx=max(tempus))
  
  file<-paste("predict_", i,".tif" ,sep="")
  writeRaster(pcntras, filename=file, overwrite=TRUE)
  print(i)
}


list<-list.files(path = getwd(), pattern = "predict.*\\.tif$")

br<-stack(list)

br_av<-mean(br)


plot(br_av, xlab="Land Use Change Distance", ylab="Annual Mean Temperature Change", main="Average Annual Population Change (%) - Birds")

writeRaster(br_av, "Global_Birds_Prediction_Average_Raster_Lambda_mean.tif", overwrite=TRUE)

####
AIC_df<-data.frame(cbind(AIC_m1,AIC_m1a, AIC_m1b, AIC_m1c, AIC_mnull))

AIC_del<-AIC_df[,c(1:4)] - AIC_df$AIC_mnull

colMeans(AIC_del)

AIC_int<-data.frame(AIC=AIC_m1)
AIC_luc<-data.frame(AIC=AIC_m1b)
AIC_mtc<-data.frame(AIC=AIC_m1c)
AIC_null<-data.frame(AIC=AIC_mnull)

AIC_int$model<-'Interacting'
AIC_luc$model<- 'Land Use Change'
AIC_mtc$model<- 'Mean Temp Change'
AIC_null$model<-'Null'
#and combine into your new data frame vegLengths
AIClengths <- rbind(AIC_int, AIC_null)
#AIClengths <- rbind(AIC_mtc, AIC_null)
#AIClengths <- rbind(AIC_luc, AIC_null)

library(ggplot2)
#now make your lovely plot
ggplot(AIClengths, aes(AIC, fill = model)) + geom_density(alpha = 0.2)




#####

Low<-(R+1)/40
High<-(R+1)-(R+1)/40 

mean_av<- c(mean(LUC_av),mean(MTC_av),mean(LUC_MTC_av), mean(LUC_avus), mean(MTC_avus))
lowCI_av<-c(sort(LUC_av)[Low], sort(MTC_av)[Low], sort(LUC_MTC_av)[Low], sort(LUC_avus)[Low], sort(MTC_avus)[Low])
highCI_av<-c(sort(LUC_av)[High], sort(MTC_av)[High], sort(LUC_MTC_av)[High],sort(LUC_avus)[High], sort(MTC_avus)[High] )

Variable<-c("LUC", "MTC", "LUC*MTC", "LUC_unscale", "MTC_unscale")

conf_av<-data.frame(rbind( lowCI_av, mean_av, highCI_av))
colnames(conf_av)<-Variable
conf_av

######



mean(marg_Rsq_m1)
mean(marg_Rsq_m1a)
mean(marg_Rsq_m1b)
mean(marg_Rsq_m1c)
mean(cond_Rsq_m1)
mean(cond_Rsq_m1a)
mean(cond_Rsq_m1b)
mean(cond_Rsq_m1c)
mean(cond_Rsq_mnull)
