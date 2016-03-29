
temp<-read.csv("All_LPI_Mean_Temp_Slope.csv")
luc<-read.csv("LUC_distance_all.csv")
LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")

pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_18.csv")
EurHil<-read.csv("Europe_HILDA_5_year_pops.csv")  # data from Euro-centric analysis

temp<-temp[,c("ID", "Estimate")]
LPI<-LPI[,c("ID","Binomial","Common_name","Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat")]

df<-merge(merge(temp,luc, by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

nrow(df)

df2<-subset(df, !is.na(Estimate)&r_sq >= 0.5  & !is.na(LUC_dist)&length_time >=5 & System!="Marine" &Specific_location == 1 )

nrow(df2)

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

library(lme4)
library(MuMIn)
source("rsquaredglmm.R")

R=999
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


#dt<-data.table(Euro)

land<-seq(min(dt$change_rate_scale),max(dt$change_rate_scale),length.out=1000)
temp<-seq(min(dt$mean_slope_scale),max(dt$mean_slope_scale),length.out=1000)

df<-expand.grid(landus,tempus)
colnames(df)<-c("change_rate_scale", "mean_slope_scale")

#unscaling the data
landus<-(land  * scale_luc) + centre_luc
tempus<-(temp * scale_temp) + centre_temp

dfus<-expand.grid(landus,tempus)
colnames(dfus)<-c("change_rate_scale", "mean_slope_scale")


for (i in 1:R) {
  dt2<-data.frame(dt[, ID[sample.int(.N, 1, TRUE)], by = loc_id])     #.N     signifies the number of rows when using data.table
  colnames(dt2)[2]<-"ID"
  sp_dups_df2<-sp_df_scale[sp_df_scale$ID %in% dt2$ID,]
  
  m1<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1T<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_dups_df2)
  
  m1a<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1aT<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=sp_dups_df2)
  
  m1b<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1bT<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=sp_dups_df2)
  
  m1c<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1cT<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=sp_dups_df2)
  
  mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_dups_df2, REML=F)
  mnullT<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_dups_df2)
  
  #AIC
  # AIC_m1[i]<-AIC(m1)
  # AIC_m1a[i]<-AIC(m1a)
  # AIC_m1b[i]<-AIC(m1b)
  # AIC_m1c[i]<-AIC(m1c)
  # AIC_mnull[i]<-AIC(mnull)
  # 
  # #Weights
  # 
  # msAICc <- model.sel(m1,m1a,m1b,m1c,mnull)
  # msAICc$model<-rownames(msAICc)
  # msAICc<-data.frame(msAICc)
  # 
  # m1_w[i]<-subset(msAICc, model=="m1")$weight
  # m1a_w[i]<-subset(msAICc, model=="m1a")$weight
  # m1b_w[i]<-subset(msAICc, model=="m1b")$weight
  # m1c_w[i]<-subset(msAICc, model=="m1c")$weight
  # mnull_w[i]<-subset(msAICc, model=="mnull")$weight
  # 
  #Rsq
  models_list<-list(m1,m1a,m1b,m1c,mnull)
  # modelsR<-lapply(models_list,rsquared.glmm)
  # modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)
  # 
  # marg_Rsq_m1[i]<-modelsRsq[1,4]
  # marg_Rsq_m1a[i]<-modelsRsq[2,4]
  # marg_Rsq_m1b[i]<-modelsRsq[3,4]
  # marg_Rsq_m1c[i]<-modelsRsq[4,4]
  # cond_Rsq_m1[i]<-modelsRsq[1,5]
  # cond_Rsq_m1a[i]<-modelsRsq[2,5]
  # cond_Rsq_m1b[i]<-modelsRsq[3,5]
  # cond_Rsq_m1c[i]<-modelsRsq[4,5]
  # cond_Rsq_mnull[i]<-modelsRsq[5,5]
  # 
  # var_imp<-summary(model.avg(models_list))
  # MTC_i[i]<-var_imp$importance["mean_slope_scale"]
  # LUC_i[i]<-var_imp$importance["change_rate_scale"]
  # LUC_MTC_i[i]<-var_imp$importance["change_rate_scale:mean_slope_scale"]
  # #BM_i[i]<-var_imp$importance["Bodymass"]
  # 
  # int_av[i]<-var_imp$coefmat.subset["(Intercept)","Estimate"]
  # MTC_av[i]<-var_imp$coefmat.subset["mean_slope_scale","Estimate"]
  # LUC_av[i]<-var_imp$coefmat.subset["change_rate_scale","Estimate"]
  # LUC_MTC_av[i]<-var_imp$coefmat.subset["change_rate_scale:mean_slope_scale","Estimate"]
  # 
  mav<-model.avg(models_list)
  pred<-predict(mav, dfus, re.form=NA)
  pdf<-data.frame(dfus, pred)
  pred2<-matrix(pred, ncol=length(landus), nrow=length(tempus), byrow=T)
  #image(land,temp,pred2)
  # head(pred2)
  # pras<-raster(pred2, xmn=min(landus), xmx=max(landus), ymn=min(tempus), ymx=max(tempus))
  # 
  pred_pcnt<-(10^pred2) - 1
  pcntras<-raster(pred_pcnt, xmn=min(landus), xmx=max(landus), ymn=min(tempus), ymx=max(tempus))
  
  file<-paste("predict_", i,".tif" ,sep="")
  writeRaster(pcntras, filename=file, overwrite=TRUE)
  print(i)
}


list<-list.files(path = getwd(), pattern = "predict.*\\.tif$")

br<-stack(list)

#plot(br)

br_av<-mean(br)

plot(br_av, xlab="Land Use Change Distance", ylab="Mean Temperature Change", main="Percentage Population Change")


# plot(pras)
# head(pred2)
# lev<-levelplot(pred ~ change_rate_scale + mean_slope_scale, data=pdf, useRaster = TRUE)
# 