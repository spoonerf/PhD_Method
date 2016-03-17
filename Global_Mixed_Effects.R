temp<-read.csv("All_LPI_Mean_Temp_Slope.csv")
luc<-read.csv("LUC_distance_all.csv")
LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")
pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_16.csv")

head(pop)
temp<-temp[,c("ID", "Estimate")]
LPI<-LPI[,c("ID","Binomial","Common_name","Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude")]


df<-merge(merge(temp,luc, by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

df2<-subset(df, !is.na(Estimate)&r_sq >= 0.5  & !is.na(LUC_dist)&length_time >=5 & System!="Marine" &Specific_location == 1)

nrow(df2)


library(plyr)
#counting duplicates at each location
sp_dups<-data.frame(ddply(df2,.(Longitude,Latitude),nrow))
sp_dups$loc_id<-1:length(sp_dups$Longitude)
sp_dups_df<-merge(sp_dups, df2, by=c("Longitude","Latitude"))

library(data.table)
dt = as.data.table(sp_dups_df)

parm_df<-sp_dups_df[,c("ID","Estimate", "LUC_dist")]  ##ID, land use, and climate

parm_mat<-as.matrix(parm_df)
parm_scale<-scale(parm_mat[,c("Estimate", "LUC_dist")])       #use the scaling factors at the bottom of these to scale the rasters

parm_id<-parm_mat[,"ID"]

parm_df_scale<-data.frame(parm_id,parm_scale)

colnames(parm_df_scale)<-c("ID","mean_slope_scale", "change_rate_scale")

sp_df_scale<-merge(sp_dups_df, parm_df_scale, by="ID")


library(lme4)
library(MuMIn)
source("rsquaredglmm.R")

R=9999
AIC_m1= numeric(R)
AIC_m1a= numeric(R)
AIC_m1b= numeric(R)
AIC_m1c= numeric(R)
AIC_mnull= numeric(R)

marg_Rsq_m1 = numeric(R)
marg_Rsq_m1a = numeric(R)
marg_Rsq_m1b = numeric(R)
marg_Rsq_m1c = numeric(R)

m1_w = numeric(R)
m1a_w = numeric(R)
m1b_w = numeric(R)
m1c_w = numeric(R)
mnull_w = numeric(R)

for (i in 1:R) {
  dt2<-data.frame(dt[, ID[sample.int(.N, 1, TRUE)], by = loc_id])     #.N     signifies the number of rows when using data.table
  colnames(dt2)[2]<-"ID"
  sp_dups_df2<-sp_df_scale[sp_df_scale$ID %in% dt2$ID,]
  
  m1<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  
  m1a<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  
  m1b<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  
  m1c<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  
  mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_dups_df2, REML=F)
#AIC
  AIC_m1[i]<-AIC(m1)
  AIC_m1a[i]<-AIC(m1a)
  AIC_m1b[i]<-AIC(m1b)
  AIC_m1c[i]<-AIC(m1c)
  AIC_mnull[i]<-AIC(mnull)
#Rsq
  models_list<-list(m1,m1a,m1b,m1c)
  modelsR<-lapply(models_list,rsquared.glmm)
  modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)
  
  marg_Rsq_m1[i]<-modelsRsq[1,4]
  marg_Rsq_m1a[i]<-modelsRsq[2,4]
  marg_Rsq_m1b[i]<-modelsRsq[3,4]
  marg_Rsq_m1c[i]<-modelsRsq[4,4]
#Weights
  msAICc <- model.sel(m1,m1a,m1b,m1c,mnull)
  msAICc$model<-rownames(msAICc)
  msAICc<-data.frame(msAICc)
  
  m1_w[i]<-subset(msAICc, model=="m1")$weight
  m1a_w[i]<-subset(msAICc, model=="m1a")$weight
  m1b_w[i]<-subset(msAICc, model=="m1b")$weight
  m1c_w[i]<-subset(msAICc, model=="m1c")$weight
  mnull_w[i]<-subset(msAICc, model=="mnull")$weight
  
  print(i)
}


AIC_df<-data.frame(cbind(AIC_m1,AIC_m1a, AIC_m1b, AIC_m1c, AIC_mnull))

AIC_del<-AIC_df[,c(1:4)] - AIC_df$AIC_mnull


colMeans(AIC_del)

mean(marg_Rsq_m1)
mean(marg_Rsq_m1a)
mean(marg_Rsq_m1b)
mean(marg_Rsq_m1c)

mean(m1_w)
mean(m1a_w)
mean(m1b_w)
mean(m1c_w)
mean(mnull_w)

Euro<-sp_df_scale[sp_df_scale$ID %in% Inc5MammalsTemp$ID,] 

minte<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=Euro, REML=F)
mnulle<-lmer(lambda_sum ~ 1+(1|Binomial),data=Euro, REML=F)

AIC(minte)
AIC(mnulle)

sp_df_scaleb<-subset(sp_df_scale, Class=="Aves")
sp_df_scalem<-subset(sp_df_scale, Class=="Mammalia")

mint<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_df_scale, REML=F)
madd<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=sp_df_scale, REML=F)
mlu<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=sp_df_scale, REML=F)
mcl<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=sp_df_scale, REML=F)
mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_df_scale, REML=F)

AIC(mint)
AIC(madd)
AIC(mlu)
AIC(mcl)
AIC(mnull)

anova(mint, mnull)

models_list<-list(mint,madd,mlu,mcl,mnull)
modelsR<-lapply(models_list,rsquared.glmm)
modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)




mintb<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_df_scaleb, REML=F)
mnullb<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_df_scaleb, REML=F)

AIC(mintb)
AIC(mnullb)

mintm<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_df_scalem, REML=F)
mnullm<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_df_scalem, REML=F)

AIC(mintm)
AIC(mnullm)



