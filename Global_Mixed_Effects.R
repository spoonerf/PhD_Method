temp<-read.csv("All_LPI_Mean_Temp_Slope.csv")
luc<-read.csv("LUC_distance_all.csv")
LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")
#pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_17.csv")
pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_18.csv")

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

luc_est = numeric(R)
clm_est = numeric(R)
int_est = numeric(R)

dt<-data.table(Euro)

for (i in 1:R) {
  dt2<-data.frame(dt[, ID[sample.int(.N, 1, TRUE)], by = loc_id])     #.N     signifies the number of rows when using data.table
  colnames(dt2)[2]<-"ID"
  sp_dups_df2<-Euro[Euro$ID %in% dt2$ID,]
  
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
  
  #Weights
  msAICc <- model.sel(m1,m1a,m1b,m1c,mnull)
  msAICc$model<-rownames(msAICc)
  msAICc<-data.frame(msAICc)
  
  m1_w[i]<-subset(msAICc, model=="m1")$weight
  m1a_w[i]<-subset(msAICc, model=="m1a")$weight
  m1b_w[i]<-subset(msAICc, model=="m1b")$weight
  m1c_w[i]<-subset(msAICc, model=="m1c")$weight
  mnull_w[i]<-subset(msAICc, model=="mnull")$weight
  
  #estimates from most complex model
  
  luc_est[i]<-fixef(m1)["change_rate_scale"]
  clm_est[i]<-fixef(m1)["mean_slope_scale"]
  int_est[i]<-fixef(m1)["change_rate_scale:mean_slope_scale"]
  print(i)

  }


AIC_df<-data.frame(cbind(AIC_m1,AIC_m1a, AIC_m1b, AIC_m1c, AIC_mnull))

AIC_del<-AIC_df[,c(1:4)] - AIC_df$AIC_mnull

colMeans(AIC_del)


hist(luc_est)
hist(clm_est)
hist(int_est)

luc_est_glb<-luc_est
clm_est_glb<-clm_est
int_est_glb<-int_est

luc_est_hil<-luc_est
clm_est_hil<-clm_est
int_est_hil<-int_est

mean(marg_Rsq_m1)
mean(marg_Rsq_m1a)
mean(marg_Rsq_m1b)
mean(marg_Rsq_m1c)

mean(cond_Rsq_m1)
mean(cond_Rsq_m1a)
mean(cond_Rsq_m1b)
mean(cond_Rsq_m1c)
mean(cond_Rsq_mnull)

mean(m1_w)
mean(m1a_w)
mean(m1b_w)
mean(m1c_w)
mean(mnull_w)

Euro<-sp_df_scale[sp_df_scale$ID %in% Inc5MammalsTemp$ID,] 

mint<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=Euro, REML=F)
madd<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=Euro, REML=F)
mlu<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=Euro, REML=F)
mcl<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=Euro, REML=F)
mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=Euro, REML=F)

AIC(mint)
AIC(madd)
AIC(mlu)
AIC(mcl)
AIC(mnull)

sp_dups<-data.frame(ddply(Inc5MammalsTemp,.(Longitude,Latitude),nrow))
sp_dups$loc_id<-1:length(sp_dups$Longitude)
sp_dups_df<-merge(sp_dups, Inc5MammalsTemp, by=c("Longitude","Latitude"))

####counting the number of mammals or birds at each location

result<-data.frame()

for (i in 1:length(sp_dups$loc_id)){
  
  sub<-subset(sp_dups_df, loc_id==i)
  bird_loc<-sum(sub$Class == "Aves")
  mamm_loc<-sum(sub$Class == "Mammalia")
  amph_loc<-sum(sub$Class == "Amphibia")
  count_loc<-cbind(i,bird_loc, mamm_loc,amph_loc)
  print(count_loc)
  result<-rbind(count_loc,result)
}

colnames(result)[1]<-"loc_id"

sp_dups_df2<-merge(sp_dups_df, result, by="loc_id")

###bootstrapping

library(data.table)
dt = as.data.table(sp_dups_df)

parm_df<-sp_dups_df[,c("ID","MnSlope50k", "mean_res_quart",  "change_rate_49")]  ##ID, land use, and climate

parm_mat<-as.matrix(parm_df)
parm_scale<-scale(parm_mat[,c("MnSlope50k", "mean_res_quart", "change_rate_49")])       #use the scaling factors at the bottom of these to scale the rasters

parm_id<-parm_mat[,"ID"]

parm_df_scale_inc<-data.frame(parm_id,parm_scale)

colnames(parm_df_scale_inc)<-c("ID","mean_slope_scale", "mean_var_scale", "change_rate_scale")

sp_df_scale_inc<-merge(sp_dups_df, parm_df_scale_inc, by="ID")

Euro_Inc<-sp_df_scale_inc[sp_df_scale_inc$ID %in% sp_df_scale$ID,] 

minti<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=Euro_Inc, REML=F)
maddi<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=Euro_Inc, REML=F)
mlui<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=Euro_Inc, REML=F)
mcli<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=Euro_Inc, REML=F)
mnulli<-lmer(lambda_sum ~ 1+(1|Binomial),data=Euro_Inc, REML=F)

AIC(minti)
AIC(maddi)
AIC(mlui)
AIC(mcli)
AIC(mnulli)




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

modelsRsq


mintb<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_df_scaleb, REML=F)
mnullb<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_df_scaleb, REML=F)

AIC(mintb)
AIC(mnullb)

mintm<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=sp_df_scalem, REML=F)
mnullm<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_df_scalem, REML=F)

AIC(mintm)
AIC(mnullm)


#Africa
Europe<-subset(sp_df_scale, Region=="Europe")

mint<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=Europe, REML=F)
madd<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=Europe, REML=F)
mlu<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=Europe, REML=F)
mcl<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=Europe, REML=F)
mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=Europe, REML=F)

AIC(mint)
AIC(madd)
AIC(mlu)
AIC(mcl)
AIC(mnull)

#######

Latin<-subset(sp_df_scale, Region=="Latin America and Caribbean")

mint<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=Latin, REML=F)
madd<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=Latin, REML=F)
mlu<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=Latin, REML=F)
mcl<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=Latin, REML=F)
mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=Latin, REML=F)

AIC(mint)
AIC(madd)
AIC(mlu)
AIC(mcl)
AIC(mnull)


#######

Africa<-subset(sp_df_scale, Region=="Africa")

mint<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=Africa, REML=F)
madd<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=Africa, REML=F)
mlu<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=Africa, REML=F)
mcl<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=Africa, REML=F)
mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=Africa, REML=F)

AIC(mint)
AIC(madd)
AIC(mlu)
AIC(mcl)
AIC(mnull)

########

NAm<-subset(sp_df_scale, Region=="North America")

mint<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=NAm, REML=F)
madd<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=NAm, REML=F)
mlu<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=NAm, REML=F)
mcl<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=NAm, REML=F)
mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=NAm, REML=F)

AIC(mint)
AIC(madd)
AIC(mlu)
AIC(mcl)
AIC(mnull)

########

Asia<-subset(sp_df_scale, Region=="Asia")

mint<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=Asia, REML=F)
madd<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=Asia, REML=F)
mlu<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=Asia, REML=F)
mcl<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=Asia, REML=F)
mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=Asia, REML=F)

AIC(mint)
AIC(madd)
AIC(mlu)
AIC(mcl)
AIC(mnull)

########
#########
Oce<-subset(sp_df_scale, Region=="Oceania")

mint<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial),data=Oce, REML=F)
madd<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial),data=Oce, REML=F)
mlu<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial),data=Oce, REML=F)
mcl<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=Oce, REML=F)
mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=Oce, REML=F)

AIC(mint)
AIC(madd)
AIC(mlu)
AIC(mcl)
AIC(mnull)

