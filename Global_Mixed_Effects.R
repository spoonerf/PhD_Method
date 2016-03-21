temp<-read.csv("All_LPI_Mean_Temp_Slope.csv")
luc<-read.csv("LUC_distance_all.csv")
LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")

#pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_17.csv")
pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_18.csv")
EurHil<-read.csv("Europe_HILDA_5_year_pops.csv")

head(pop)
temp<-temp[,c("ID", "Estimate")]
LPI<-LPI[,c("ID","Binomial","Common_name","Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat")]

df<-merge(merge(temp,luc, by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

nrow(df)

hab<-subset(df, Primary_threat == "Habitat degradation/change" | Primary_threat =="Habitat loss" | Secondary_threat =="Habitat degradation/change"| Secondary_threat =="Habitat loss"| Tertiary_threat =="Habitat degradation/change"| Tertiary_threat =="Habitat loss" )
nohab<-df[!(df$ID %in% hab$ID),]

hab$Habitat_Loss_Threat <- 1
nohab$Habitat_Loss_Threat<- 0 

df<-rbind(hab,nohab)

df2<-subset(df, !is.na(Estimate)&r_sq >= 0.5  & !is.na(LUC_dist)&length_time >=5 & System!="Marine" &Specific_location == 1 )

#EurGlo<-subset(df, !is.na(Estimate)&r_sq >= 0.5  & !is.na(LUC_dist)&length_time >=5 & System!="Marine" &Specific_location == 1)

#df2<-EurGlo[EurGlo$ID %in% EurHil$ID, ]  #Global data

# nrow(df2)
# 
# hab<-subset(df2, Primary_threat == "Habitat degradation/change" | Primary_threat =="Habitat loss" | Secondary_threat =="Habitat degradation/change"| Secondary_threat =="Habitat loss"| Tertiary_threat =="Habitat degradation/change"| Tertiary_threat =="Habitat loss" )
# 
# nohab<-df2[!(df2$ID %in% hab$ID),]
# 
# nrow(hab)
# nrow(nohab)
# 
# EUdf<-merge(EurHil, df, by="ID")
# 
# habE<-subset(EUdf, Primary_threat == "Habitat degradation/change" | Primary_threat =="Habitat loss" | Secondary_threat =="Habitat degradation/change"| Secondary_threat =="Habitat loss"| Tertiary_threat =="Habitat degradation/change"| Tertiary_threat =="Habitat loss" )
# 
# nohabE<-EUdf[!(EUdf$ID %in% habE$ID),]
# 
# nrow(habE)
# nrow(nohabE)

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

parm_id<-parm_mat[,"ID"]

parm_df_scale<-data.frame(parm_id,parm_scale)

colnames(parm_df_scale)<-c("ID","mean_slope_scale", "change_rate_scale")

sp_df_scale<-merge(sp_dups_df, parm_df_scale, by="ID")

dt<-data.table(sp_df_scale)


#sp_df_scale<-EurHil

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

for (i in 1:R) {
  dt2<-data.frame(dt[, ID[sample.int(.N, 1, TRUE)], by = loc_id])     #.N     signifies the number of rows when using data.table
  colnames(dt2)[2]<-"ID"
  sp_dups_df2<-sp_df_scale[sp_df_scale$ID %in% dt2$ID,]
  
  m1<-lmer(lambda_sum ~ Habitat_Loss_Threat+mean_slope_scale+Habitat_Loss_Threat:mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1T<-lmer(lambda_sum ~ Habitat_Loss_Threat+mean_slope_scale+Habitat_Loss_Threat:mean_slope_scale+(1|Binomial),data=sp_dups_df2)
  
  m1a<-lmer(lambda_sum ~ Habitat_Loss_Threat+mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1aT<-lmer(lambda_sum ~ Habitat_Loss_Threat+mean_slope_scale+(1|Binomial),data=sp_dups_df2)
  
  m1b<-lmer(lambda_sum ~ Habitat_Loss_Threat+(1|Binomial),data=sp_dups_df2, REML=F)
  m1bT<-lmer(lambda_sum ~ Habitat_Loss_Threat+(1|Binomial),data=sp_dups_df2)
  
  m1c<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=sp_dups_df2, REML=F)
  m1cT<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial),data=sp_dups_df2)
 
  mnull<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_dups_df2, REML=F)
  mnullT<-lmer(lambda_sum ~ 1+(1|Binomial),data=sp_dups_df2)
  
  #AIC
  AIC_m1[i]<-AIC(m1)
  AIC_m1a[i]<-AIC(m1a)
  AIC_m1b[i]<-AIC(m1b)
  AIC_m1c[i]<-AIC(m1c)
  AIC_mnull[i]<-AIC(mnull)

  #Weights
  
  msAICc <- model.sel(m1,m1a,m1b,m1c,mnull)
  msAICc$model<-rownames(msAICc)
  msAICc<-data.frame(msAICc)
  
  m1_w[i]<-subset(msAICc, model=="m1")$weight
  m1a_w[i]<-subset(msAICc, model=="m1a")$weight
  m1b_w[i]<-subset(msAICc, model=="m1b")$weight
  m1c_w[i]<-subset(msAICc, model=="m1c")$weight
  mnull_w[i]<-subset(msAICc, model=="mnull")$weight
  
  #Rsq
  models_list<-list(m1T,m1aT,m1bT,m1cT,mnullT)
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

  
  #estimates from most complex model
  var_imp<-summary(model.avg(models_list))
  MTC_i[i]<-var_imp$importance["mean_slope_scale"]
  LUC_i[i]<-var_imp$importance["Habitat_Loss_Threat"]
  LUC_MTC_i[i]<-var_imp$importance["Habitat_Loss_Threat:mean_slope_scale"]
  #BM_i[i]<-var_imp$importance["Bodymass"]

  int_av[i]<-var_imp$coefmat.subset["(Intercept)","Estimate"]
  MTC_av[i]<-var_imp$coefmat.subset["mean_slope_scale","Estimate"]
  LUC_av[i]<-var_imp$coefmat.subset["Habitat_Loss_Threat","Estimate"]
  LUC_MTC_av[i]<-var_imp$coefmat.subset["Habitat_Loss_Threat:mean_slope_scale","Estimate"]
  
  print(i)
  }


AIC_df<-data.frame(cbind(AIC_m1,AIC_m1a, AIC_m1b, AIC_m1c, AIC_mnull))

AIC_del<-AIC_df[,c(1:4)] - AIC_df$AIC_mnull

colMeans(AIC_del)


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

mean_av<- c(mean(LUC_av),mean(MTC_av),mean(LUC_MTC_av))
lowCI_av<-c(sort(LUC_av)[25], sort(MTC_av)[25], sort(LUC_MTC_av)[25])
highCI_av<-c(sort(LUC_av)[975], sort(MTC_av)[975], sort(LUC_MTC_av)[975])

Variable<-c("LUC", "MTC", "LUC*MTC")

conf_av<-data.frame(rbind( lowCI_av, mean_av, highCI_av))
colnames(conf_av)<-Variable
conf_av

library(plotrix)

plotCI(1:3, mean_av, (highCI_av-mean_av), (mean_av-lowCI_av), ylab="Coefficient (95% C.I.)", xlab="" ,xaxt = "n", 
       main="Variable Coefficients", lwd=1, ylim=c(min(lowCI_av*1.1), max(highCI_av*1.1)))
axis(1, at=1:3, labels=colnames(conf_i), las=2)
abline(h=0, col="red", lty =2)



#############
Euro<-Inc5MammalsTemp[!(Inc5MammalsTemp$ID %in% sp_df_scale$ID),] 
Euro<-subset(sp_df_scale, Region=="Europe")
Afr<-subset(sp_df_scale, Region=="Africa")
LAm<-subset(sp_df_scale, Region=="Latin America and Caribbean")
Asia<-subset(sp_df_scale, Region=="Asia")
NAm<-subset(sp_df_scale, Region=="North America")
Oce<-subset(sp_df_scale, Region=="Oceania")

mintE<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial)+(1|Country),data=Euro, REML=F)
maddE<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial)+(1|Country),data=Euro, REML=F)
mluE<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial)+(1|Country),data=Euro, REML=F)
mchE<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial)+(1|Country),data=Euro, REML=F)
mnullE<-lmer(lambda_sum ~ 1+(1|Binomial)+(1|Country),data=Euro, REML=F)

AIC(mintE)
AIC(maddE)
AIC(mluE)
AIC(mchE)
AIC(mnullE)

mintA<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial)+(1|Country),data=Afr, REML=F)
mnullA<-lmer(lambda_sum ~ 1+(1|Binomial)+(1|Country),data=Afr, REML=F)


AIC(mintA)
AIC(mnullA)


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

