
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

nrow(df)
#nrow(dfa)
nrow(dfd) 


df2<-subset(dfd, !is.na(Estimate)  & r_sq >= 0.4999999 &length_time >=5& System!="Marine"
            &Specific_location == 1 &!is.na(both_change) & !is.na(Log_Body_Mass_g)
            & ( Class=="Aves") & Protected_status != "Unknown"  & Protected_status != "Both")
nrow(df2)

# df2<-subset(dfd, !is.na(Estimate) &length_time >=5& System!="Marine"
#             &Specific_location == 1 &!is.na(both_change) & !is.na(Log_Body_Mass_g)
#             & ( Class=="Mammalia") & Protected_status != "Unknown"  & Protected_status != "Both")
# nrow(df2)


df2$Protected_status[df2$Protected_status == "No (area surrounding PA)"] <- "No"
df2$Protected_status[df2$Protected_status == "No (large survey area)"] <- "No"

table(df2$Class)


df2[is.na(df2$lambda_mean),]$lambda_mean<-0

# # #mammals
# diet<-read.csv("mammaldiet.csv")
# diet$Binomial<-paste(diet$Genus, diet$Species, sep= "_")
# df2<-merge(df2, diet[,c(6,24,31)], by="Binomial")
# 
# nrow(df2)
# library(dplyr)
# #birds 
# diet<-read.csv("birddiet.csv")
# diet$Binomial<-gsub(" ", "_", diet$Scientific)
# df2<-merge(df2, diet[,c(20,41)], by="Binomial")

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

#write.csv(dt, "GCB_Data.csv")

#removing atlantic forest populations

#dt<-dt[dt$loc_id != 109 ,]


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

  mav_aic6<-model.avg(msAICc, subset =  delta < 6)
  
  #mav<-mav_aic6
  
  smav<-summary(mav)
  
  coef_av<-smav$coefmat.subset[,"Estimate"]
  coef_df<-data.frame(coef_av)
  coef_df$lowCI<-confint(mav)[,1]
  coef_df$highCI<-confint(mav)[,2]
  coef_df
  
  
  coef_pcnt<-data.frame(((10^coef_df) - 1)*100)
  coef_pcnt

  

  coef_pcnt$Var_name<-rownames(coef_pcnt)

library(coefplot)
library(ggplot2)

coef_pcnt$val<-1:6
coef_pcnt$var_name <- factor(coef_pcnt$Var_name, levels = coef_pcnt$Var_name[order(coef_pcnt$val)])
  
  p1<-ggplot(coef_pcnt)
  p1<- p1 + geom_hline(yintercept = 0, colour=gray(1/2), lty=2)
  p1<- p1 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=1.5, position = position_dodge(width=1/2))
  p1<- p1 + geom_pointrange(aes(x= Var_name, y=coef_av, ymin=lowCI, ymax=highCI), lwd=1, position=position_dodge(width=1/2), shape=21, fill="White")
  p1<- p1 + scale_y_continuous(breaks=seq(-8, 4, 2)) +theme_bw() + labs(y = "Annual Population Change (%)", x = "Variable") + theme(legend.title=element_blank(), text = element_text(size=20),axis.title.x = element_text(margin = unit(c(5, 5, 0, 0), "mm")))
  print(p1)
  
  
  #coef_pcnt$Class<-"Birds"  
  coef_pcnt$Class<-"Mammals"  
  
  #coef_pcntb<-coef_pcnt
  coef_pcntm<-coef_pcnt  
  
  
  
  smav6<-summary(mav_aic6)
  
  coef_av6<-smav6$coefmat.subset[,"Estimate"]
  coef_df6<-data.frame(coef_av6)
  coef_df6$lowCI<-confint(mav_aic6)[,1]
  coef_df6$highCI<-confint(mav_aic6)[,2]
  coef_df6
  
  
  coef_pcnt6<-data.frame(((10^coef_df6) - 1)*100)
  coef_pcnt6
  
  coef_pcnt6$Var_name<-rownames(coef_pcnt6)
  
  library(coefplot)
  library(ggplot2)
  
  coef_pcnt6$val<-1:6
  coef_pcnt6$var_name <- factor(coef_pcnt6$Var_name, levels = coef_pcnt6$Var_name[order(coef_pcnt6$val)])
  
  p16<-ggplot(coef_pcnt6)
  p16<- p16 + geom_hline(yintercept = 0, colour=gray(1/2), lty=2)
  p16<- p16 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=1.5, position = position_dodge(width=1/2))
  p16<- p16 + geom_pointrange(aes(x= Var_name, y=coef_av6, ymin=lowCI, ymax=highCI), lwd=1, position=position_dodge(width=1/2), shape=21, fill="White")
  p16<- p16 + scale_y_continuous(breaks=seq(-8, 4, 2)) +theme_bw() + labs(y = "Annual Population Change (%)", x = "Variable") + theme(legend.title=element_blank(), text = element_text(size=20),axis.title.x = element_text(margin = unit(c(5, 5, 0, 0), "mm")))
  print(p16)
  
  
  #coef_pcnt$Class<-"Birds"  
  coef_pcnt6$Class<-"Mammals"  
  
  #coef_pcntb<-coef_pcnt
  coef_pcntm6<-coef_pcnt6
  
  

#####Birds


df2<-subset(dfd, !is.na(Estimate)  & r_sq >= 0.4999999 &length_time >=5& System!="Marine" 
            &Specific_location == 1 &!is.na(both_change) & !is.na(Log_Body_Mass_g)
            & ( Class=="Aves") & Protected_status != "Unknown"  & Protected_status != "Both")

nrow(df2)

df2<-subset(dfd, !is.na(Estimate)  &length_time >=5& System!="Marine" 
            &Specific_location == 1 &!is.na(both_change) & !is.na(Log_Body_Mass_g)
            & ( Class=="Aves") & Protected_status != "Unknown"  & Protected_status != "Both")

nrow(df2)

df2$Protected_status[df2$Protected_status == "No (area surrounding PA)"] <- "No"
df2$Protected_status[df2$Protected_status == "No (large survey area)"] <- "No"

table(df2$Class)

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

#write.csv(dt, "GCB_Data.csv")

#removing atlantic forest populations

#dt<-dt[dt$loc_id != 109 ,]


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
#m1cr<-lmer(lambda_mean ~ mean_slope_scale+WWF_REALM2+(1|Binomial)+(1|loc_id),data=dt, REML=F)
#m1ce<-lmer(lambda_mean ~ mean_slope_scale+Elevation + (1|Binomial)+(1|loc_id),data=dt, REML=F)
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

AIC(m0,m0a,m0b,m0c,m1,m1a,m1b,m1c,mnull)
AIC(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull,m0f,m0af,m0bf,m0cf,m0df,m1f,m1af,m1bf,m1cf)

#Rsq
models_list<-list(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull)
models_list<-list(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull,m0f,m0af,m0bf,m0cf,m0df,m1f,m1af,m1bf,m1cf)

modelsR<-lapply(models_list,rsquared.glmm)
modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)
rownames(modelsRsq)<-c("m0","m0a","m0b","m0c","m0d","m1","m1a","m1b","m1c","mnull","m0f","m0af","m0bf","m0cf","m0df","m1f","m1af","m1bf","m1cf")
modelsRsq

library(MuMIn)
var_imp<-summary(model.avg(msAICc))


mav<-model.avg(msAICc, subset =  cumsum(weight) <= .95)

mav_aic6<-model.avg(msAICc, subset =  delta < 6)

#mav<-mav_aic6

smav<-summary(mav)

coef_av<-smav$coefmat.subset[,"Estimate"]
coef_df<-data.frame(coef_av)
coef_df$lowCI<-confint(mav)[,1]
coef_df$highCI<-confint(mav)[,2]
coef_df


coef_pcnt<-data.frame(((10^coef_df) - 1)*100)
coef_pcnt


coef_pcnt$Var_name<-rownames(coef_pcnt)


library(coefplot)
library(ggplot2)

coef_pcnt$val<-1:6
coef_pcnt$var_name <- factor(coef_pcnt$Var_name, levels = coef_pcnt$Var_name[order(coef_pcnt$val)])

p1<-ggplot(coef_pcnt)
p1<- p1 + geom_hline(yintercept = 0, colour=gray(1/2), lty=2)
p1<- p1 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=1.5, position = position_dodge(width=1/2))
p1<- p1 + geom_pointrange(aes(x= Var_name, y=coef_av, ymin=lowCI, ymax=highCI), lwd=1, position=position_dodge(width=1/2), shape=21, fill="White")
p1<- p1 + scale_y_continuous(breaks=seq(-8, 4, 2)) +theme_bw() + labs(y = "Annual Population Change (%)", x = "Variable") + theme(legend.title=element_blank(), text = element_text(size=20),axis.title.x = element_text(margin = unit(c(5, 5, 0, 0), "mm")))
print(p1)


coef_pcnt$Class<-"Birds"  
#coef_pcnt$Class<-"Mammals"  

coef_pcntb<-coef_pcnt
#coef_pcntm<-coef_pcnt  

smav6<-summary(mav_aic6)

coef_av6<-smav6$coefmat.subset[,"Estimate"]
coef_df6<-data.frame(coef_av6)
coef_df6$lowCI<-confint(mav_aic6)[,1]
coef_df6$highCI<-confint(mav_aic6)[,2]
coef_df6


coef_pcnt6<-data.frame(((10^coef_df6) - 1)*100)
coef_pcnt6

coef_pcnt6$Var_name<-rownames(coef_pcnt6)

library(coefplot)
library(ggplot2)

coef_pcnt6$val<-1:6
coef_pcnt6$var_name <- factor(coef_pcnt6$Var_name, levels = coef_pcnt6$Var_name[order(coef_pcnt6$val)])

p16<-ggplot(coef_pcnt6)
p16<- p16 + geom_hline(yintercept = 0, colour=gray(1/2), lty=2)
p16<- p16 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=1.5, position = position_dodge(width=1/2))
p16<- p16 + geom_pointrange(aes(x= Var_name, y=coef_av6, ymin=lowCI, ymax=highCI), lwd=1, position=position_dodge(width=1/2), shape=21, fill="White")
p16<- p16 + scale_y_continuous(breaks=seq(-8, 4, 2)) +theme_bw() + labs(y = "Annual Population Change (%)", x = "Variable") + theme(legend.title=element_blank(), text = element_text(size=20),axis.title.x = element_text(margin = unit(c(5, 5, 0, 0), "mm")))
print(p16)


coef_pcnt6$Class<-"Birds"  
#coef_pcnt$Class<-"Mammals"  

coef_pcntb6<-coef_pcnt6
#coef_pcntm<-coef_pcnt


coef_both<-rbind(coef_pcntb[,c(1,2,3,4,7)], coef_pcntm[,c(1,2,3,4,7)])

coef_both6<-rbind(coef_pcntb6[,c(1,2,3,4,7)], coef_pcntm6[,c(1,2,3,4,7)])

coef_both$Var_name
coef_both$Var_name[coef_both$Var_name == "(Intercept)"] <- "aIntercept"
coef_both$Var_name[coef_both$Var_name == "mean_slope_scale"] <- "bMTC"
coef_both$Var_name[coef_both$Var_name == "change_rate_scale"] <- "cLUC"
coef_both$Var_name[coef_both$Var_name == "change_rate_scale:mean_slope_scale"] <- "dMTC*LUC"
coef_both$Var_name[coef_both$Var_name == "Bodymass_scale"] <- "eBM"
coef_both$Var_name[coef_both$Var_name == "Protected_statues"] <- "fPA"


coef_both6$Var_name
coef_both6$Var_name[coef_both6$Var_name == "(Intercept)"] <- "aIntercept"
coef_both6$Var_name[coef_both6$Var_name == "mean_slope_scale"] <- "bMTC"
coef_both6$Var_name[coef_both6$Var_name == "change_rate_scale"] <- "cLUC"
coef_both6$Var_name[coef_both6$Var_name == "change_rate_scale:mean_slope_scale"] <- "dMTC*LUC"
coef_both6$Var_name[coef_both6$Var_name == "Bodymass_scale"] <- "eBM"
coef_both6$Var_name[coef_both6$Var_name == "Protected_statues"] <- "fPA"


#write.csv(coef_both, "Model_Average_coefs_No_Atlantic_Forest.csv")

#write.csv(coef_both, "Model_Average_coefs_aic6.csv")
#write.csv(coef_both, "Model_Average_coefs4.csv")
coef_both<-read.csv("Model_Average_coefs4.csv")
#coef_both<-read.csv("Model_Average_coefs_No_Atlantic_Forest.csv")
coef_both6<-read.csv("Model_Average_coefs_aic6.csv")


ele<-readPNG("elephant.png")
pel<-readPNG("pelican.png")
e<-rasterGrob(ele, interpolate = FALSE)
p<-rasterGrob(pel, interpolate = FALSE)

# coef_old<-coef_both

library(ggplot2)
p1<-ggplot(coef_both6, aes(colour=Class))
p1<- p1 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=2.5, position = position_dodge(width=2/3))
p1<- p1 + geom_pointrange(aes(x= Var_name, y=coef_av6, ymin=lowCI, ymax=highCI), lwd=2, position=position_dodge(width=2/3), shape=21, fill="White")
p1<- p1 + scale_y_continuous(breaks=seq(-10, 14, 4), limits=(c(-10,15)))
p1<-p1 + theme_bw() + labs(y = "Population Change (%)", x = "")
p1<- p1 + theme(legend.position="none",text=element_text(size=20),axis.text.x=element_text(size=20) , axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p1<- p1 + scale_color_manual(values=c("black", "black"))
p1<-p1 + theme_bw() 
p1<-p1+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p1<- p1+ geom_hline(yintercept = 0, linetype=2)
p1<- p1+ geom_vline(xintercept = 1.5, linetype=2)

p1<-p1+ scale_x_discrete(labels=c("aIntercept" ="Intercept", "bMTC" ="Rate of climate \nwarming \n(RCW)",  "cLUC" = "Rate of \nconversion to \nanthropogenic \nland use \n(RCA)", "dMTC*LUC"= "RCW:RCA", "eBM"= "Body mass", "Protected_statusYes" = "Inside protected \narea"))
p1<- p1+annotation_custom(e, xmin=0.2, xmax=2.1, ymin=3.33, ymax=5.23)+
annotation_custom(e, xmin=1.2, xmax=3.1, ymin=-0.305, ymax=1.605)+
annotation_custom(e, xmin=2.2, xmax=4.1, ymin= 1.025, ymax=2.925)+
annotation_custom(e, xmin=3.2, xmax=5.1, ymin= 0.327, ymax=2.227)+
annotation_custom(e, xmin=4.2, xmax=6.1, ymin= 3.078, ymax=4.978)+
annotation_custom(e, xmin=5.2, xmax=7.1, ymin= 5.072, ymax=6.907)+
annotation_custom(p, xmin=-0.115, xmax=1.735, ymin= 0.504, ymax=2.404)+
annotation_custom(p, xmin=0.885, xmax=2.735, ymin= -2.776, ymax=-0.876)+
annotation_custom(p, xmin=1.885, xmax=3.735, ymin= 2.672, ymax=4.572)+
annotation_custom(p, xmin=2.885, xmax=4.735, ymin= 1.445, ymax=3.345)+
annotation_custom(p, xmin=3.885, xmax=5.735, ymin= 2.756, ymax=4.656)+
annotation_custom(p, xmin=4.885, xmax=6.735, ymin= 13.213, ymax=15.113)
p1

png(filename="Figure3_1000.png",  units="in", width=12, height=8, pointsize=12, res=1000)

print(p1)

dev.off()


print(p1)



((10^msAICc[,1:5])-1)*100


coef_diff<-data.frame((coef_both$coef_av - coef_both6$coef_av6),(coef_both$lowCI - coef_both6$lowCI),  (coef_both$highCI - coef_both6$highCI), coef_both$Var_name, coef_both$Class)
colnames(coef_diff)<-colnames(coef_both)


library(ggplot2)
p1<-ggplot(coef_diff, aes(colour=Class))
p1<- p1 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=2.5, position = position_dodge(width=2/3))
p1<- p1 + geom_pointrange(aes(x= Var_name, y=coef_av, ymin=lowCI, ymax=highCI), lwd=2, position=position_dodge(width=2/3), shape=21, fill="White")
#p1<- p1 + scale_y_continuous(breaks=seq(-10, 14, 4), limits=(c(-10,15)))
p1<-p1 + theme_bw() + labs(y = "Population Change (%)", x = "")
p1<- p1 + theme(legend.position="none",text=element_text(size=20),axis.text.x=element_text(size=20) , axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p1<- p1 + scale_color_manual(values=c("black", "black"))
p1<-p1 + theme_bw() 
p1<-p1+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p1<- p1+ geom_hline(yintercept = 0, linetype=2)
p1<- p1+ geom_vline(xintercept = 1.5, linetype=2)
print(p1)



library(ggplot2)
dt_pred<-data.frame(dt, predict(m1c))


# ggplot(dt_pred, aes(x = Longitude, y=Latitude, colour =10^predict.m1c. ))+
#   geom_point()+scale_color_gradient2(midpoint=1, low="red",
#                                      high="green" )
# 
# 
# CR30s<-brick("cru_ts3.23.1931.1940.tmp.dat.nc")
# CR40s<-brick("cru_ts3.23.1941.1950.tmp.dat.nc")
# CR50s<-brick("cru_ts3.23.1951.1960.tmp.dat.nc")
# CR60s<-brick("cru_ts3.23.1961.1970.tmp.dat.nc")
# CR70s<-brick("cru_ts3.23.1971.1980.tmp.dat.nc")
# CR80s<-brick("cru_ts3.23.1981.1990.tmp.dat.nc")
# CR90s<-brick("cru_ts3.23.1991.2000.tmp.dat.nc")
# CR00s<-brick("cru_ts3.23.2001.2010.tmp.dat.nc")
# CR10s<-brick("cru_ts3.23.2011.2014.tmp.dat.nc")
# CR<-stack(CR30s,CR40s,CR50s,CR60s,CR70s,CR80s,CR90s,CR00s,CR10s)
# 
# X <- cbind(1, time)
# 
# ## pre-computing constant part of least squares
# invXtX <- solve(t(X) %*% X) %*% t(X)
# 
# ## much reduced regression model; [2] is to get the slope
# quickfun <- function(y) (invXtX %*% y)[2]
# x4 <- calc(CR, quickfun) 
# map.p <- rasterToPoints(x4)
# 
# #Make the points a dataframe for ggplot
# df <- data.frame(map.p)
# #Make appropriate column headings
# colnames(df) <- c("Longitude", "Latitude", "MAP")
# 
# ggplot(data=df, aes(y=Latitude, x=Longitude)) +
#   geom_raster(aes(fill=MAP)) +
#      geom_point(data=dt_pred, aes(x=Longitude, y=Latitude, colour =10^predict.m1c.))+
#   scale_color_gradient2(midpoint=1, low="red",high="green" )+
#   coord_equal()
# 
# 
# #######fitted and predicted values
# 
# #dt$fitted<-fitted(m1c)
# dt$fitted<-predict(m1c, re.form=NA)
# rm1c<-ranef(m1c)
# rm1cb<-rm1c$Binomial
# rm1cb$Binomial<-rownames(rm1cb)
# colnames(rm1cb)[1]<-"Binomial_ranef"
# dt<-merge(dt, rm1cb, by="Binomial")
# 
# rm1cl<-rm1c$loc_id
# rm1cl$loc_id<-rownames(rm1cl)
# colnames(rm1cl)
# colnames(rm1cl)[1]<-"loc_id_ranef"
# rm1cl$loc_id<-as.numeric(rm1cl$loc_id)
# dt<-merge(dt, rm1cl, by="loc_id")
# dt$predicted<-dt$fitted+dt$Binomial_ranef+dt$loc_id_ranef
# 
# pred_fit<-data.frame(dt$ID, dt$lambda_mean, dt$fitted, dt$predicted)
# 
# 
# plot(pred_fit$dt.predicted, pred_fit$dt.lambda_mean)
# points(pred_fit$dt.fitted, pred_fit$dt.lambda_mean, col="red")
# 
# pred_melt<-melt(pred_fit, id = c("dt.ID", "dt.lambda_mean" ))
# 
# library(ggplot2)
# 
# ggplot(pred_melt, aes(x = 10^value, y = 10^dt.lambda_mean, colour = variable))+
#   geom_point(size = 3,  alpha = 0.3 )+
#   geom_smooth(method = "lm", se=FALSE)+
#   labs( x = "Observed Population Growth Rate", y = "Predicted Growth Rate", color = "")+
#   scale_color_manual(labels = c("Climate Effects Only","Including Random Effects"), values = c("Red", "Black")) +
#   theme(axis.title=element_text(size=14))+ theme(legend.text=element_text(size=14))+
#   #coord_equal( xlim=c(0, 2.2), ylim=c(0, 2.2))+
#   geom_abline(linetype = "dotted",slope=1, intercept=0)
# 
# 
# 
# pred_fit_table<-data.frame(dt$Binomial, dt$Country, dt$lambda_mean, dt$Estimate, dt$both_change, dt$fitted, dt$predicted)
# colnames(pred_fit_table)<-c("Binomial", "Country", "Average Population Growth Rate", "Rate of Climate Change", "Rate of Land Use Change", "Growth Rate Predicted From Climate Change", "Growth Rate Predicted From Climate Change, Species and Location")
# pred_fit_table$`Average Population Growth Rate`<-10^pred_fit_table$`Average Population Growth Rate`
# pred_fit_table$`Growth Rate Predicted From Climate Change`<-10^pred_fit_table$`Growth Rate Predicted From Climate Change`
# pred_fit_table$`Growth Rate Predicted From Climate Change, Species and Location`<-10^pred_fit_table$`Growth Rate Predicted From Climate Change, Species and Location`
# 
# #write.csv(pred_fit_table, "Fitted_values_climate_mammals.csv")
# 
# 
# 
# 
# 
# sp_ran<-data.frame(ranef(m0)[2])
# LPI<-read.csv("LPI_pops_20160523_edited.csv")
# head(LPI)
# ord<-data.frame(LPI$Binomial, LPI$Family, LPI$Order)
# colnames(ord)<-c("Binomial", "Family", "Order")
# sp_ran$Binomial<-rownames(sp_ran)
# sp_ran<-data.frame(sp_ran)
# colnames(sp_ran)<-c("Intercept","Binomial")
# 
# unique(merge(ord, sp_ran, by="Binomial"))
# sp_ord<-unique(merge(ord, sp_ran, by="Binomial"))
# sp_ord<-sp_ord[sp_ord$Family != "",]
# plot(sp_ord$Intercept, col=sp_ord$Family)
# fam_int<-aggregate(sp_ord$Intercept,by=list(sp_ord$Family) ,FUN="mean")
# fam_int$pgr<-10^fam_int$x
# fam_int
# colnames(fam_int)<-c("Family", "Intercept", "Lambda")
# ord_int<-aggregate(sp_ord$Intercept,by=list(sp_ord$Order) ,FUN="mean")
# ord_int$pgr<-10^ord_int$x
# ord_int
# colnames(ord_int)<-c("Order", "Intercept", "Lambda")
# 
# boxplot(10^sp_ord$Intercept ~ sp_ord$Order)
# 
# which.min(fitted.values(m0))
# 
# dt[285,]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
