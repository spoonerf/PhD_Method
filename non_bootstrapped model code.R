
#temp<-read.csv("All_LPI_Mean_Temp_Slope.csv")
#temp<-read.csv("All_LPI_Mean_Temp_Slope_nobuff.csv")
temp<-read.csv("All_LPI_All_Years_Nobuff_1931_moreLPI_end2005.csv")
#temp<-read.csv("All_LPI_All_Years_Nobuff_1931_moreLPI_min_mean_temp.csv")
#temp<-read.csv("All_LPI_All_Years_Nobuff_1931.csv")
#temp<-read.csv("All_LPI_All_Years_Nobuff.csv")
body<-read.csv("bird_and_mammal_traits2.csv")
body2<-read.csv("Bird_and_Mammal_BM.csv")
# body2<-read.csv("LPI_traits.csv")
# body<-read.csv("LPI_traits.csv")
# body3<-rbind(body, body2)
# 
# 
# body4<-subset(body3, !duplicated(ID))


body3<-rbind(body, body2)

hyde<-read.csv("Hyde_crop_pasture_annual_change.csv")
###
#hyde<-read.csv("Sensitivity_Land_Use_morepop.csv")

#LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")
LPI<-read.csv("LPI_pops_20160523_edited.csv")

#Realm<-read.csv("selected_pops_Ecoregion.csv")
Realm<-read.csv("Realm.csv", na.strings="")
Realm<-Realm[,c("ID", "WWF_REALM2")]
Realm2<-read.csv("hub_realms.csv")
colnames(Realm2)<-c("ID", "WWF_REALM2")

Realm3<-rbind(Realm,Realm2) #sort this has duplicates
Realm<-Realm[!is.na(Realm3$WWF_REALM2),]

#EurHil<-read.csv("Europe_HILDA_5_year_pops.csv")  # data from Euro-centric analysis
pop<-read.csv("Global_Population_Trends_Rsq_Lambda_07_10_16.csv")

temp<-temp[,c("ID", "Estimate")]

LPI<-LPI[,c("ID","Binomial","Common_name", "Order", "Protected_status", "Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat", "Migratory", "Forest")]

df<-merge(merge(temp,Realm, by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

dfb<-merge(df, body3[,c(3:5)], by="ID", all=TRUE)     #41 pops bodysizes missing for birds

#dfc<-merge(dfb, forest, by="ID", all=TRUE)
dfd<-merge(dfb, hyde[,c(-1,-3)], by="ID")


nrow(df)
#nrow(dfa)
nrow(dfd) 

df2<-subset(dfd, !is.na(Estimate) & r_sq >= 0.4999999  &length_time >=5 & System!="Marine" 
            &Specific_location == 1 &!is.na(both_change) & !is.na(Bodymass_g))

#df2<-subset(dfd, !is.na(Estimate) & r_sq >= 0.4999999  &!is.na(change_rate_49)&length_time >=5 & 
#             System!="Marine" &Specific_location == 1 & !is.na(Bodymass)&Class=="Aves")#& Class=="Mammalia"&((Primary_threat =="Habitat 

#df2bm<-subset(dfd, !is.na(Estimate) & r_sq >= 0.4999999  &length_time >=5 & System!="Marine" 
#            &Specific_location == 1 &!is.na(both_change)& !is.na(Bodymass_g)&(Class=="Aves"|Class=="Mammalia"))


# select_bm<-!df2$ID %in% df2bm$ID
# 
# df_nobm<-df2[select_bm,]
# unique(df_nobm$Binomial.x)
# 
# nrow(df2bm) -nrow(df2)
# 
# write.csv(df2, "bodymass_missing.csv")
#  

# df2<-subset(dfd, !is.na(Estimate) & r_sq >= 0.4999999  &length_time >=5 & System!="Marine"
#             &Specific_location == 1 & !is.na(Bodymass)&!is.na(both_change) &((Primary_threat =="Habitat degradation/change"|
#             Primary_threat=="Habitat loss"|Primary_threat=="Climate change")|
#             (Secondary_threat =="Habitat degradation/change"| Secondary_threat=="Habitat loss"|Secondary_threat=="Climate change")|
#             (Tertiary_threat == "Habitat degradation/change"| Tertiary_threat=="Habitat loss"|Tertiary_threat=="Climate change")))
# # 
# nrow(df2)
# 
# df2<-subset(dfd, !is.na(Estimate) & r_sq >= 0.4999999  &length_time >=5 & System!="Marine"
#             &Specific_location == 1 &!is.na(Bodymass)&!is.na(both_change) &((Primary_threat!="Disease"
#             & Primary_threat!="Exploitation"
#             &Primary_threat!="Invasive spp/genes"&Primary_threat!="Pollution") & (Secondary_threat!="Disease"&
#             Secondary_threat!="Exploitation"&Secondary_threat!="Invasive spp/genes"&Secondary_threat!="Pollution")
#             & (Tertiary_threat!="Disease"&Tertiary_threat!="Exploitation"&Tertiary_threat!="Invasive spp/genes"
#             &Tertiary_threat!="Pollution")))


df2$Estimate_sum<-df2$Estimate * df2$length_time
df2$both_change_sum<-df2$both_change * df2$length_time


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

parm_df<-sp_dups_df[,c("ID","Estimate", "both_change", "Bodymass_g")]  ##ID, land use, and climate  use "LUC_dist" or "Nat_change" for purely annual change in summed primary, secondary and other
parm_df<-sp_dups_df[,c("ID","Estimate_sum", "both_change_sum", "Bodymass_g")]  ##ID, land use, and climate  use "LUC_dist" or "Nat_change" for purely annual change in summed primary, secondary and other
#parm_df<-sp_dups_df[,c("ID","Estimate", "change_rate_49", "Bodymass_g")]

parm_mat<-as.matrix(parm_df)
parm_scale<-scale(parm_mat[,c("Estimate", "both_change", "Bodymass_g")])       #use the scaling factors at the bottom of these to scale the rasters
parm_scale<-scale(parm_mat[,c("Estimate_sum", "both_change_sum", "Bodymass_g")])       #use the scaling factors at the bottom of these to scale the rasters
#parm_scale<-scale(parm_mat[,c("Estimate", "change_rate_49", "Bodymass_g")])  

parm_id<-parm_mat[,"ID"]

parm_df_scale<-data.frame(parm_id,parm_scale)

colnames(parm_df_scale)<-c("ID","mean_slope_scale", "change_rate_scale", "Bodymass_scale")

sp_df_scale<-merge(sp_dups_df, parm_df_scale, by="ID")

dt<-data.table(sp_df_scale)

length(unique(dt$loc_id))

source("rsquaredglmm.R")

  library(lme4) 
  
  m0<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)

  m0a<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
  m0b<-lmer(lambda_sum ~ change_rate_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
  m0c<-lmer(lambda_sum ~ mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
  m0d<-lmer(lambda_sum ~ Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
  m1<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
  m1a<-lmer(lambda_sum ~ change_rate_scale+mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
  m1b<-lmer(lambda_sum ~ change_rate_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
  m1c<-lmer(lambda_sum ~ mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
  mnull<-lmer(lambda_sum ~ 1+(1|Binomial)+(1|loc_id),data=dt, REML=F)
  
 
  
  # #Weights
  library(MuMIn)
  
  #msAICc <- model.sel(m1,m1a,m1b,m1c,mnull)
  msAICc <- model.sel(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull)
  #msAICc <- model.sel(m1,m1a,m1b,m1c,mnull)
  msAICc$model<-rownames(msAICc)
  msAICc<-data.frame(msAICc)
  msAICc
  
  #Rsq
  models_list<-list(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull)
  #models_list<-list(m1,m1a,m1b,m1c,mnull)
  modelsR<-lapply(models_list,rsquared.glmm)
  modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)
  rownames(modelsRsq)<-c("m0","m0a","m0b","m0c","m0d","m1","m1a","m1b","m1c","mnull")

  modelsRsq
    
  library(MuMIn)
  var_imp<-summary(model.avg(models_list))

  mav<-model.avg(models_list, subset = cumsum <= 0.95)
  
  smav<-summary(mav)
  
  coef_av<-smav$coefmat.subset[2:5,"Estimate"]
  coef_df<-data.frame(coef_av)
  coef_df$lowCI<-confint(mav)[2:5,1]
  coef_df$highCI<-confint(mav)[2:5,2]
  coef_df
  coef_pcnt<-data.frame(((10^coef_df) - 1)*100)
  coef_pcnt
  
  # cnames<-c("MTC", "LUC", "LUC*MTC", "Bodymass")
  # 
  # rownames(coef_pcnt)<-cnames
  # 
  coef_pcnt$Var_name<-rownames(coef_pcnt)
  library(plotrix)
  
  plotCI(1:4, y=coef_pcnt$coef_av, ui=coef_pcnt$highCI, li=coef_pcnt$lowCI, ylab="Annual Population Change (%)", xlab="" ,xaxt = "n", 
         main="Birds and Mammals", lwd=1, ylim=c(min(coef_pcnt$lowCI*1.1), max(coef_pcnt$highCI*1.2)))
  axis(1, at=1:4, labels=rownames(coef_pcnt), las=2)
  abline(h=0, col="red", lty =2)

  #AIC
  AICs<-c(AIC(m0),AIC(m0a),AIC(m0b),AIC(m0c),AIC(m0d),AIC(m1),AIC(m1a),AIC(m1b),AIC(m1c),AIC(mnull))
  mnames<-c("LUC*MTC+BM","LUC+MTC+BM", "LUC+BM", "MTC+BM","BM","LUC*MTC","LUC+MTC", "LUC", "MTC", "NULL")
  
  AIC_diff<-AICs - AIC(mnull)
  del_AIC_df<-data.frame(AIC_diff, mnames)
  del_AIC_df<-del_AIC_df[order(del_AIC_df$AIC_diff),]
  del_AIC_df
  
  
library(coefplot)
library(ggplot2)

coef_pcnt$val<-1:4
coef_pcnt$var_name <- factor(coef_pcnt$Var_name, levels = coef_pcnt$Var_name[order(coef_pcnt$val)])
  
  p1<-ggplot(coef_pcnt)
  p1<- p1 + geom_hline(yintercept = 0, colour=gray(1/2), lty=2)
  p1<- p1 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=1.5, position = position_dodge(width=1/2))
  p1<- p1 + geom_pointrange(aes(x= Var_name, y=coef_av, ymin=lowCI, ymax=highCI), lwd=1, position=position_dodge(width=1/2), shape=21, fill="White")
  p1<- p1 + scale_y_continuous(breaks=seq(-8, 4, 2)) +theme_bw() + labs(y = "Annual Population Change (%)", x = "Variable") + theme(legend.title=element_blank(), text = element_text(size=20),axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
  print(p1)
  
  
coef_pcnt$Class<-"Birds"  
#coef_pcnt$Class<-"Mammals"  

coef_pcntb<-coef_pcnt
#coef_pcntm<-coef_pcnt  



coef_both<-rbind(coef_pcntb[,c(1,2,3,4,7)], coef_pcntm[,c(1,2,3,4,7)])

coef_both$Var_name
coef_both$Var_name[coef_both$Var_name == "mean_slope_scale"] <- "MTC"
coef_both$Var_name[coef_both$Var_name == "change_rate_scale"] <- "LUC"
coef_both$Var_name[coef_both$Var_name == "change_rate_scale:mean_slope_scale"] <- "MTC*LUC"
coef_both$Var_name[coef_both$Var_name == "Bodymass_scale"] <- "BM"

#write.csv(coef_both, "Model_Average_coefs3.csv")
coefs_both<-read.csv("Model_Average_coefs3.csv")

p1<-ggplot(coefs_both, aes(colour=Class))
p1<- p1 + geom_hline(yintercept = 0, colour=gray(1/2), lty=2)
p1<- p1 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=1.5, position = position_dodge(width=1/2))
p1<- p1 + geom_pointrange(aes(x= Var_name, y=coef_av, ymin=lowCI, ymax=highCI), lwd=1, position=position_dodge(width=1/2), shape=21, fill="White")
p1<- p1 + coord_flip()+ scale_y_continuous(breaks=seq(-8, 4, 2)) +theme_bw() + labs(y = "Annual Population Change (%)", x = "Variable") + theme(legend.title=element_blank(), text = element_text(size=20),axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
print(p1)

((10^msAICc[,2:5])-1)*100



  