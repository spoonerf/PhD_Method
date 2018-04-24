library(dplyr)
library(lme4)

temp<-read.csv("All_LPI_All_Years_Nobuff_1931_moreLPI_end2005.csv")

body<-read.csv("bird_and_mammal_traits2.csv")
body2<-read.csv("Bird_and_Mammal_BM.csv")

body3<-rbind(body, body2)

body4<-read.csv("LPI_BodyMass_Amniote_Database_edit.csv")
body4<-body4[,-3]
hyde<-read.csv("Hyde_crop_pasture_annual_change.csv")

hyde2<-read.csv("Hyde_crop_pasture_annual_change_sum.csv")

temp<-merge(temp, hyde[,c(2,3)], by="ID")

LPI<-read.csv("LPI_pops_20160523_edited.csv")
LPI_bin<-LPI[,c("ID", "Binomial")]

body4<-merge(body4[,c(2:4)], LPI_bin, by="Binomial")

body4<-unique(body4)


Realm<-read.csv("Realm.csv", na.strings="")
Realm<-Realm[,c("ID", "WWF_REALM2")]
Realm2<-read.csv("hub_realms.csv")
colnames(Realm2)<-c("ID", "WWF_REALM2")

Realm3<-rbind(Realm,Realm2) #sort this has duplicates
Realm<-Realm[!is.na(Realm3$WWF_REALM2),]

pop<-read.csv("Global_Population_Trends_Rsq_Lambda_07_10_16.csv")

temp<-temp[,c("ID", "Estimate")]

LPI<-LPI[,c("ID","Binomial","Confidential","Common_name", "Order","Family", "Protected_status", "Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat", "Migratory", "Forest", "Confidential")]

df<-merge(merge(temp,body4[,c(2:4)], by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)
dfd<-merge(df, hyde[,c(-1,-3)], by="ID")

head(dfd)

nrow(df)
nrow(dfd) 

#Birds
df2<-subset(dfd, !is.na(Estimate)  & r_sq >= 0.4999999 &length_time >=5& System!="Marine" 
            &Specific_location == 1 &!is.na(both_change) & !is.na(Log_Body_Mass_g)
            & (Class=="Aves") & Protected_status != "Unknown"  & Protected_status != "Both")
nrow(df2)

df2$Protected_status[df2$Protected_status == "No (area surrounding PA)"] <- "No"
df2$Protected_status[df2$Protected_status == "No (large survey area)"] <- "No"

table(df2$Class)


df2[is.na(df2$lambda_mean),]$lambda_mean<-0

library(dplyr)
#birds 
diet<-read.csv("birddiet.csv")
diet$Binomial<-gsub(" ", "_", diet$Scientific)
df2<-merge(df2, diet[,c(20,41)], by="Binomial")

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

library(lme4) 

m0<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0di<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0f<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0a<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0ad<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Bodymass_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0af<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0b<-lmer(lambda_mean ~ change_rate_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0bd<-lmer(lambda_mean ~ change_rate_scale+Bodymass_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0bf<-lmer(lambda_mean ~ change_rate_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0c<-lmer(lambda_mean ~ mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0cd<-lmer(lambda_mean ~ mean_slope_scale+Bodymass_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0cf<-lmer(lambda_mean ~ mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0d<-lmer(lambda_mean ~ Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0dd<-lmer(lambda_mean ~ Bodymass_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0df<-lmer(lambda_mean ~ Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1d<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1f<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1fd<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Protected_status+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1a<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1ad<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1af<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1afd<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Protected_status+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1b<-lmer(lambda_mean ~ change_rate_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1bd<-lmer(lambda_mean ~ change_rate_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1bf<-lmer(lambda_mean ~ change_rate_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)  
m1bfd<-lmer(lambda_mean ~ change_rate_scale+Protected_status+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)  

m1c<-lmer(lambda_mean ~ mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1cd<-lmer(lambda_mean ~ mean_slope_scale+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1cf<-lmer(lambda_mean ~ mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1cfd<-lmer(lambda_mean ~ mean_slope_scale+Protected_status+Diet.5Cat+(1|Binomial)+(1|loc_id),data=dt, REML=F)

mnull<-lmer(lambda_mean ~ 1+(1|Binomial)+(1|loc_id),data=dt, REML=F)



# #Weights
library(MuMIn)

msAICc <- model.sel(m0,m0di, m0f,m0a,m0ad,m0af,m0b,m0bd,m0bf, m0c,m0cd, m0cf, m0d, m0dd, m0df,m1, m1d,m1f,m1fd,m1a,m1ad,m1af, m1afd,m1b, m1bd, m1bf, m1bfd,m1c, m1cd, m1cf,m1cfd,mnull)
msAICc$model<-rownames(msAICc)
msAICc<-data.frame(msAICc)
msAICc

((10^msAICc[,c(1:5)]) - 1)*100

models_list<-list(m0,m0di, m0f,m0a,m0ad,m0af,m0b,m0bd,m0bf, m0c,m0cd, m0cf, m0d, m0dd, m0df,m1, m1d,m1f,m1fd,m1a,m1ad,m1af, m1afd,m1b, m1bd, m1bf, m1bfd,m1c, m1cd, m1cf,m1cfd,mnull)
modelsR<-lapply(models_list,rsquared.glmm)
modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)
rownames(modelsRsq)<-c("m0","m0di","m0f","m0a", "m0ad", "m0af","m0b", "m0bd", "m0bf","m0c","m0d", "m0cf","m0cd", "m0dd", "m0df","m1", "m1d", "m1f","m1fd", "m1a", "m1ad", "m1af", "m1afd","m1b","m1bd", "m1bf", "m1bfd","m1c", "m1cd", "m1cf", "m1cfd","mnull")
modelsRsq

library(MuMIn)
var_imp<-summary(model.avg(models_list))

mav<-model.avg(models_list, subset =  cumsum(weight) <= .95)

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

coef_pcnt$val<-1:length(coef_pcnt$Var_name)
coef_pcnt$var_name <- factor(coef_pcnt$Var_name, levels = coef_pcnt$Var_name[order(coef_pcnt$val)])

p1<-ggplot(coef_pcnt)
p1<- p1 + geom_hline(yintercept = 0, colour=gray(1/2), lty=2)
p1<- p1 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=1.5, position = position_dodge(width=1/2))
p1<- p1 + geom_pointrange(aes(x= Var_name, y=coef_av, ymin=lowCI, ymax=highCI), lwd=1, position=position_dodge(width=1/2), shape=21, fill="White")
p1<- p1 + scale_y_continuous(breaks=seq(-8, 4, 2)) +theme_bw() + labs(y = "Annual Population Change (%)", x = "Variable") + theme(legend.title=element_blank(), text = element_text(size=20),axis.title.x = element_text(margin = unit(c(5, 5, 0, 0), "mm")))
print(p1)



coef_pcnt$Class<-"Birds"  
coef_pcntb<-coef_pcnt




#mammals


df2<-subset(dfd, !is.na(Estimate)  & r_sq >= 0.4999999 &length_time >=5& System!="Marine" 
            &Specific_location == 1 &!is.na(both_change) & !is.na(Log_Body_Mass_g)
            & (Class=="Mammalia") & Protected_status != "Unknown"  & Protected_status != "Both")
nrow(df2)

df2$Protected_status[df2$Protected_status == "No (area surrounding PA)"] <- "No"
df2$Protected_status[df2$Protected_status == "No (large survey area)"] <- "No"

table(df2$Class)


df2[is.na(df2$lambda_mean),]$lambda_mean<-0

# #mammals
diet<-read.csv("mammaldiet.csv")
diet$Binomial<-paste(diet$Genus, diet$Species, sep= "_")
df2<-merge(df2, diet[,c(6,24,31)], by="Binomial")

nrow(df2)
library(dplyr)


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
m0di<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0f<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0a<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0ad<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Bodymass_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0af<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0b<-lmer(lambda_mean ~ change_rate_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0bd<-lmer(lambda_mean ~ change_rate_scale+Bodymass_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0bf<-lmer(lambda_mean ~ change_rate_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0c<-lmer(lambda_mean ~ mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0cd<-lmer(lambda_mean ~ mean_slope_scale+Bodymass_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0cf<-lmer(lambda_mean ~ mean_slope_scale+Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m0d<-lmer(lambda_mean ~ Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0dd<-lmer(lambda_mean ~ Bodymass_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m0df<-lmer(lambda_mean ~ Bodymass_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1d<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1f<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1fd<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Protected_status+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1a<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1ad<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1af<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1afd<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+Protected_status+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)

m1b<-lmer(lambda_mean ~ change_rate_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1bd<-lmer(lambda_mean ~ change_rate_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1bf<-lmer(lambda_mean ~ change_rate_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)  
m1bfd<-lmer(lambda_mean ~ change_rate_scale+Protected_status+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)  

m1c<-lmer(lambda_mean ~ mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1cd<-lmer(lambda_mean ~ mean_slope_scale+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1cf<-lmer(lambda_mean ~ mean_slope_scale+Protected_status+(1|Binomial)+(1|loc_id),data=dt, REML=F)
m1cfd<-lmer(lambda_mean ~ mean_slope_scale+Protected_status+TrophicLevel+(1|Binomial)+(1|loc_id),data=dt, REML=F)

mnull<-lmer(lambda_mean ~ 1+(1|Binomial)+(1|loc_id),data=dt, REML=F)



# #Weights
library(MuMIn)

msAICc <- model.sel(m0,m0di, m0f,m0a,m0ad,m0af,m0b,m0bd,m0bf, m0c,m0cd, m0cf, m0d, m0dd, m0df,m1, m1d,m1f,m1fd,m1a,m1ad,m1af, m1afd,m1b, m1bd, m1bf, m1bfd,m1c, m1cd, m1cf,m1cfd,mnull)
msAICc$model<-rownames(msAICc)
msAICc<-data.frame(msAICc)
msAICc

((10^msAICc[,c(1:5)]) - 1)*100

models_list<-list(m0,m0di, m0f,m0a,m0ad,m0af,m0b,m0bd,m0bf, m0c,m0cd, m0cf, m0d, m0dd, m0df,m1, m1d,m1f,m1fd,m1a,m1ad,m1af, m1afd,m1b, m1bd, m1bf, m1bfd,m1c, m1cd, m1cf,m1cfd,mnull)
modelsR<-lapply(models_list,rsquared.glmm)
modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)
rownames(modelsRsq)<-c("m0","m0di","m0f","m0a", "m0ad", "m0af","m0b", "m0bd", "m0bf","m0c","m0d", "m0cf","m0cd", "m0dd", "m0df","m1", "m1d", "m1f","m1fd", "m1a", "m1ad", "m1af", "m1afd","m1b","m1bd", "m1bf", "m1bfd","m1c", "m1cd", "m1cf", "m1cfd","mnull")
modelsRsq

library(MuMIn)
var_imp<-summary(model.avg(models_list))


mav<-model.avg(models_list, subset =  cumsum(weight) <= .95, fit=TRUE)

mav<-model.avg(models_list, subset =  delta < 2, fit=TRUE)

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

coef_pcnt$val<-1:length(coef_pcnt$Var_name)
coef_pcnt$var_name <- factor(coef_pcnt$Var_name, levels = coef_pcnt$Var_name[order(coef_pcnt$val)])

p1<-ggplot(coef_pcnt)
p1<- p1 + geom_hline(yintercept = 0, colour=gray(1/2), lty=2)
p1<- p1 + geom_linerange(aes(x=Var_name, ymin=lowCI, ymax=highCI), lwd=1.5, position = position_dodge(width=1/2))
p1<- p1 + geom_pointrange(aes(x= Var_name, y=coef_av, ymin=lowCI, ymax=highCI), lwd=1, position=position_dodge(width=1/2), shape=21, fill="White")
p1<- p1 + scale_y_continuous(breaks=seq(-8, 4, 2)) +theme_bw() + labs(y = "Annual Population Change (%)", x = "Variable") + theme(legend.title=element_blank(), text = element_text(size=20),axis.title.x = element_text(margin = unit(c(5, 5, 0, 0), "mm")))
print(p1)


coef_pcnt$Class<-"Mammals"  
coef_pcntm<-coef_pcnt  

