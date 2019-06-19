library(ggplot2)
library(dplyr)

ccf_cervus<-read.table("ccf_melt_cervus.csv")
ccf_melt_short<-read.table("ccf_melt_ursus2_short.csv")
ccf_melt<-read.table("ccf_melt_ursus2.csv")
ccf_bear<-rbind(ccf_melt, ccf_melt_short)

ccf_ibex<-read.table("ccf_melt_capra_new.csv")
ccf_ibex<-ccf_ibex[,-c(8,9)]

ccf_cervus$CName<-"Red deer"
ccf_bear$CName<-"Brown bear"
ccf_ibex$CName<-"Alpine ibex"

ccf_melt<-rbind(ccf_cervus, ccf_bear, ccf_ibex)

ccf_melt$ID<-as.factor(ccf_melt$ID)

ccf_melt$ID<- factor(ccf_melt$ID,levels = c('10695','10713','10717','10714','539','10694',
                                            '10692','10710','10718','10696','3478','3480','3479',
                                            '4357','4358','543','6555','11176','11178','11180'),ordered = TRUE)

ccf_melt$ldd<-as.factor(ccf_melt$ldd)
ccf_melt$SD<-as.factor(ccf_melt$SD)

head(ccf_melt)

get_lag<-function(x){
  lag_out<-as.numeric(strsplit(as.character(ccf_melt$variable[x]), " ")[[1]][3])
  return(lag_out)
}

lags<-lapply(1:length(ccf_melt$variable), get_lag)
ccf_melt$lag<-do.call("rbind", lags)


ccf_melt_cnd<-ccf_melt[grepl("cnd", ccf_melt$variable),]

ccf_melt_sdm<-ccf_melt[grepl("sdm", ccf_melt$variable),]

ccf_melt_sdm<-ccf_melt_sdm[,c("ID", "variable", "value", "CName", "lag")]

colnames(ccf_melt_sdm)<-c("ID", "variable", "sdm_value", "CName", "lag")

ccf_melt_sdm<-unique(ccf_melt_sdm)

ccf_melt<-merge(ccf_melt_cnd, ccf_melt_sdm, by = c("ID", "CName", "lag"))


ccf_melt<-ccf_melt[ccf_melt$ID != 4357 & ccf_melt$ID != 4358 & ccf_melt$ID != 11178,]

ccf_melt<-ccf_melt[complete.cases(ccf_melt),]

ggplot(ccf_melt[ccf_melt$lag == 0 ,],aes(x=CName, y=value, group = ID, fill = CName)) + 
  geom_boxplot(position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("#d8b365", "#f5f5f5","#5ab4ac"))+
  #geom_point(size = 2)+
  geom_point(position = position_dodge(width = 0.9), aes(y = sdm_value, group =ID), size = 5, colour = "red")+
  geom_hline(yintercept=0, linetype = "dashed")+
  geom_vline(xintercept=1.5, linetype = "dashed")+
  geom_vline(xintercept=2.5, linetype = "dashed")+
  theme_bw()+
  ylim(-1,1)+
  xlab("")+
  ylab("Correlation with Trends in Observed \nPopulation Growth Rates")+
  annotate("text", x = 0.55, y = 0.9, label = "A", fontface = 2, size = 20)+
  theme(legend.position="none",axis.text=element_text(size=20),
        axis.title=element_text(size=20))

ccf_melt$ldd_sd<-paste(ccf_melt$ldd, ccf_melt$SD, sep = "_")

ccf_max<-ccf_melt%>%
  group_by(ID,ldd_sd)%>%
  mutate(max_value = max(value), max_sdm_value = max(sdm_value))%>%
  dplyr:::select(ID,rep_id,N_used,lag,value,max_value, sdm_value,max_sdm_value, CName, ldd, SD)%>%
  distinct()

ccf_max<-ccf_max[ccf_max$value == ccf_max$max_value,]

ccf_max%>%
  group_by(ID)%>%
  mutate(mean_lag = mean(abs(lag)), sd_lag = sd(abs(lag)))%>%
  select(ID, mean_lag, sd_lag)%>%
  distinct()





ccf_max_sdm<-ccf_max[ccf_max$sdm_value == ccf_max$max_sdm_value,]
ccf_max_sdm<-unique(ccf_max_sdm[,c("ID", "lag","max_sdm_value")])


ccf_max_sdm<-ccf_max_sdm[,c("ID", "max_sdm_value")]

ccf_max<-ccf_max[complete.cases(ccf_max),]



ggplot(ccf_max,aes(x=CName, y=value, group = ID, fill = CName)) + 
  geom_boxplot(position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("#d8b365", "#f5f5f5","#5ab4ac"))+
  #geom_point(size = 2)+
  geom_point(position = position_dodge(width = 0.9), aes(y = max_sdm_value, group =ID), size = 5, colour = "red")+
  geom_hline(yintercept=0, linetype = "dashed")+
  geom_vline(xintercept=1.5, linetype = "dashed")+
  geom_vline(xintercept=2.5, linetype = "dashed")+
  theme_bw()+
  ylim(-1,1)+
  xlab("")+
  ylab("Maximum Correlation Coefficient")+
  annotate("text", x = 0.55, y = 0.9, label = "B", fontface = 2, size = 20)+
  theme(legend.position="none",axis.text=element_text(size=20),
        axis.title=element_text(size=20))





#####plotting lambdas

ursus_cnd_lambdas<-read.csv("melt_lambda_short_ursus.csv")
cervus_cnd_lambdas<-read.csv("melt_lambda_short_cervus.csv")
capra_cnd_lambdas<-read.csv("melt_lambda_short_capra.csv")

ursus_obs<-read.csv("all_year_ab_ursus.csv")
cervus_obs<-read.csv("all_year_ab_cervus.csv")
capra_obs<-read.csv("all_year_ab_capra.csv")

library(dplyr)
ursus_years<-ursus_obs%>%
  group_by(ID)%>%
  mutate(year_min = min(Year), year_max = max(Year))%>%
  select(ID, year_min, year_max)%>%
  distinct()

capra_years<-capra_obs%>%
  group_by(ID)%>%
  mutate(year_min = min(Year), year_max = max(Year))%>%
  select(ID, year_min, year_max)%>%
  distinct()


ursus_cnd_lambdas<-merge(ursus_cnd_lambdas, ursus_years)

capra_cnd_lambdas<-merge(capra_cnd_lambdas, capra_years)

# ursus_cnd_lambdas<-ursus_cnd_lambdas[,-2]
# cervus_cnd_lambdas<-cervus_cnd_lambdas[,-1]
# cervus_cnd_lambdas<-cervus_cnd_lambdas[,-8]
# capra_cnd_lambdas<-capra_cnd_lambdas[,-2]

ursus_cnd_lambdas$CName <-"Brown bear"
cervus_cnd_lambdas$CName <-"Red deer"
capra_cnd_lambdas$CName<-"Alpine ibex"

cnd_lambdas<-rbind(ursus_cnd_lambdas, cervus_cnd_lambdas, capra_cnd_lambdas)

#write.csv(cnd_lambdas, "cnd_lambdas_all.csv")


ursus_obs$CName<-"Brown bear"
cervus_obs$CName<-"Red deer"
capra_obs$CName<-"Alpine ibex"

obs_lambdas<-rbind(ursus_obs, cervus_obs, capra_obs)

#write.csv(obs_lambdas, "obs_lambdas.csv")
obs_lambdas<-read.csv("obs_lambdas.csv")
obs_lambdas<-obs_lambdas[,-c(1,2)]
cnd_lambdas<-read.csv("cnd_lambdas_all.csv")
cnd_lambdas<-cnd_lambdas[,-1]



library(dplyr)
sum_lambdas_obs<-obs_lambdas%>%
  group_by(ID)%>%
  mutate(sum_lambdas_obs = sum(Lambdas), mean_lambdas_obs = mean(Lambdas))%>%
  select(ID, mean_lambdas_obs)%>%
  distinct()

sum_lambdas_cnd<-cnd_lambdas%>%
  group_by(ID, ldd,SD, rep_id)%>%
  filter(Year >= year_min & Year <= year_max)%>%
  mutate(sum_lambda = sum(Lambdas), mean_lambda = mean(Lambdas))%>%
  ungroup()


sum_lambdas_cnd$ID<-as.factor(sum_lambdas_cnd$ID)
sum_lambdas_obs$ID<-as.factor(sum_lambdas_obs$ID)


sum_lambdas<-merge(sum_lambdas_cnd, sum_lambdas_obs, by="ID")

sum_lambdas$ID<- factor(sum_lambdas$ID,levels = c('10695','10713','10717','10714','539','10694',
                                            '10692','10710','10718','10696','3478','3480','3479',
                                            '4357','4358','543','6555','11176','11178','11180'),ordered = TRUE)

sum_lambdas<-sum_lambdas[complete.cases(sum_lambdas),]


#getting rid of the short brown bear time-series (they aren't in Ch.4 either)
sum_lambdas<-sum_lambdas[sum_lambdas$ID != 4357 &sum_lambdas$ID != 4358 &sum_lambdas$ID != 11178,]


sum_lambdas%>%
  filter(CName == "Brown bear")%>%
  group_by(ID)%>%
  summarize(median_av_lambda = median(mean_lambda),sd_av_lambda = sd(mean_lambda) ,obs_lambda = unique(mean_lambdas_obs))%>%
  summarize(mean_med_lambdas_pred = mean(median_av_lambda),sd_mean_lambdas_pred = sd(median_av_lambda),mean_mean_lambdas_obs = mean(obs_lambda),sd_mean_lambdas_obs = sd(obs_lambda))


# ggplot(data = sum_lambdas,aes(x=CName, y=sum_lambda, group = ID, fill = CName)) + 
#   geom_boxplot(position = position_dodge(width = 0.9))+
#   scale_fill_manual(values=c("#d8b365", "#f5f5f5","#5ab4ac"))+
#   #geom_point(size = 2)+
#   geom_point(position = position_dodge(width = 0.9), aes(y = sum_lambdas_obs, group =ID), size = 5, colour = "black")+
#   geom_hline(yintercept=0, linetype = "dashed")+
#   geom_vline(xintercept=1.5, linetype = "dashed")+
#   geom_vline(xintercept=2.5, linetype = "dashed")+
#   theme_bw()+
#   xlab("")+
#   ylab("Summed Lambdas")+
#   theme(legend.position="none",axis.text=element_text(size=16),
#         axis.title=element_text(size=20)) 
# 
# 

ch2_fit<-read.csv("fitted_values_mammals.csv")
ch2_fit<-ch2_fit[,-1]

ch2_vals<-ch2_fit[ch2_fit$ID %in% sum_lambdas$ID, c("ID", "fitted")]
ch2_vals<-dplyr::select(ch2_vals, ID, fitted)

ch2_noranef<-read.csv("ch5_pops_ch2_pred_no_ranef.csv")
ch2_noranef<-ch2_noranef[,-1]
colnames(ch2_noranef)<-c("ID", "preds_no_ranef")

ch2_noranef<-read.csv("ch5_pops_ch2_pred_no_ranef.csv")
ch2_noranef<-ch2_noranef[,-1]
ch2_noranef<-select(ch2_noranef, ID, preds_no_ranef)


sum_lambdas<-merge(sum_lambdas, ch2_vals, by="ID")
sum_lambdas<-merge(sum_lambdas, ch2_noranef, by="ID")


sum_lambdas$fitted_perc<-((10^sum_lambdas$fitted)-1)*100
sum_lambdas$mean_lambda_obs_perc<-((10^sum_lambdas$mean_lambdas_obs)-1)*100
sum_lambdas$noranef_perc<-((10^sum_lambdas$preds_no_ranef)-1)*100


sum_lambdas<-sum_lambdas %>%
  group_by(ID, ldd, SD)%>%
  mutate(median_lambda_cnd  = median(mean_lambda))


sum_lambdas$sign_same <- sign(sum_lambdas$mean_lambdas_obs) == sign(sum_lambdas$median_lambda_cnd)


chk<-sum_lambdas %>%
  select(ID, ldd, SD, mean_lambdas_obs, median_lambda_cnd, sign_same) %>%
  distinct()

sum(chk$sign_same)/487


library(ggplot2)

ggplot(data = sum_lambdas,aes(x=CName, y=((10^mean_lambda)-1)*100, group = ID, fill = CName)) + 
  geom_boxplot(position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("#d8b365", "#f5f5f5","#5ab4ac"))+
  #geom_point(size = 2)+
  geom_point(position = position_dodge(width = 0.9), aes(y = ((10^mean_lambdas_obs)-1)*100, group =ID), size = 5, colour = "black")+
  geom_point(position = position_dodge(width = 0.9), aes(y = ((10^fitted)-1)*100, group =ID), size = 5, colour = "red", shape = 2)+
 geom_point(position = position_dodge(width = 0.9), aes(y = ((10^preds_no_ranef)-1)*100, group =ID), size = 5, colour = "blue", shape = 3)+
  geom_hline(yintercept=0, linetype = "dashed")+
  geom_vline(xintercept=1.5, linetype = "dashed")+
  geom_vline(xintercept=2.5, linetype = "dashed")+
  theme_bw()+
  #ylim(c(-0.3,0.2))+
  xlab("")+
  ylab("Average Annual Rate\nof Population Change (%)")+
  #annotate("text", x = 0.55, y = 0.175, label = "A", fontface = 2, size = 20)+
  theme(legend.position="none",axis.text=element_text(size=20),
        axis.title=element_text(size=20))





m0_out<-sum_lambdas %>%
  group_by(ID)%>%
  mutate(diff_fit = abs(mean_lambda_obs_perc - fitted_perc), diff_noranef = abs(mean_lambda_obs_perc - noranef_perc)) %>%
  select(ID,CName, diff_fit, diff_noranef)%>%
  distinct()





ggplot(data = sum_lambdas,aes(x=CName, y=mean_lambda, group = ID, fill = CName)) + 
  geom_boxplot(position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("#d8b365", "#f5f5f5","#5ab4ac"))+
  #geom_point(size = 2)+
  geom_point(position = position_dodge(width = 0.9), aes(y = mean_lambdas_obs, group =ID), size = 5, colour = "black")+
  geom_point(position = position_dodge(width = 0.9), aes(y = fitted, group =ID), size = 5, colour = "red", shape = 2)+
  geom_hline(yintercept=0, linetype = "dashed")+
  geom_vline(xintercept=1.5, linetype = "dashed")+
  geom_vline(xintercept=2.5, linetype = "dashed")+
  theme_bw()+
  #ylim(c(-0.3,0.2))+
  xlab("")+
  ylab("Average Annual Rate\nof Population Change (%)")+
  #annotate("text", x = 0.55, y = 0.175, label = "A", fontface = 2, size = 20)+
  theme(legend.position="none",axis.text=element_text(size=20),
        axis.title=element_text(size=20))








sum_lambdas_bear<-sum_lambdas[sum_lambdas$CName == "Brown bear",]


sum_lambdas_bear$ldd_sd<-sum_lambdas_bear$SD* sum_lambdas_bear$ldd

library(ggplot2)
ggplot(data = sum_lambdas_bear, aes(x = ldd, y = mean_lambda, group = ID,colour = ID ,fill= ID))+
  geom_point()


library(ggplot2)
ggplot(data = sum_lambdas_bear, aes(x = SD, y = mean_lambda, group = ID,colour = ID ,fill= ID))+
  geom_point()



plot(sum_lambdas_bear$ldd, sum_lambdas_bear$mean_lambda)
plot(sum_lambdas_bear$ldd_sd, sum_lambdas_bear$mean_lambda)


pops<-unique(sum_lambdas$ID)
library(akima)

names<-unique(sum_lambdas$CName)

#for (x in pops){
for (x in names){
sum_lambdas_bear<-sum_lambdas[sum_lambdas$CName == x,]

CName = unique(sum_lambdas_bear$CName)
  mean_mean_lambdas_bear<-sum_lambdas_bear%>%
  group_by(ldd, SD)%>%
  summarise(mean_lambda = mean(mean_lambda))

im <- with(mean_mean_lambdas_bear,interp(ldd,SD,mean_lambda))
vals<-t(im$z[,seq(from=ncol(im$z), to = 1, by = -1)])
ras<-raster(ncol = 40, nrow = 40)
extent(ras)<-extent(0,0.5, 0,0.5)
values(ras)<-((10^vals)-1)*100

test_spdf <- as(ras, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8)+
  scale_fill_viridis(name = "Average Population \nGrowth Rate (%)")+
  theme_bw()+
  xlab("Dispersal Rate")+
  ylab("Transition Matrix Stochasticity")+
  labs(title = paste(x, sep=""))+ theme(axis.text=element_text(size=20),
                                        axis.title=element_text(size=20), title = element_text(size= 20),legend.text=element_text(size=15))



}


test_spdf <- as(ras, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

library(viridis)

ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8)+
  scale_fill_viridis(name = "Average Population \nGrowth Rate (%)")+
  theme_bw()+
  xlab("Dispersal Rate")+
  ylab("Transition Matrix Stochasticity")+
  main(x)



plot(ras,xlab = "Dispersal Rate", ylab = "Stochasticity",  col=rev(terrain.colors(255))) 
#plot(ras, smallplot=c(.65, .17, .5, .85), legend.only=TRUE) 

pairs(sum_lambdas_bear$mean_lambda~ sum_lambdas_bear$ldd + sum_lambdas_bear$SD)


###all pops
sum_lambdas$ldd_sd<-sum_lambdas$SD* sum_lambdas$ldd


plot(sum_lambdas$ldd, sum_lambdas$mean_lambda)
plot(sum_lambdas$SD, sum_lambdas$mean_lambda)
plot(sum_lambdas$ldd_sd, sum_lambdas$mean_lambda)

sum_lambdas$CName<-as.factor(sum_lambdas$CName)
sum_lambdas$ID<-as.factor(sum_lambdas$ID)

library(ggplot2)
ggplot(data = sum_lambdas, aes(x = ldd, y = mean_lambda, group = ID,colour = ID ,fill= ID))+
  geom_point()


library(lme4)


m1<-lmer(sum_lambdas$mean_lambda~ sum_lambdas$ldd*sum_lambdas$SD + (1|sum_lambdas$ID))

pairs(sum_lambdas$mean_lambda~ sum_lambdas$ldd + sum_lambdas$SD)



mean_mean_lambdas<-sum_lambdas%>%
  filter(CName == "Brown bear")%>%
  group_by(ldd, SD)%>%
  summarise(mean_lambda = mean(mean_lambda))

library(akima)
im <- with(mean_mean_lambdas,interp(ldd,SD,mean_lambda))
im_out<-graphics:::image(im$x,im$y,im$z, col = topo.colors(10, alpha = 1), useRaster = T, xlab = "ldd", ylab = "SD")



#mean median correlations across pops

ccf_melt0 <- ccf_melt[ccf_melt$lag == 0,]

ccf_melt0 %>%
  group_by(ID)%>%
  filter(CName == "Red deer")%>%
    summarise(median_cor = median(value), median_hsm_cor = mean(sdm_value))%>%
  summarise(mean_median_cor = mean(median_cor),sd_median_cor = sd(median_cor), mean_median_hsm_cor = mean(median_hsm_cor),sd_median_hsm_cor = sd(median_hsm_cor))


ccf_melt0 %>%
  group_by(ID)%>%
  mutate(mean_cnd_cor = mean(value), sd_cnd_cor = sd(value) )%>%
  select(ID,sdm_value, mean_cnd_cor,sd_cnd_cor)%>%
  distinct()


ccf_max %>%
  group_by(ID)%>%
  summarise(mean_cor = mean(max_value), sd_cnd = sd(max_value), mean_hsm_cor = mean(max_sdm_value))



species_list<-list(c(10692, 10694))
species_list<-unlist(species_list)

dplyr::filter(ccf_max,ID %in% species_list)


years<-1950:2005

capra_sdm<-stack(paste("D:/Fiona/Git_Method/Git_Method/Legion/snow_capra_bias/hyde_weighted_ensemble_sdm_", years, ".tif", sep=""))
cervus_sdm<-stack(paste("D:/Fiona/Git_Method/Git_Method/Legion/snow_cervus_bias_faster/hyde_weighted_ensemble_sdm_", years, ".tif", sep=""))
bear_sdm<-stack(paste("D:/Fiona/Git_Method/Git_Method/Legion/snow_bear_bias_faster/hyde_weighted_ensemble_sdm_", years, ".tif", sep=""))

names(capra_sdm)<-years

capra_sdm[is.na(capra_sdm[])] <- 0 
cervus_sdm[is.na(cervus_sdm[])] <- 0 
bear_sdm[is.na(bear_sdm[])] <- 0 

ai<- as.matrix(capra_sdm)
aim<-melt(ai)
aim<-aim[,-1]
aim<-cbind("Alpine ibex", aim)
colnames(aim)<-c("CName", "Year", "HS")


names(cervus_sdm)<-years

ce<- as.matrix(cervus_sdm)
cem<-melt(ce)
cem<-cem[,-1]
cem<-cbind("Cervus elaphus", cem)
colnames(cem)<-c("CName", "Year", "HS")


names(bear_sdm)<-years

ua<- as.matrix(bear_sdm)
uam<-melt(ua)
uam<-uam[,-1]
uam<-cbind("Ursus arctos", uam)
colnames(uam)<-c("CName", "Year", "HS")


all_hs<-rbind(aim, cem, uam)

all_hs$Year<-(as.numeric(gsub("X", "",all_hs$Year)))


ggplot(all_hs, aes(x = Year, y =HS, group = CName, colour = CName))+
  geom_smooth()


model_diff<-sum_lambdas %>%
  group_by(ID) %>%
  mutate(median_cnd = median(mean_lambda))%>%
  summarize(cnd_diff = unique(abs(median_cnd - mean_lambdas_obs)), lme_diff = unique(abs(fitted -mean_lambdas_obs)))

model_diff<-sum_lambdas %>%
  group_by(ID) %>%
  summarize(cnd_diff = min(abs(mean_lambda - mean_lambdas_obs)), lme_diff = unique(abs(fitted -mean_lambdas_obs)))





#https://stackoverflow.com/questions/38173544/how-to-calculate-p-values-from-cross-correlation-function-in-r


ccf_melt0$sig<-(2 * (1 - pnorm(abs(ccf_melt0$value), mean = 0, sd = 1/sqrt(ccf_melt0$N_used)))<0.05)*1

ccf_melt0$sdm_sig<-(2 * (1 - pnorm(abs(ccf_melt0$sdm_value), mean = 0, sd = 1/sqrt(ccf_melt0$N_used)))<0.05)*1


library(dplyr)
ccf_melt0%>%
  group_by(CName)%>%
  filter(value >0 &sig ==1)%>%
  mutate(num_sig = sum(sig))%>%
  select(CName, num_sig)%>%
  distinct()
# 


ccf_max$sig<-(2 * (1 - pnorm(abs(ccf_max$max_value), mean = 0, sd = 1/sqrt(ccf_max$N_used)))<0.05)*1


ccf_max$sdm_sig<-(2 * (1 - pnorm(abs(ccf_max$max_sdm_value), mean = 0, sd = 1/sqrt(ccf_max$N_used)))<0.05)*1



library(dplyr)
ccf_max%>%
  group_by(CName)%>%
  filter(max_value >0 &sig ==1)%>%
  mutate(num_sig = sum(sig))%>%
  select(CName, num_sig)%>%
  distinct()
# 


library(dplyr)
ccf_max%>%
  group_by(CName)%>%
  select(CName, sdm_sig)%>%
  distinct()

  filter(max_value >0 &sig ==1)%>%
  mutate(num_sig = sum(sig))%>%
  
  
# 


cor_diff<-ccf_melt0 %>%
  group_by(ID) %>%
  mutate(median_cnd = median(value))%>%
  select(ID, median_cnd, sdm_value)%>%
  distinct()


cor_diff_max<-ccf_max%>%
  group_by(ID) %>%
  mutate(median_cnd_max = median(max_value))%>%
  select(ID, median_cnd_max, max_sdm_value)%>%
  distinct()


df<-merge(cor_diff, cor_diff_max)

colMeans(df[,2:5])


best_cnd<-sum_lambdas %>%
  group_by(ID)%>%
  mutate(diff = abs(mean_lambdas_obs - mean_lambda))%>%
  mutate(min_diff = min(diff))%>%
  mutate(mean_lambda_perc = ((10^mean_lambda)-1)*100)%>%
  select(ID, diff, min_diff,mean_lambda_perc)%>%
  distinct()

best_cnd[best_cnd$diff == best_cnd$min_diff,]




sum_lambdas %>%
  group_by(ID)%>%
  filter(CName == "Alpine ibex")%>%
  mutate(grow_cnd = sum(mean_lambda>0))%>%
  select(ID, grow_cnd)%>%
  distinct()


sum_lambdas %>%
  group_by(ID)%>%
  mutate(mean_trans = (((10^mean_lambdas_obs)-1)*100), median_lambda_cnd = ((10^median(mean_lambda)-1)*100), sd_lambda_cnd = sd(mean_lambda))%>%
  select(ID, mean_lambdas_obs, fitted_perc,mean_trans, median_lambda_cnd)%>%
    distinct()
  
  



