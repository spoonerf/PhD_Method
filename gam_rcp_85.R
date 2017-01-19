df_bird<-read.csv("bird_data_for_prediction.csv")
df_mammal<-read.csv("mammal_data_for_prediction.csv")

df_bird[is.na(df_bird$lambda_mean),]$lambda_mean<-0
df_mammal[is.na(df_mammal$lambda_mean),]$lambda_mean<-0

library(plyr)
#counting duplicates at each location
sp_dups_bird<-data.frame(ddply(df_bird,.(Longitude,Latitude),nrow))
sp_dups_mamm<-data.frame(ddply(df_mammal,.(Longitude,Latitude),nrow))

sp_dups_bird$loc_id<-1:length(sp_dups_bird$Longitude)
sp_dups_mamm$loc_id<-1:length(sp_dups_mamm$Longitude)

sp_dups_df_bird<-merge(sp_dups_bird, df_bird, by=c("Longitude","Latitude"))
sp_dups_df_mamm<-merge(sp_dups_mamm, df_mammal, by=c("Longitude","Latitude"))

library(data.table)
dt_bird = as.data.table(sp_dups_df_bird)
dt_mamm = as.data.table(sp_dups_df_mamm)

parm_df_bird<-sp_dups_df_bird[,c("ID","Estimate", "both_change", "Bodymass_g")]  ##ID, land use, and climate  use "LUC_dist" or "Nat_change" for purely annual change in summed primary, secondary and other
parm_df_mamm<-sp_dups_df_mamm[,c("ID","Estimate", "both_change", "Bodymass_g")]  ##ID, land use, and climate  use "LUC_dist" or "Nat_change" for purely annual change in summed primary, secondary and othe

parm_mat_bird<-as.matrix(parm_df_bird)
parm_mat_mamm<-as.matrix(parm_df_mamm)

parm_scale_bird<-scale(parm_mat_bird[,c("Estimate", "both_change", "Bodymass_g")])       #use the scaling factors at the bottom of these to scale the rasters
parm_scale_mamm<-scale(parm_mat_mamm[,c("Estimate", "both_change", "Bodymass_g")])       #use the scaling factors at the bottom of these to scale the rasters

parm_id_bird<-parm_mat_bird[,"ID"]
parm_id_mamm<-parm_mat_mamm[,"ID"]

parm_df_scale_bird<-data.frame(parm_id_bird,parm_scale_bird)
parm_df_scale_mamm<-data.frame(parm_id_mamm,parm_scale_mamm)

colnames(parm_df_scale_bird)<-c("ID","mean_slope_scale", "change_rate_scale", "Bodymass_scale")
colnames(parm_df_scale_mamm)<-c("ID","mean_slope_scale", "change_rate_scale", "Bodymass_scale")

sp_df_scale_bird<-merge(sp_dups_df_bird, parm_df_scale_bird, by="ID")
sp_df_scale_mamm<-merge(sp_dups_df_mamm, parm_df_scale_mamm, by="ID")

dt_bird<-data.table(sp_df_scale_bird)
dt_mamm<-data.table(sp_df_scale_mamm)

m1c_b<-lmer(lambda_mean ~ mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt_bird, REML=F)
m1c_m<-lmer(lambda_mean ~ mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt_mamm, REML=F)

m1_m<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+(1|Binomial)+(1|loc_id),data=dt_mamm, REML=F)


###############################################
library(raster)

CR40s<-brick("cru_ts3.23.1941.1950.tmp.dat.nc")
CR50s<-brick("cru_ts3.23.1951.1960.tmp.dat.nc")
CR60s<-brick("cru_ts3.23.1961.1970.tmp.dat.nc")
CR70s<-brick("cru_ts3.23.1971.1980.tmp.dat.nc")
CR80s<-brick("cru_ts3.23.1981.1990.tmp.dat.nc")
CR90s<-brick("cru_ts3.23.1991.2000.tmp.dat.nc")
CR00s<-brick("cru_ts3.23.2001.2010.tmp.dat.nc")

obs<-stack(CR40s[[109:120]], CR50s, CR60s, CR70s, CR80s, CR90s,CR00s[[1:48]] )

###85 rcp

######GISS E2

GISS85_06_25<-brick("ts_Amon_GISS-E2-R-CC_rcp85_r1i1p1_200601-202512.nc")
GISS85_26_50<-brick("ts_Amon_GISS-E2-R-CC_rcp85_r1i1p1_202601-205012.nc")
GISS85_51_75<-brick("ts_Amon_GISS-E2-R-CC_rcp85_r1i1p1_205101-207512.nc")
GISS85_76_00<-brick("ts_Amon_GISS-E2-R-CC_rcp85_r1i1p1_207601-210012.nc")

GISS85_06_00<-stack(GISS85_06_25,GISS85_26_50,GISS85_51_75,GISS85_76_00)

modelg_cc<-rotate(GISS85_06_00)

obs_g<-resample(obs, modelg_cc)
obs_gk<-obs_g+273.15

mask_mod_g<-mask(modelg_cc, obs_gk[[1]])

obs_mod_g<-stack(obs_gk, mask_mod_g)

obs_mod_mon_mean_g<-matrix(cellStats(obs_mod_g, mean))

plot(obs_mod_mon_mean_g)

#####HadGEM ES

ts85_05_30<-brick("ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_200512-203011.nc")
ts85_30_55<-brick("ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_203012-205511.nc")
ts85_55_80<-brick("ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_205512-208011.nc")
ts85_80_00<-brick("ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_208012-210011.nc")

ts85_05_00<-stack(ts85_05_30[[2:300]],ts85_30_55,ts85_55_80,ts85_80_00)
#ts_05_00<-brick(ts_05_00)

modelcc<-rotate(ts85_05_00)

obs_r<-resample(obs, modelcc)

obs_rk<-obs_r+273.15

mask_mod<-mask(modelcc, obs_rk[[1]])

#obs_mod<-stack(obs_rk, modelcc)

obs_mod<-stack(obs_rk, mask_mod)

obs_mod_mon_mean<-matrix(cellStats(obs_mod, mean))

plot(obs_mod_mon_mean)

#############################################

tsl85_05_30<-brick("tsl_Lmon_HadGEM2-ES_rcp85_r1i1p1_200512-203011.nc")
tsl85_30_55<-brick("tsl_Lmon_HadGEM2-ES_rcp85_r1i1p1_203012-205511.nc")
tsl85_55_80<-brick("tsl_Lmon_HadGEM2-ES_rcp85_r1i1p1_205512-208011.nc")
tsl85_80_99<-brick("tsl_Lmon_HadGEM2-ES_rcp85_r1i1p1_208012-209912.nc")

tsl85_05_00<-stack(tsl85_05_30[[2:300]],tsl85_30_55,tsl85_55_80,tsl85_80_99)
#ts_05_00<-brick(ts_05_00)

modelcc_tsl<-rotate(tsl85_05_00)

obs_r_tsl<-resample(obs, modelcc_tsl)

obs_rk_tsl<-obs_r_tsl+273.15

mask_mod_tsl<-mask(modelcc_tsl, obs_rk_tsl[[1]])

#obs_mod<-stack(obs_rk, modelcc)

obs_mod_tsl<-stack(obs_rk_tsl, mask_mod_tsl)

obs_mod_mon_mean_tsl<-matrix(cellStats(obs_mod_tsl, mean))

plot(obs_mod_mon_mean_tsl)

##################################################
year_mon<-rep(1:(nlayers(obs_mod)/12), each=12)
year_mon_g<-rep(1:(nlayers(obs_mod_g)/12), each=12)
year_mon_tsl<-rep(1:(nlayers(obs_mod_tsl)/12), each=12)

#changing from monthly to annual data
obs_mod_ann<-stackApply(obs_mod, indices=year_mon, fun=mean, na.rm=TRUE)
obs_mod_ann_g<-stackApply(obs_mod_g, indices=year_mon_g, fun=mean, na.rm=TRUE)
obs_mod_ann_tsl<-stackApply(obs_mod_tsl, indices=year_mon_tsl, fun=mean, na.rm=TRUE)

obs_mat_ann<-matrix(cellStats(obs_mod_ann, mean))
obs_mat_ann_g<-matrix(cellStats(obs_mod_ann_g, mean))
obs_mat_ann_tsl<-matrix(cellStats(obs_mod_ann_tsl, mean))
plot(obs_mat_ann_tsl)


########################################
x = c(1:149)
g = mgcv:::gam(obs_mat_ann~s(x, k=8), fx=TRUE)
plot(g)

g_pred = predict(g)
g_pred = g_pred[55:148]
d = diff(g_pred)
plot(g_pred)
plot(d)
####################################

x<-c(1:150)
gg = mgcv:::gam(obs_mat_ann_g~s(x, k=8), fx=TRUE)
plot(gg)

gg_pred = predict(gg)
gg_pred = gg_pred[55:150]
dg = diff(gg_pred)
plot(gg_pred)
plot(dg)
####################################

x = c(1:149)
gt = mgcv:::gam(obs_mat_ann_tsl~s(x, k=8), fx=TRUE)
plot(gt)

g_predt = predict(gt)
g_predt = g_predt[55:148]
dt = diff(g_predt)
plot(g_predt)
plot(dt)


#bird values
centre_tempb<-0.04526127
scale_tempb<-0.07006285

#mammal values
centre_tempm<-
scale_tempm<-

#both values

centre_temp<-0.3605663
scale_temp<-0.6856779


gam_diff_b<-scale(d, center=centre_tempb, scale=scale_tempb)
plot(gam_diff_b)

gam_diff_m<-scale(d, center=centre_tempm, scale=scale_tempm)
plot(gam_diff_m)

gam_diff_bg<-scale(dg, center=centre_tempb, scale=scale_tempb)
plot(gam_diff_b)

gam_diff_mg<-scale(dg, center=centre_tempm, scale=scale_tempm)
plot(gam_diff_m)

gam_diff_bt<-scale(dt, center=centre_tempb, scale=scale_tempb)
plot(gam_diff_bt)

gam_diff_mt<-scale(dt, center=centre_tempm, scale=scale_tempm)
plot(gam_diff_mt)


######################################
ave_mod_b<-data.frame(gam_diff_b, rep(0, length(gam_diff_b)), rep(0, length(gam_diff_b)))
colnames(ave_mod_b)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred_b<-predict(m1c_b, ave_mod_b, re.form=NA, se.fit=T)

ave_mod_m<-data.frame(gam_diff_m, rep(0, length(gam_diff_m)), rep(0, length(gam_diff_m)))
colnames(ave_mod_m)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred_m<-predict(m1_m, ave_mod_m, re.form=NA, se.fit=T)

######################################
ave_mod_bg<-data.frame(gam_diff_bg, rep(0, length(gam_diff_b)), rep(0, length(gam_diff_b)))
colnames(ave_mod_bg)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred_bg<-predict(m1c_b, ave_mod_bg, re.form=NA, se.fit=T)

ave_mod_mg<-data.frame(gam_diff_mg, rep(0, length(gam_diff_m)), rep(0, length(gam_diff_m)))
colnames(ave_mod_mg)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred_mg<-predict(m1_m, ave_mod_mg, re.form=NA, se.fit=T)

########################################
ave_mod_bt<-data.frame(gam_diff_bt, rep(0, length(gam_diff_bt)), rep(0, length(gam_diff_bt)))
colnames(ave_mod_bt)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred_bt<-predict(m1c_b, ave_mod_bt, re.form=NA, se.fit=T)

ave_mod_mt<-data.frame(gam_diff_mt, rep(0, length(gam_diff_mt)), rep(0, length(gam_diff_mt)))
colnames(ave_mod_mt)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred_mt<-predict(m1_m, ave_mod_mt, re.form=NA, se.fit=T)



plot(pop_pred_b$fit, ylim=c(-0.045, 0.005), main="Birds", type="l")
lines(pop_pred_b$fit + 1.96*pop_pred_b$se.fit, col="red", lty=2)
lines(pop_pred_b$fit - 1.96*pop_pred_b$se.fit, col="red", lty=2)

plot(pop_pred_m$fit, ylim=c(-0.015, 0.01), main="Mammals",  type="l")
lines(pop_pred_m$fit + 1.96*pop_pred_m$se.fit, col="red", lty=2)
lines(pop_pred_m$fit - 1.96*pop_pred_m$se.fit, col="red", lty=2)

#####################################
plot(pop_pred_bg$fit, ylim=c(-0.045, 0.005), main="Birds", type="l")
lines(pop_pred_bg$fit + 1.96*pop_pred_bg$se.fit, col="red", lty=2)
lines(pop_pred_bg$fit - 1.96*pop_pred_bg$se.fit, col="red", lty=2)

plot(pop_pred_mg$fit, ylim=c(-0.015, 0.01), main="Mammals",  type="l")
lines(pop_pred_mg$fit + 1.96*pop_pred_mg$se.fit, col="red", lty=2)
lines(pop_pred_mg$fit - 1.96*pop_pred_mg$se.fit, col="red", lty=2)

#####################################
plot(pop_pred_bt$fit, ylim=c(-0.045, 0.005), main="Birds", type="l")
lines(pop_pred_bt$fit + 1.96*pop_pred_bt$se.fit, col="red", lty=2)
lines(pop_pred_bt$fit - 1.96*pop_pred_bt$se.fit, col="red", lty=2)

plot(pop_pred_mt$fit, ylim=c(-0.015, 0.01), main="Mammals",  type="l")
lines(pop_pred_mt$fit + 1.96*pop_pred_mt$se.fit, col="red", lty=2)
lines(pop_pred_mt$fit - 1.96*pop_pred_mt$se.fit, col="red", lty=2)

#######################################
10^sum(pop_pred_b$fit + 1.96*pop_pred_b$se.fit)
10^sum(pop_pred_b$fit - 1.96*pop_pred_b$se.fit)
10^sum(pop_pred_b$fit)

10^sum(pop_pred_m$fit + 1.96*pop_pred_m$se.fit)
10^sum(pop_pred_m$fit - 1.96*pop_pred_m$se.fit)
10^sum(pop_pred_m$fit)
######################################
10^sum(pop_pred_bg$fit + 1.96*pop_pred_bg$se.fit)
10^sum(pop_pred_bg$fit - 1.96*pop_pred_bg$se.fit)
10^sum(pop_pred_bg$fit)

10^sum(pop_pred_mg$fit + 1.96*pop_pred_mg$se.fit)
10^sum(pop_pred_mg$fit - 1.96*pop_pred_mg$se.fit)
10^sum(pop_pred_mg$fit)
#######################################
10^sum(pop_pred_bt$fit + 1.96*pop_pred_bt$se.fit)
10^sum(pop_pred_bt$fit - 1.96*pop_pred_bt$se.fit)
10^sum(pop_pred_bt$fit)

10^sum(pop_pred_mt$fit + 1.96*pop_pred_mt$se.fit)
10^sum(pop_pred_mt$fit - 1.96*pop_pred_mt$se.fit)
10^sum(pop_pred_mt$fit)


#################################
pop_pred_b$fit2<-c(0, pop_pred_b$fit)
pop_pred_b$se.fit2<-c(0, pop_pred_b$se.fit)

pop_pred_m$fit2<-c(0, pop_pred_m$fit)
pop_pred_m$se.fit2<-c(0, pop_pred_m$se.fit)

###########################################
pop_pred_bg$fit2<-c(0, pop_pred_bg$fit)
pop_pred_bg$se.fit2<-c(0, pop_pred_bg$se.fit)

pop_pred_mg$fit2<-c(0, pop_pred_mg$fit)
pop_pred_mg$se.fit2<-c(0, pop_pred_mg$se.fit)

############################################

pop_pred_bt$fit2<-c(0, pop_pred_bt$fit)
pop_pred_bt$se.fit2<-c(0, pop_pred_bt$se.fit)

pop_pred_mt$fit2<-c(0, pop_pred_mt$fit)
pop_pred_mt$se.fit2<-c(0, pop_pred_mt$se.fit)

##########################################

plot(y=10^cumsum(pop_pred_b$fit2), x=2006:2099, ylim=c(0, 1.5), type="l", ylab="Bird Population Index", xlab="Year", main="Predicted Bird Population Index Under Scenario RCP 8.5")
lines(y=10^cumsum(pop_pred_b$fit2 - 1.96*pop_pred_b$se.fit2), x=2006:2099, lty=2, col="Red")
lines(y=10^cumsum(pop_pred_b$fit2 + 1.96*pop_pred_b$se.fit2), x=2006:2099, lty=2, col="Red")

y<-10^cumsum(pop_pred_b$fit2)
x<-2006:2099
lci<-10^cumsum(pop_pred_b$fit2 - 1.96*pop_pred_b$se.fit2)
uci<-10^cumsum(pop_pred_b$fit2 + 1.96*pop_pred_b$se.fit2)
  
df<-data.frame(x,y, lci, uci)

ggplot(df, aes(x,y))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.3, colour=NA)+
  labs(y = "Bird Population Index", x = "")+
  scale_y_continuous(breaks=seq(0, 1, 0.25))+
  theme_bw()+
  theme(text = element_text(size=20), axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")))

##########GISS


plot(y=10^cumsum(pop_pred_bg$fit2), x=2005:2100, ylim=c(0, 1.5), type="l", ylab="Bird Population Index", xlab="Year", main="Predicted Bird Population Index Under Scenario RCP 8.5")
lines(y=10^cumsum(pop_pred_bg$fit2 - 1.96*pop_pred_bg$se.fit2), x=2005:2100, lty=2, col="Red")
lines(y=10^cumsum(pop_pred_bg$fit2 + 1.96*pop_pred_bg$se.fit2), x=2005:2100, lty=2, col="Red")

y<-10^cumsum(pop_pred_bg$fit2)
x<-2005:2100
lci<-10^cumsum(pop_pred_bg$fit2 - 1.96*pop_pred_bg$se.fit2)
uci<-10^cumsum(pop_pred_bg$fit2 + 1.96*pop_pred_bg$se.fit2)

dfg<-data.frame(x,y, lci, uci)

ggplot(dfg, aes(x,y))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.3, colour=NA)+
  labs(y = "Bird Population Index", x = "")+
  scale_y_continuous(breaks=seq(0, 1, 0.25))+
  theme_bw()+
  theme(text = element_text(size=20), axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")))

##########################################tsl

plot(y=10^cumsum(pop_pred_bt$fit2), x=2006:2099, ylim=c(0, 1.5), type="l", ylab="Bird Population Index", xlab="Year", main="Predicted Bird Population Index Under Scenario RCP 8.5")
lines(y=10^cumsum(pop_pred_bt$fit2 - 1.96*pop_pred_bt$se.fit2), x=2006:2099, lty=2, col="Red")
lines(y=10^cumsum(pop_pred_bt$fit2 + 1.96*pop_pred_bt$se.fit2), x=2006:2099, lty=2, col="Red")

y<-10^cumsum(pop_pred_bt$fit2)
x<-2006:2099
lci<-10^cumsum(pop_pred_bt$fit2 - 1.96*pop_pred_bt$se.fit2)
uci<-10^cumsum(pop_pred_bt$fit2 + 1.96*pop_pred_bt$se.fit2)

dft<-data.frame(x,y, lci, uci)

ggplot(dft, aes(x,y))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.3, colour=NA)+
  labs(y = "Bird Population Index", x = "")+
  scale_y_continuous(breaks=seq(0, 1, 0.25))+
  theme_bw()+
  theme(text = element_text(size=20), axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")))
###########################################tsl mammals


plot(y=10^cumsum(pop_pred_mt$fit2), x=2006:2099, ylim=c(0, 1.5), type="l", ylab="Mammal Population Index", xlab="Year", main="Predicted Mammal Population Index Under Scenario RCP 8.5")
lines(y=10^cumsum(pop_pred_mt$fit2 - 1.96*pop_pred_mt$se.fit2), x=2006:2099, lty=2, col="Red")
lines(y=10^cumsum(pop_pred_mt$fit2 + 1.96*pop_pred_mt$se.fit2), x=2006:2099, lty=2, col="Red")

y<-10^cumsum(pop_pred_mt$fit2)
x<-2006:2099
lci<-10^cumsum(pop_pred_mt$fit2 - 1.96*pop_pred_mt$se.fit2)
uci<-10^cumsum(pop_pred_mt$fit2 + 1.96*pop_pred_mt$se.fit2)

dft<-data.frame(x,y, lci, uci)

library(ggplot2)
ggplot(dft, aes(x,y))+
  geom_line(size=2)+
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.3, colour=NA)+
  labs(y = "Mammal Population Index", x = "")+
  scale_y_continuous(breaks=seq(0, 2, 0.25))+
  theme_bw()+
  theme(text = element_text(size=20), axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")))
###########################################







plot(y=10^cumsum(pop_pred_m$fit2), x=2006:2099, ylim=c(0, 3), type="l", ylab="Mammal Population Index", xlab="Year", main="Predicted Mammal Population Trends Under RCP 8.5")
lines(y=10^cumsum(pop_pred_m$fit2 - 1.96*pop_pred_m$se.fit2), x=2006:2099, lty=2, col="Red")
lines(y=10^cumsum(pop_pred_m$fit2 + 1.96*pop_pred_m$se.fit2), x=2006:2099, lty=2, col="Red")

x = c(1:149)

fun<-function(y){
  
  z<-matrix(y)
  if (sum(is.na(z))<1){
    g<-mgcv:::gam(matrix(z) ~ s(x, k=8), fx=TRUE)
    g_pred<-predict(g)
    g_pred<-g_pred[55:148]
    #g_pred<-g_pred[55:100] #2050
    d<-diff(g_pred)
    gam_diff<-scale(d, center = centre_tempm, scale = scale_tempm)
    ave_mod<-data.frame(gam_diff, rep(0, length(gam_diff)), rep(0, length(gam_diff)))
    colnames(ave_mod)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
    pop_pred<-predict(m1_m, ave_mod, re.form=NA, se.fit=T)
  } else{
    
    #pop_pred$fit<-rep(NA,45) #2050
    pop_pred$fit<-rep(NA,93)
  }
  return(pop_pred$fit)
}

#gam_pred_rast_2050<-calc(obs_mod_ann, fun)
#gam_pred_rast<-calc(obs_mod_ann, fun)

gam_pred_rast_tsl<-calc(obs_mod_ann_tsl, fun)

plot(gam_pred_rast_tsl[[93]])

gam_pred_rast_2100<-calc(obs_mod_ann, fun)


#plot(gam_pred_rast_2100[[93]])

pop_index_2100<-function(x){10^cumsum(x)}


index_2100_tsl<-calc(gam_pred_rast_tsl, pop_index_2100)

#index_2100<-index_2100_tsl[[93]]

plot(index_2100_tsl[[93]])



# pop_2100<-function(x){10^sum(x)}
# 
# 
# 
# pred_2100<-calc(gam_pred_rast_2100, pop_2100)
# 
# pred_2050<-calc(gam_pred_rast_2050, pop_2100)

library(RColorBrewer)
library(rgdal)
library(ggplot2)
library(rgeos)
world<-readOGR(dsn=getwd(), layer="ne_10m_land")
e<-c(-180, 180, -63, 83.6341)
world_c<-crop(world, e)
plot(world_c)
world_df<-fortify(world_c)
colnames(world_df)[c(1:2)]<-c("Longitude", "Latitude")


pred_2100_df<-rasterToPoints(index_2100_tsl[[93]])
#pred_2100_df<-rasterToPoints(log10(pred_2100))
pred_2100_df2<-data.frame(pred_2100_df)
colnames(pred_2100_df2)<-c("Longitude", "Latitude", "Population_Decline")


theme_opts<-list(theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_rect(fill = 'white', colour = NA),
                       plot.background = element_rect(),
                       axis.line = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank()))

ggplot(data=pred_2100_df2, aes(Longitude, Latitude))+
   geom_raster(aes(fill=Population_Decline))+
  scale_fill_gradientn(colors=brewer.pal(7, "RdYlGn"), guide=guide_colourbar(title=""))+ 
  #labels=c("-99%", "-50%", "0%", "+50%", "+100%"), breaks=c(0.01,0.5, 1, 1.5, 2))+   
  #labels=c("-99.9999%", "-99.99%", "-99%","0%"))+
  geom_path(data=world_df, aes(Longitude, Latitude, group=group))+
  coord_fixed(ratio = 1)+
  theme_bw()+
  labs(title="Mammal Population Index 2100 - (2005 = 1)")+
  theme(text = element_text(size=20))
  




#   brk <- c(1, 0, -1, -2, -3, -4)
# plot(log10(pred_2100), col=brewer.pal(7, "RdYlGn"))
# 
# #plot(world, lwd=0.1, add=T)
# 
# legend("right",inset=F, title="Predicted Population Decline",
#        c("0%","90%","99%", "99.9%", "99.99%", "99.999%", "99.9999%"), fill=terrain.colors(7))
# 
# 
# 






plot(log10(pred_2050))
hist(pred_2050)



plot((pred_2100 - 1)*100,col=terrain.colors(12), breaks=brk)



#install.packages("merTools")
library(merTools)
preds <- predictInterval(m1c, newdata = ave_mod, n.sims = 999)



