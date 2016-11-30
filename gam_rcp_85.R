  CR40s<-brick("cru_ts3.23.1941.1950.tmp.dat.nc")
  CR50s<-brick("cru_ts3.23.1951.1960.tmp.dat.nc")
  CR60s<-brick("cru_ts3.23.1961.1970.tmp.dat.nc")
  CR70s<-brick("cru_ts3.23.1971.1980.tmp.dat.nc")
  CR80s<-brick("cru_ts3.23.1981.1990.tmp.dat.nc")
  CR90s<-brick("cru_ts3.23.1991.2000.tmp.dat.nc")
  CR00s<-brick("cru_ts3.23.2001.2010.tmp.dat.nc")
  
  obs<-stack(CR40s[[109:120]], CR50s, CR60s, CR70s, CR80s, CR90s,CR00s[[1:48]] )
  
  
  ###85 rcp
  
  ts85_05_30<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_200512-203011.nc")
  ts85_30_55<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_203012-205511.nc")
  ts85_55_80<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_205512-208011.nc")
  ts85_80_00<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/HadGEM2_ES_RCP26/ts_Amon_HadGEM2-ES_esmrcp85_r1i1p1_208012-210011.nc")
  
  ts85_05_00<-stack(ts85_05_30[[2:300]],ts85_30_55,ts85_55_80,ts85_80_00[[1:229]])
  #ts_05_00<-brick(ts_05_00)
  
  modelcc<-rotate(ts85_05_00)
  
  obs_r<-resample(obs, modelcc)
  
  obs_rk<-obs_r+273.15
  
  mask_mod<-mask(modelcc, obs_rk[[1]])
  
  #obs_mod<-stack(obs_rk, modelcc)
  
  obs_mod<-stack(obs_rk, mask_mod)
  
  obs_mod_mon_mean<-matrix(cellStats(obs_mod, mean))
  
  plot(obs_mod_mon_mean)
  ##################################################
  year_mon<-rep(1:(nlayers(obs_mod)/12), each=12)
  
  #changing from monthly to annual data
  obs_mod_ann<-stackApply(obs_mod, indices=year_mon, fun=mean, na.rm=TRUE)
  
  obs_mat_ann<-matrix(cellStats(obs_mod_ann, mean))
  plot(obs_mat_ann)
  
  x = c(1:149)
  g = mgcv:::gam(obs_mat_ann~s(x, k=8), fx=TRUE)
  plot(g)
  
  g_pred = predict(g)
  g_pred = g_pred[55:148]
d = diff(g_pred)




centre_temp<-0.04526127
scale_temp<-0.07075669

gam_diff<-scale(d, center=centre_temp, scale=scale_temp)

ave_mod<-data.frame(gam_diff, rep(0, length(gam_diff)), rep(0, length(gam_diff)))
colnames(ave_mod)<-c("mean_slope_scale", "change_rate_scale", "Bodymass_scale")
pop_pred<-predict(m1c, ave_mod, re.form=NA, se.fit=T)

plot(pop_pred$fit, ylim=c(-0.04, 0.01))
lines(pop_pred$fit + 1.96*pop_pred$se.fit)
lines(pop_pred$fit - 1.96*pop_pred$se.fit)
10^sum(pop_pred$fit + 1.96*pop_pred$se.fit)
10^sum(pop_pred$fit - 1.96*pop_pred$se.fit)
10^sum(pop_pred$fit)

plot(10^cumsum(pop_pred$fit), ylim=c(0, 1.5))
lines(10^cumsum(pop_pred$fit - 1.96*pop_pred$se.fit))
lines(10^cumsum(pop_pred$fit + 1.96*pop_pred$se.fit))



#install.packages("merTools")
library(merTools)
preds <- predictInterval(m1c, newdata = ave_mod, n.sims = 999)



