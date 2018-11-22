

#There are some NAs in rep_id - not sure why..... find out!
#need to take out the population trends which crashed out and make a note of which ones these were
#need to talk to damaris about optimising the matrix

library(reshape2)
library(sp)
library(raster)
wd<-getwd()
lpi<-read.csv("LPI_pops_20160523_edited.csv")

spin_years<-1940:1949
years<-1950:2005

binomial = "Ursus_arctos"
#demoniche_folder<-"C:/Users/Fiona/Documents/PhD/PhD_Method/Legion/cervus_output"
#demoniche_folder<-"D:/Fiona/Git_Method/Git_Method/Legion/snow_cervus_bias/output_new"

demoniche_folder<-paste(wd, "/Legion/snow_bear_bias_faster/output_smooth", sep="")
#demoniche_folder<-paste(wd, "/Legion/snow_cervus_test/output_test", sep="")
l<-list.files(demoniche_folder)
nf<-length(list.files(paste(demoniche_folder, l[1], sep="/")))



highfoldernames<-list.files(demoniche_folder)
#lowfoldernames<-rep(1:nf, each=length(l))

#foldernames<-paste(highfoldernames, lowfoldernames, sep="/")

sp_lpi<-lpi[lpi$Binomial == binomial & lpi$Specific_location ==1 & lpi$Region == "North America",]


xy<-data.frame(sp_lpi$Longitude, sp_lpi$Latitude)
coordinates(xy)<-c("sp_lpi.Longitude", "sp_lpi.Latitude")
proj4string(xy) <- CRS("+init=epsg:4326") # WGS 84
#CRS.new <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
#xy <- spTransform(xy, CRS.new)

convert_pop_out<-function(foldername){
  
  if(length(dir(paste(demoniche_folder, foldername, sep="/"))) !=0){
    pop_out<-read.csv(paste(demoniche_folder ,foldername, "Reference_matrix_pop_output.csv", sep="/"), header = TRUE)
    pop_out<-pop_out[,-1]
    coordinates(pop_out) <- ~ X + Y
    gridded(pop_out) <- TRUE
    rasterDF <- raster:::stack(pop_out)
    crs(rasterDF)<- CRS("+proj=longlat +datum=WGS84")
    trends<-raster:::extract(rasterDF,xy, buffer = 50000, fun = mean, na.rm = TRUE)
    # SD<-strsplit(foldername, "[/_]")[[1]][2]
    # sdd<-strsplit(foldername, "[/_]")[[1]][3]
    # ldd<-strsplit(foldername, "[/_]")[[1]][4]
    # dens<-strsplit(foldername, "[/_]")[[1]][5]
    # link<-strsplit(foldername, "[/_]")[[1]][6]
    # med_disp<-strsplit(foldername, "[/_]")[[1]][7]
    # max_disp<-strsplit(foldername, "[/_]")[[1]][8]
    # rep_id<-strsplit(foldername, "[/_]")[[1]][9]
    SD<-strsplit(foldername, "[/_]")[[1]][6]
    sdd<-strsplit(foldername, "[/_]")[[1]][7]
    ldd<-strsplit(foldername, "[/_]")[[1]][8]
    dens<-strsplit(foldername, "[/_]")[[1]][9]
    link<-strsplit(foldername, "[/_]")[[1]][10]
    med_disp<-strsplit(foldername, "[/_]")[[1]][11]
    max_disp<-strsplit(foldername, "[/_]")[[1]][12]
    rep_id<-strsplit(foldername, "[/_]")[[1]][13]
    trends_df<-data.frame(sp_lpi$ID,med_disp,sdd,ldd,SD,dens,link,rep_id,trends)
  }
}
demoniche_pop_out<-lapply(highfoldernames, convert_pop_out)
df <- do.call("rbind", demoniche_pop_out)

#write.csv(df, "ursus_arctos_results_all.csv")
#save(df, file ="ursus_arctos_results_all.RData")

df<-read.csv("ursus_arctos_results_all.csv")
df<-df[,-1]

dfm<-as.matrix(df)

# lambda<-function(x){
#   
#   l10<-diff(log1p(as.numeric(x[10:length(x)])))
#   #l10<-10^diff(log1p(as.numeric(x[20:length(x)])))
#   
# }
# 
# dft<-t(apply(dfm,1,lambda))
# 
# df_lambda<-data.frame(dfm[,1:8],dft)
# 
# colnames(df_lambda)[9:ncol(df_lambda)]<-colnames(dfm)[9:ncol(df_lambda)]


melt_df<-melt(df, id=1:8)
melt_df$year<-as.numeric(gsub("Year_", "", melt_df$variable))

melt_short<-melt_df[melt_df$year>spin_years[length(spin_years)] ,]
melt_short$sp_lpi.ID<-as.factor(melt_short$sp_lpi.ID)


library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=interaction(ldd, SD, rep_id), colour= sp_lpi.ID))+
  geom_line()+
  #geom_smooth()+
  facet_grid(~ sp_lpi.ID)


ggplot(melt_short, aes(x= year, y=c(diff(log10(value)),NA), group=interaction(ldd, SD, rep_id), colour= sp_lpi.ID))+
  geom_line()+
  #geom_smooth()+
  facet_grid(~ sp_lpi.ID)

ggplot(melt_short, aes(x= year, y=c(diff(log10(value)),NA), group=interaction(ldd, SD), colour= sp_lpi.ID))+
  #geom_line()+
  geom_smooth()+
  facet_grid(~ sp_lpi.ID)


melt_short$ldd_sd<-paste(melt_short$ldd, melt_short$SD, sep="_")

dfl<-split(melt_short, list(melt_short$sp_lpi.ID, melt_short$ldd_sd, melt_short$rep_id))
gam_smooth<-function(x){
  
  ind<-seq(length(unique(x$rep_id)),(length(unique(x$rep_id)) * length(unique(x$year))* length(unique(x$sp_lpi.ID))), by=length(unique(x$rep_id)))
  id<-as.numeric(as.character(x$sp_lpi.ID[1]))
  ldd<-as.numeric(as.character(x$ldd[1]))
  sd<-as.numeric(as.character(x$SD[1]))
  rep_id<-as.numeric(as.character(x$rep_id[1]))
  
  if(nrow(x)>0){
    mg<-mgcv:::gam(log10(value+1)~s(year, k = round(nrow(x)/2)), data = x)
    fmg<-fitted.values(mg)
    smooth_lambdas<-diff(fmg[ind])
    
    lambdas<-diff(log10(x$value +1))
    out<-cbind(id,ldd,sd,rep_id,lambdas,smooth_lambdas, unique(x$year)[-length(unique(x$year))])
    gam_out<-data.frame(out)
    colnames(gam_out)<-c("ID","ldd","SD","rep_id","Lambdas","Smooth_Lambdas" ,"Year")
    print(paste(id, ldd, sd, rep_id))
    
  } else {
    
    out<-cbind(NA,NA,NA,NA,NA,NA,NA)
    gam_out<-data.frame(out)
    colnames(gam_out)<-c("ID","ldd","SD","rep_id","Lambdas", "Smooth_Lambdas","Year")
    print(paste(id, ldd, sd, rep_id))
    
  }
  
  return(gam_out)
  
}

df_out<-lapply(dfl, gam_smooth)
melt_lambda_short<-do.call( "rbind", df_out)

melt_lambda_short$ID<-as.factor(melt_lambda_short$ID)

melt_lambda_short<-melt_lambda_short[complete.cases(melt_lambda_short),]

#write.csv(melt_lambda_short, "melt_lambda_short_ursus.csv")

melt_lambda_short<-read.csv("melt_lambda_short_ursus.csv")


ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x= Year, y=Lambdas, group=interaction(ldd, SD), colour= ID), method="loess", SE=FALSE)+
  facet_grid(~ ID)

library(zoo)
library(taRifx)
library(plyr)
library(mgcv)


pops<-sp_lpi[,c(1,65:120)]

pop_counts<-sp_lpi[,c(65:120)]
pop_counts <- (pop_counts !="NULL")
points_per_pop1950_2005 = rowSums(pop_counts)
#if(sum(points_per_pop1950_2005<5)>=1){
# pyr<-pyr[-which(points_per_pop1950_2005<5),]
# pops<-pops[-which(points_per_pop1950_2005<5),]
#}

colnames(pops)<-c("ID", 1950:2005)

colnames(pops)[2:ncol(pops)]<-paste("Year", 1950:2005, sep="_")
pops[pops=="NULL"]<-NA
pops$rep_id<-"Observed"
pops$md_id<-"Observed"

popsm<-as.matrix(pops)

gam_lpi<-function(x){
  #subsetting the population data by each population
  spid = x[2:(length(x)-2)]                     #subsetting only the dates
  names(spid)<-1950:2005              #renaming the date column names as R doesn't like numbered column names
  spid<-as.numeric(spid)
  pop_datab <- (!is.na(spid) )
  points = sum(pop_datab)
  id<-x[1]
  Date<-1950:2005
  spidt<-destring(t(spid))
 
   if (points > 1){
    time<-length(min(which(!is.na(spidt))):max(which(!is.na(spidt))))
    missing<-time-points
  
    Year<-Date[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
    Population<-spidt[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
    Population[Population == 0] <- mean(Population, na.rm=TRUE)*0.01 #if a population is zero one year thhis is replaced with 1% of the average population estimate - because you can log zeros
  
    df<-data.frame(Year,Population)
  
  #not sure what this does - adding a constant of 1 so that logging doesn't go weird?
    if (sum(na.omit(df$Population<1))>0) {
      df$Population<-df$Population+1
    }
  
    PopN = log10(df$Population)
  

    
    if (points >=6) {
    
      if (length(na.omit(PopN)) >=6) {
      SmoothParm = round(length(na.omit(PopN))/2)
        } else {
      SmoothParm=3
      }
    
    mg2<-mgcv:::gam(PopN ~ s(Year, k=SmoothParm), fx=TRUE)
    pv2 <- predict(mg2,df,type="response",se=TRUE)
    R_sq2<-summary(mg2)$r.sq
    model<-1
    pv2$fit[pv2$fit <= 0] <- NA
    
    lambda2<-diff(pv2$fit)
    
    ial<-data.frame(id, Year[-length(Year)], lambda2)
    colnames(ial)<-c("ID", "Year", "Lambdas")
  } else {
    
    lint<-na.approx(PopN)
    lint[lint<=0] <- NA
    lambda2<-diff(lint)
    ial<-data.frame(id, Year[-length(Year)], lambda2)
    colnames(ial)<-c("ID", "Year", "Lambdas")
  }
    } else{
    ial<-data.frame(id, NA,NA)
    colnames(ial)<-c("ID", "Year", "Lambdas")
  }
  print(ial)
  return(ial)
}

gam_lpi_r<-apply(popsm,  1, gam_lpi)
gam_r<-do.call( "rbind", gam_lpi_r)

gam_r<-gam_r[gam_r$Year <2005,]

fill<-data.frame(rep(pops$ID, each=length(1950:2004)), 1950:2004)
colnames(fill)<-c("ID", "Year")

all_year_ab<-join(fill, gam_r, type="right")

#write.csv(all_year_ab, "all_year_ab_ursus.csv")

colnames(melt_lambda_short)[1]<-"ID"

all_year_ab$ID<-as.numeric(as.character(all_year_ab$ID))
melt_lambda_short$ID<-as.numeric(as.character(melt_lambda_short$ID))

#melt_lambda_short<-melt_lambda_short[melt_lambda_short$ID %in% all_year_ab$ID,]


ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x = Year, y= Lambdas, group=interaction(ldd, SD)), colour = "black", alpha = 0.3, method = "loess")+
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red", method = "loess")+
  #geom_line(data =  gam_r_lambda, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ID)



smooth_vals_pred = predict(loess(Lambdas~Year,melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[1],]), melt_lambda_short$year[melt_lambda_short$ID == melt_lambda_short$ID[1]])

smooth_vals_obs = predict(loess(Lambdas~Year,all_year_ab[all_year_ab$ID == all_year_ab$ID[1],]), all_year_ab$Year[all_year_ab$ID == all_year_ab$ID[1]])

####need to match the years up 
#rmse(smooth_vals_pred, smooth_vals_obs)



###sdm trends

species_directory<-paste(wd, "/",binomial, "_bias", sep="")
sdm_folder<-paste(species_directory, "SDM_folder_bias", sep = "/")

sdm_folder<-paste(wd, "/Legion/snow_bear_bias_faster/", sep="")


sdm_stack<-stack(paste(sdm_folder, "/hyde_weighted_ensemble_sdm_", years,".tif", sep=""))
patch_stack<-stack(paste(sdm_folder, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years,".tif", sep=""))

pyr<-subset(lpi, Binomial ==binomial & Specific_location==1& (Region!="Europe"& Region!="Asia"))    #record 11470 had wrong longitude - in Russia!
pops<-pyr[,c(65:120)]
pop_counts <- (pops !="NULL")
points_per_pop1950_2005 = rowSums(pop_counts)
#if(sum(points_per_pop1950_2005<5)>=1){
#  pyr<-pyr[-which(points_per_pop1950_2005<5),]
#  pops<-pops[-which(points_per_pop1950_2005<5),]
#}
pops<-data.frame(pyr$ID, pops)
colnames(pops)<-c("ID", 1950:2005)


#formatting the lpi data for use in demoniche
pyrs<-pyr[,c("ID","Binomial","Longitude","Latitude")]
xy_lpi<-data.frame(pyrs$Longitude, pyrs$Latitude)

sdm_lpi<-raster:::extract(sdm_stack,xy_lpi, buffer = 50000, fun = mean, na.rm = TRUE)

colnames(sdm_lpi)<-1950:2005
sdm_lpi_melt<-melt(sdm_lpi)
colnames(sdm_lpi_melt)<-c("ID", "Year", "HSI")
sdm_lpi_melt$ID<-rep(pyrs$ID, length(1950:2005))


ggplot(sdm_lpi_melt, aes(x = Year, y= HSI, group = ID))+
  geom_smooth()+
  facet_grid(.~ID)


sdm_lpi_melt %>%
  group_by(ID)%>%
  mutate(log_hsi = log10(HSI))%>%
  mutate(HSI_Lambdas = c(diff(log_hsi), NA)) %>%
  filter(ID == 3478)

  




sdl<-split(sdm_lpi_melt, list(sdm_lpi_melt$ID))
sdm_smooth<-function(x){
  
  #ind<-seq(length(unique(x$rep_id)),(length(unique(x$rep_id)) * length(unique(x$year))* length(unique(x$sp_lpi.ID))), by=length(unique(x$rep_id)))
  id<-as.numeric(as.character(x$ID[1]))
  mg<-gam(log10(HSI)~s(Year, bs="cs", k = 28), data = x)
 # mg<-gam(log10(HSI)~s(Year, bs="cs", k = -1), data = x)
  fmg<-fitted.values(mg)
  smooth_lambdas<-diff(fmg)

  lambdas<-diff(log10(x$HSI))
  out<-cbind(id,lambdas,smooth_lambdas, unique(x$Year)[-length(unique(x$Year))])
  gam_out<-data.frame(out)
  colnames(gam_out)<-c("ID","HSI_Lambdas","HSI_Lambdas_Smooth", "Year")
  return(gam_out)
}

df_out<-lapply(sdl, sdm_smooth)
sdm_lambdas_melt<-do.call( "rbind", df_out)



#with loess
ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x = Year, y= Lambdas, group=interaction(ldd, SD)), colour = "black", fill="grey", method = "loess")+   #demoniche
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red", fill="lightcoral",method = "loess")+    #real
  geom_smooth(data = sdm_lambdas_melt, aes(x = Year, y= HSI_Lambdas, group=ID), colour = "blue", fill="light blue",method = "loess")+   #sdm trend
  facet_grid(.~ID,labeller=label_both)

#with gam
ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x = Year, y= Lambdas, group=interaction(ldd, SD)), colour = "black", fill="grey", method = "gam", formula = y ~ s(x, bs = "cs", k = -1))+   #demoniche
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red", fill="lightcoral",method = "gam", formula = y ~ s(x, bs = "cs", k =-1))+    #real
  geom_smooth(data = sdm_lambdas_melt, aes(x = Year, y= HSI_Lambdas, group=ID), colour = "blue", fill="light blue",method = "gam", formula = y ~ s(x, bs = "cs", k = -1))+   #sdm trend
  facet_grid(.~ID,labeller=label_both)


#need start year and end year


library(Metrics)

melt_lambda_short$ldd_sd<-paste(melt_lambda_short$ldd, melt_lambda_short$SD, sep="_")

vg<-expand.grid(unique(melt_lambda_short$ID), unique(melt_lambda_short$ldd_sd), unique(melt_lambda_short$rep_id))
vg$Var2<-as.character(vg$Var2)

vg$ldd<-matrix(as.numeric(unlist(strsplit(vg$Var2, "_"))), ncol = 2, byrow = T)[,1]
vg$SD<-matrix(as.numeric(unlist(strsplit(vg$Var2, "_"))), ncol = 2, byrow = T)[,2]
colnames(vg)<-c("ID", "ldd_SD", "rep_id", "ldd", "SD")
rmse_get_sdm<-function(x){
  
  cnd_x<-melt_lambda_short[melt_lambda_short$ID == vg[x,"ID"]& melt_lambda_short$ldd == vg[x,"ldd"] & melt_lambda_short$SD == vg[x,"SD"]& melt_lambda_short$rep_id == vg[x,"rep_id"],]
  obs_x<-all_year_ab[all_year_ab$ID == vg[x,"ID"],]
  sdm_x<-sdm_lambdas_melt[sdm_lambdas_melt$ID == vg[x,"ID"],]
  ID<-vg[x,"ID"]
  rep_id<-vg[x,"rep_id"]
  ldd<-vg[x,"ldd"]
  sd<-vg[x,"SD"]
  
  if (nrow(cnd_x)>1){
    
    #need to get separate line per iteration of parameters for cnd
    # cnd_gam = mgcv:::gam(Lambdas~s(Year, bs="cs", k = -1),data = cnd_x )
    # obs_gam = gam(Lambdas~s(Year, bs="cs", k = -1),data = obs_x)
    # sdm_gam = gam(HSI_Lambdas~s(Year, bs="cs", k =-1),data =sdm_x)
    # 
    # smooth_vals_cnd = predict(cnd_gam,newdata = cnd_x)
    # smooth_vals_obs = predict(obs_gam,newdata = obs_x)
    # smooth_vals_sdm = predict(sdm_gam,newdata = sdm_x)    

    start_cnd<-which(unique(cnd_x$Year) == min(obs_x$Year))
    end_cnd<-which(unique(cnd_x$Year) == max(obs_x$Year))
    
    start_obs<-which(unique(obs_x$Year) == min(obs_x$Year))
    end_obs<-which(unique(obs_x$Year) == max(obs_x$Year))
    
    start_sdm<-which(unique(sdm_x$Year) == min(obs_x$Year))
    end_sdm<-which(unique(sdm_x$Year) == max(obs_x$Year))
    
    cnd_rmse<-rmse(cnd_x$Lambdas[start_cnd:end_cnd], obs_x$Lambdas[start_obs:end_obs])
    sdm_rmse<-rmse(sdm_x$HSI_Lambdas[start_sdm:end_sdm], obs_x$Lambdas[start_obs:end_obs])
    
    cnd_sdm_rmse<-rmse(cnd_x$Lambdas[start_cnd:end_cnd], sdm_x$HSI_Lambdas[start_sdm:end_sdm])
    
    rmse_out<-data.frame(ID, rep_id, ldd, sd, cnd_rmse,sdm_rmse, cnd_sdm_rmse)
    colnames(rmse_out)<-c("ID","rep_id","ldd","sd" ,"cnd", "sdm", "cnd_sdm")
    # rmse_out<-data.frame(cnd_rmse)
    # colnames(rmse_out)<-c("cnd")
    print((x/nrow(vg))*100)
    return(rmse_out)
  }
}
rmse_scores<-lapply(1:nrow(vg), rmse_get_sdm)
rmse_out<-do.call( "rbind", rmse_scores)

rmse_out$diff<-rmse_out$cnd - rmse_out$sdm


#write.csv(rmse_out, "rmse_out_ursus_all.csv")

hist(rmse_out$diff, main = "RMSE difference between CND and SDM")

rmse_out$ID<-as.factor(rmse_out$ID)

ggplot(rmse_out,aes(x=ID, y=cnd, group = ID)) + 
  geom_violin(linetype = "solid", size = 1.5)+
  #geom_point(size = 2)+
  geom_point(data = rmse_out, aes(x = ID, y = sdm), size = 5, colour = "red")+
   theme_bw()+
  xlab("Population")+
  ylab("Root Mean Square Error")+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=20))


#vg<-vg[vg$ID == 4357 |vg$ID == 4358, ]



ccf_get<-function(x){
  
  cnd_x<-melt_lambda_short[melt_lambda_short$ID == vg[x,"ID"]& melt_lambda_short$ldd == vg[x,"ldd"]& melt_lambda_short$SD== vg[x,"SD"]& melt_lambda_short$rep_id == vg[x,"rep_id"],]
  obs_x<-all_year_ab[all_year_ab$ID == vg[x,"ID"],]
  obs_x<-obs_x[obs_x$Year >= 1951,]
  sdm_x<-sdm_lambdas_melt[sdm_lambdas_melt$ID == vg[x,"ID"],]
  
  if (nrow(cnd_x)>0){
    
    start_cnd<-which(unique(cnd_x$Year) == min(obs_x$Year))
    end_cnd<-which(unique(cnd_x$Year) == max(obs_x$Year))
    
    start_obs<-which(unique(obs_x$Year) == min(obs_x$Year))
    end_obs<-which(unique(obs_x$Year) == max(obs_x$Year))
    
    start_sdm<-which(unique(sdm_x$Year) == min(obs_x$Year))
    end_sdm<-which(unique(sdm_x$Year) == max(obs_x$Year))
    
    cnd_ccf<-ccf(cnd_x$Lambdas[start_cnd:end_cnd], obs_x$Lambdas[start_obs:end_obs], type="correlation")
    sdm_ccf<-ccf(sdm_x$HSI_Lambdas[start_sdm:end_sdm], obs_x$Lambdas[start_obs:end_obs], type="correlation")
    # cnd_ccf_smooth<-ccf(cnd_x$Smooth_Lambdas[start_cnd:end_cnd], obs_x$Lambdas[start_obs:end_obs], type="correlation")
    # sdm_ccf_smooth<-ccf(sdm_x$HSI_Lambdas_Smooth[start_sdm:end_sdm], obs_x$Lambdas[start_obs:end_obs], type="correlation")
    # 
    
    lag_n5_cnd<-cnd_ccf$acf[which(as.numeric(cnd_ccf$lag) == -5)]
    lag_n5_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -5)]
    # lag_n5_cnd_sm<-cnd_ccf_smooth$acf[which(as.numeric(cnd_ccf_smooth$lag) == -5)]
    # lag_n5_sdm_sm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -5)]


    lag_n4_cnd<-cnd_ccf$acf[which(as.numeric(cnd_ccf$lag) == -4)]
    lag_n4_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -4)]
    # lag_n4_cnd_sm<-cnd_ccf_smooth$acf[which(as.numeric(cnd_ccf_smooth$lag) == -4)]
    # lag_n4_sdm_sm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -4)]
    
    
    lag_n3_cnd<-cnd_ccf$acf[which(as.numeric(cnd_ccf$lag) == -3)]
    lag_n3_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -3)]
    # lag_n3_cnd_sm<-cnd_ccf_smooth$acf[which(as.numeric(cnd_ccf_smooth$lag) == -3)]
    # lag_n3_sdm_sm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -3)]
    
    
    lag_n2_cnd<-cnd_ccf$acf[which(as.numeric(cnd_ccf$lag) == -2)]
    lag_n2_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -2)]
    # lag_n2_cnd_sm<-cnd_ccf_smooth$acf[which(as.numeric(cnd_ccf_smooth$lag) == -2)]
    # lag_n2_sdm_sm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -2)]
    
    
    lag_n1_cnd<-cnd_ccf$acf[which(as.numeric(cnd_ccf$lag) == -1)]
    lag_n1_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == -1)]
    # lag_n1_cnd_sm<-cnd_ccf_smooth$acf[which(as.numeric(cnd_ccf_smooth$lag) == -1)]
    # lag_n1_sdm_sm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == -1)]
    

    lag_0_cnd<-cnd_ccf$acf[which(as.numeric(cnd_ccf$lag) == 0)]
    lag_0_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 0)]
    # lag_0_cnd_sm<-cnd_ccf_smooth$acf[which(as.numeric(cnd_ccf_smooth$lag) == 0)]
    # lag_0_sdm_sm<-sdm_ccf_smooth$acf[which(as.numeric(sdm_ccf_smooth$lag) == 0)]
    
    
    n_used<-cnd_ccf$n.used
    
    ccf_out<-data.frame(cnd_x$ID, cnd_x$ldd, cnd_x$SD, cnd_x$rep_id,n_used,lag_n5_cnd,lag_n5_sdm,lag_n4_cnd,lag_n4_sdm,lag_n3_cnd,lag_n3_sdm,lag_n2_cnd, lag_n2_sdm, lag_n1_cnd, lag_n1_sdm, lag_0_cnd, lag_0_sdm)
  #  ccf_out<-data.frame(cnd_x$ID, cnd_x$ldd, cnd_x$SD, cnd_x$rep_id,n_used,lag_n3_cnd,lag_n3_sdm,lag_n2_cnd, lag_n2_sdm, lag_n1_cnd, lag_n1_sdm, lag_0_cnd, lag_0_sdm)

    # colnames(ccf_out)<-c("ID","ldd","SD","rep_id","cnd_n5", "sdm_n5","cnd_n4", "sdm_n4","cnd_n3", "sdm_n3","cnd_n2", "sdm_n2","cnd_n1", "sdm_n1", "cnd_0", "sdm_0","cnd_1", "sdm_1","cnd_2", "sdm_2","cnd_3", "sdm_3","cnd_4", "sdm_4","cnd_5", "sdm_5")
    colnames(ccf_out)<-c("ID","ldd","SD","rep_id","N_used", "cnd_n5", "sdm_n5","cnd_n4", "sdm_n4","cnd_n3", "sdm_n3","cnd_n2", "sdm_n2","cnd_n1", "sdm_n1", "cnd_0", "sdm_0")
   # colnames(ccf_out)<-c("ID","ldd","SD","rep_id","N_used", "cnd_n3", "sdm_n3","cnd_n2", "sdm_n2","cnd_n1", "sdm_n1", "cnd_0", "sdm_0")
    
    #  ccf_out<-data.frame(cnd_x$ID, cnd_x$ldd, cnd_x$SD, cnd_x$rep_id, lag_n3_cnd,lag_n3_sdm,lag_n2_cnd, lag_n2_sdm, lag_n1_cnd, lag_n1_sdm, lag_0_cnd, lag_0_sdm,lag_1_cnd, lag_1_sdm,lag_2_cnd, lag_2_sdm,lag_3_cnd, lag_3_sdm)
    #  colnames(ccf_out)<-c("ID","ldd","SD","rep_id","cnd_n3", "sdm_n3","cnd_n2", "sdm_n2","cnd_n1", "sdm_n1", "cnd_0", "sdm_0","cnd_1", "sdm_1","cnd_2", "sdm_2","cnd_3", "sdm_3")
    
    print((x/nrow(vg))*100)
    return(ccf_out)
    
  }
}


ccf_scores<-lapply(1:nrow(vg), ccf_get)
#ccf_scores<-sapply(1:10000, ccf_get)

ccf_all<-ccf_scores[lapply(ccf_scores,length)>0] 

ccf_df<-do.call( "rbind",ccf_all)
#colnames(ccf_df)<-c("ID","ldd","SD","cnd lag -5", "sdm lag -5","cnd lag -4", "sdm lag -4","cnd lag -3", "sdm lag -3","cnd lag -2", "sdm lag -2","cnd lag -1", "sdm lag -1", "cnd lag 0", "sdm lag 0","cnd lag 1", "sdm lag 1","cnd lag 2", "sdm lag 2","cnd lag 3", "sdm lag 3","cnd lag 4", "sdm lag 4","cnd lag 5", "sdm lag 5")
#colnames(ccf_df)<-c("ID","ldd","SD","rep_id","cnd lag -3", "sdm lag -3","cnd lag -2", "sdm lag -2","cnd lag -1", "sdm lag -1", "cnd lag 0", "sdm lag 0","cnd lag 1", "sdm lag 1","cnd lag 2", "sdm lag 2","cnd lag 3", "sdm lag 3")

colnames(ccf_df)<-c("ID","ldd","SD","rep_id","N_used","cnd lag -5", "sdm lag -5","cnd lag -4", "sdm lag -4","cnd lag -3", "sdm lag -3","cnd lag -2", "sdm lag -2","cnd lag -1", "sdm lag -1", "cnd lag 0", "sdm lag 0")
#colnames(ccf_df)<-c("ID","ldd","SD","rep_id","N_used","cnd lag -3", "sdm lag -3","cnd lag -2", "sdm lag -2","cnd lag -1", "sdm lag -1", "cnd lag 0", "sdm lag 0")


ccf_melt<-melt(ccf_df, id.vars = c("ID", "ldd", "SD", "rep_id", "N_used"))
#ccf_melt<-melt(ccf_df, id.vars = c("ID", "ldd", "SD"))
ccf_melt<-unique(ccf_melt)
ccf_melt_short<-ccf_melt
#

#write.table(ccf_melt, "ccf_melt_ursus2.csv")
#write.table(ccf_melt, "ccf_melt_ursus2_short.csv")
# #save(ccf_melt, file = "ccf_melt_ursus2.RData")
ccf_melt_short<-read.table("ccf_melt_ursus2_short.csv")
ccf_melt<-read.table("ccf_melt_ursus2.csv")
 #load("ccf_melt_ursus.RData")

ccf_melt<-rbind(ccf_melt, ccf_melt_short)


# for(i in 1:length(ccf_melt$variable)){
#   
#   if (grepl( "cnd",ccf_melt$variable[i])){
#     ccf_melt$model[i]<-"cnd"  
#   } 
#   
#   if (grepl( "sdm",ccf_melt$variable[i])){
#     ccf_melt$model[i]<-"sdm"  
#   } 
#   print(i) 
# }
# 
# ccf_melt_cnd<-ccf_melt[ccf_melt$model =="cnd",]
#ccf_melt_cnd$ID<-as.factor(ccf_melt_cnd$ID)


#ccf_melt_sdm<-ccf_melt[ccf_melt$model =="sdm",]
#ccf_melt_sdm$ID<-as.factor(ccf_melt_sdm$ID)


# plot_func<-function(x){
#   ggplot(ccf_melt_cnd[ccf_melt_cnd$ID == x,], aes(x=variable, y=value)) + 
#     geom_violin()+
#     geom_jitter(data = ccf_melt_cnd[ccf_melt_cnd$ID == x,],aes(x=variable, y=value,col = ID) )
# }
# 
# lapply(unique(ccf_melt_cnd$ID), plot_func)

ccf_melt$ldd<-as.factor(ccf_melt$ldd)
ccf_melt$SD<-as.factor(ccf_melt$SD)
ccf_melt$ID<-as.factor(ccf_melt$ID)

ccf_melt$model<-NA

ccf_melt$model<-ifelse(grepl("cnd", ccf_melt$variable), "cnd", "sdm")


ggplot(ccf_melt[grepl("cnd",ccf_melt$variable) ,],aes(x=ID, y=value, group = ID, colour = ID,  shape = ldd)) + 
  #geom_violin()+
  geom_jitter(data = ccf_melt[ccf_melt$model == "cnd",],aes(x=variable, y=value,col = ID) )+
  scale_x_discrete(limits=c("cnd lag -5","cnd lag -4","cnd lag -3","cnd lag -2","cnd lag -1","cnd lag 0"))



#lag0

ccf_melt_cnd0<-ccf_melt[ccf_melt$variable == "cnd lag 0" & ccf_melt$model == "cnd",]
ccf_melt_cnd0$ldd<-as.factor(ccf_melt_cnd0$ldd)
ccf_melt_cnd0$SD<-as.factor(ccf_melt_cnd0$SD)

ggplot(ccf_melt[ccf_melt$variable == "cnd lag 0" ,], aes(x=ID, y=value, group = ID, colour = SD,  shape = ldd)) + 
  #  geom_violin()+
  geom_jitter()+
  labs(y = "Correlation")+
  theme_bw()


ggplot(ccf_melt[ccf_melt$variable == "cnd lag 0",], aes(x=ID, y=value, group = ID)) + 
  #  geom_violin()+
  geom_jitter()+
  labs(y = "Correlation")+
  theme_bw()

# ccf_melt_sdm0<-ccf_melt_sdm[ccf_melt_sdm$variable == "sdm lag 0",]
# ccf_melt_sdm0$ldd<-as.factor(ccf_melt_sdm0$ldd)
# ccf_melt_sdm0$SD<-as.factor(ccf_melt_sdm0$SD)


ccf_melt$ID<-as.factor(ccf_melt$ID)


ggplot(ccf_melt[ccf_melt$variable == "cnd lag 0" ,],aes(x=ID, y=value, group = ID)) + 
  geom_boxplot(linetype = "solid", size = 1)+
  #geom_point(size = 2)+
  geom_point(data = ccf_melt[ccf_melt$variable == "sdm lag 0",], aes(x = ID, y = value), size = 5, colour = "red")+
  geom_hline(yintercept=0, linetype = "dashed")+
  theme_bw()+
  ylim(-1,1)+
  xlab("Population")+
  ylab("Correlation Coefficient")+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=20))



get_lag<-function(x){
  lag_out<-as.numeric(strsplit(as.character(ccf_melt$variable[x]), " ")[[1]][3])
  return(lag_out)
}

lags<-lapply(1:length(ccf_melt$variable), get_lag)
ccf_melt$lag<-do.call("rbind", lags)




library(dplyr)

ccf_max<-ccf_melt%>%
  group_by(ID,rep_id ,model)%>%
   mutate(max_value = max(value))%>%
  dplyr:::select(ID,rep_id,N_used,lag,model,value,max_value)%>%
  distinct()

# 
# ccf_max_cnd[ccf_max_cnd$lag == 0 & ccf_max_cnd$model == "sdm",]
# 
# 
# 
ccf_max<-ccf_max[ccf_max$value == ccf_max$max_value,]
# 
# #filter not working well
# 
# 
# ccf_max_cnd<-ccf_max[ccf_max$model == "cnd",]
# 
# 
# ccf_max_sdm<-ccf_max[ccf_max$model == "sdm",]
# ccf_max_sdm<-ccf_max_sdm[,c(1,3,5:7)]
# ccf_max_sdm<-unique(ccf_max_sdm)

ggplot(ccf_max[ccf_max$model == "cnd",],aes(x=ID, y=max_value, group = ID)) + 
  geom_boxplot(size = 1)+
  #geom_point(size = 2)+
  geom_point(data = ccf_max[ccf_max$model == "sdm",], aes(x = ID, y = value), size = 5, colour = "red")+
  geom_hline(yintercept=0, linetype = "dashed")+
  theme_bw()+
  ylim(-1,1)+
  xlab("Population")+
  ylab("Correlation Coefficient")+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=20))



# ggplot(ccf_melt[ccf_melt$variable == "cnd lag 0" ,],aes(x=ID, y=value, group = ID)) + 
#   geom_boxplot(size = 1)+
#   #geom_point(size = 2)+
#   geom_point(data = ccf_max[ccf_max$model == "sdm",], aes(x = ID, y = value), size = 5, colour = "red")+
#   geom_hline(yintercept=0, linetype = "dashed")+
#   theme_bw()+
#   ylim(-1,1)+
#   xlab("Population")+
#   ylab("Correlation Coefficient")+
#   theme(legend.position="none",axis.text=element_text(size=16),
#         axis.title=element_text(size=20))
# 




ggplot(ccf_melt_cnd0, aes(x=ID, y=value, group = ID, colour = ldd))+#,  shape = ldd)) + 
  #  geom_violin()+
  geom_jitter(alpha = 0.7)+
  labs(y = "Correlation")+
  theme_bw()


ccf_melt_cnd0[which(ccf_melt_cnd0$value < -0.9),]


#lag3
ccf_melt_cnd3<-ccf_melt_cnd[ccf_melt_cnd$variable == "cnd lag 3",]

ccf_melt_cnd3$ldd<-as.factor(ccf_melt_cnd3$ldd)
ccf_melt_cnd3$SD<-as.factor(ccf_melt_cnd3$SD)

ggplot(ccf_melt_cnd3, aes(x=ID, y=value, group = ID))+#, colour = SD,  shape = ldd)) + 
  #  geom_violin()+
  geom_jitter(alpha = 0.7)+
  labs(y = "Correlation")+
  theme_bw()



#lag1
ccf_melt_cnd1<-ccf_melt_cnd[ccf_melt_cnd$variable == "cnd lag 1",]

ccf_melt_cnd1$ldd<-as.factor(ccf_melt_cnd1$ldd)
ccf_melt_cnd1$SD<-as.factor(ccf_melt_cnd1$SD)

ggplot(ccf_melt_cnd1, aes(x=ID, y=value, group = ID))+#, colour = SD,  shape = ldd)) + 
  #  geom_violin()+
  geom_jitter(alpha = 0.1)+
  labs(y = "Correlation")+
  theme_bw()



#lag n3
ccf_melt_cndn3<-ccf_melt_cnd[ccf_melt_cnd$variable == "cnd lag -3",]

ccf_melt_cndn3$ldd<-as.factor(ccf_melt_cndn3$ldd)
ccf_melt_cndn3$SD<-as.factor(ccf_melt_cndn3$SD)

ggplot(ccf_melt_cndn3, aes(x=ID, y=value, group = ID))+#, colour = SD,  shape = ldd)) + 
  #  geom_violin()+
  geom_jitter(alpha = 0.1)+
  labs(y = "Correlation")+
  theme_bw()


years<-all_year_ab%>%
  group_by(ID)%>%
  mutate(year_min = min(Year), year_max = max(Year))%>%
  select(ID, year_min, year_max)%>%
  distinct()

melt_lambda_short<-merge(melt_lambda_short, years)


sum_lambdas_obs<-all_year_ab%>%
  group_by(ID)%>%
  mutate(sum_lambdas_obs = sum(Lambdas), mean_lambdas_obs = mean(Lambdas))%>%
  ungroup()


#####compare lambdas


sum_lambdas<-melt_lambda_short%>%
  group_by(ID, ldd,SD, rep_id)%>%
  filter(Year >= year_min & Year <= year_max)%>%
  summarize(sum_lambda = sum(Lambdas), mean_lambda = mean(Lambdas))%>%
  ungroup()
  


sum_lambdas$ID<-as.factor(sum_lambdas$ID)
sum_lambdas_obs$ID<-as.factor(sum_lambdas_obs$ID)


ggplot(data = sum_lambdas,aes(x=ID, y=sum_lambda, group = ID)) + 
  geom_boxplot(linetype = "solid", size = 1)+
  #geom_point(size = 2)+
  geom_point(data = sum_lambdas_obs, aes(x = ID, y = sum_lambdas_obs, group = ID), size = 5, colour = "red")+
  geom_hline(yintercept=0, linetype = "dashed")+
  theme_bw()+
  xlab("Population")+
  ylab("Sum Lambdas")+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=20))


ggplot(data = sum_lambdas,aes(x=ID, y=mean_lambda, group = ID)) + 
  geom_boxplot(linetype = "solid", size = 1)+
  #geom_point(size = 2)+
  geom_point(data = sum_lambdas_obs, aes(x = ID, y = mean_lambdas_obs, group = ID), size = 5, colour = "red")+
  geom_hline(yintercept=0, linetype = "dashed")+
  theme_bw()+
  xlab("Population")+
  ylab("Average Growth Rate")+
  theme(legend.position="none",axis.text=element_text(size=16),
        axis.title=element_text(size=20))

