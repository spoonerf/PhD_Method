

#There are some NAs in rep_id - not sure why..... find out!
#need to take out the population trends which crashed out and make a note of which ones these were
#need to talk to damaris about optimising the matrix

library(reshape2)
library(sp)
library(raster)
library(mgcv)
wd<-getwd()
lpi<-read.csv("LPI_pops_20160523_edited.csv")

spin_years<-1940:1949
years<-1950:2005

binomial = "Capra_ibex"
demoniche_folder<-"C:/Users/Fiona/Documents/PhD/PhD_Method/Legion/cervus_output"
demoniche_folder<-"D:/Fiona/Git_Method/Git_Method/Legion/snow_cervus_bias/output_new"
demoniche_folder<-paste(wd, "/Legion/snow_capra_bias/output", sep="")


l<-list.files(demoniche_folder)
nf<-length(list.files(paste(demoniche_folder, l[1], sep="/")))

highfoldernames<-list.files(demoniche_folder)
lowfoldernames<-rep(1:nf, each=length(l))

foldernames<-paste(highfoldernames, lowfoldernames, sep="/")

sp_lpi<-lpi[lpi$Binomial == binomial & lpi$Specific_location ==1 ,]

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
    SD<-strsplit(foldername, "[/_]")[[1]][2]
    sdd<-strsplit(foldername, "[/_]")[[1]][3]
    ldd<-strsplit(foldername, "[/_]")[[1]][4]
    dens<-strsplit(foldername, "[/_]")[[1]][5]
    link<-strsplit(foldername, "[/_]")[[1]][6]
    med_disp<-strsplit(foldername, "[/_]")[[1]][7]
    max_disp<-strsplit(foldername, "[/_]")[[1]][8]
    rep_id<-strsplit(foldername, "[/_]")[[1]][9]
    trends_df<-data.frame(sp_lpi$ID,med_disp,sdd,ldd,SD,dens,link,rep_id,trends)
  }
}
demoniche_pop_out<-lapply(highfoldernames, convert_pop_out)
df <- do.call("rbind", demoniche_pop_out)
dfm<-as.matrix(df)

melt_df<-melt(df, id=1:8)
melt_df$year<-as.numeric(gsub("Year_", "", melt_df$variable))

melt_short<-melt_df[melt_df$year>spin_years[length(spin_years)] ,]
melt_short$sp_lpi.ID<-as.factor(melt_short$sp_lpi.ID)


library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=interaction(ldd, SD), colour= sp_lpi.ID))+
  # geom_line()+
  geom_smooth(se = FALSE)+
  facet_grid(~ sp_lpi.ID)

dfl<-split(melt_short, list(melt_short$sp_lpi.ID, melt_short$ldd, melt_short$SD))
gam_smooth<-function(x){
  if(nrow(x)>0){
    ind<-seq(length(unique(x$rep_id)),(length(unique(x$rep_id)) * length(unique(x$year))* length(unique(x$sp_lpi.ID))), by=length(unique(x$rep_id)))
    id<-as.numeric(as.character(x$sp_lpi.ID[1]))
    ldd<-as.numeric(as.character(x$ldd[1]))
    sd<-as.numeric(as.character(x$SD[1]))
    x$value[x$value ==0] <-1
    mg<-mgcv:::gam(log10(value)~s(year, bs="cs", k = -1), data = x)
    fmg<-fitted.values(mg)
    lambdas<-diff(fmg[ind])
    out<-cbind(id,ldd,sd,lambdas, unique(x$year)[-1])
    gam_out<-data.frame(out)
    colnames(gam_out)<-c("ID","ldd","SD","Lambdas", "Year")
    return(gam_out)
     }
print(paste(unique(x$sp_lpi.ID), unique(x$ldd)), sep=" " )
  }
for(i in 1:length(dfl)){

  mgcv:::gam(log10(value)~s(year, bs="cs", k = -1), data = dfl[[i]])
  print(i)

  }

df_out<-lapply(dfl, gam_smooth)
melt_lambda_short<-do.call( "rbind", df_out)

melt_lambda_short$ID<-as.factor(melt_lambda_short$ID)

ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x= Year, y=Lambdas, group=interaction(ldd, SD), colour= ID), method="loess")+
  facet_grid(~ ID)


library(taRifx)
library(plyr)
library(mgcv)

pops<-sp_lpi[,c(1,65:120)]

pops_melt<-melt(pops, id = "ID")
pops_melt$variable<-rep(1950:2005, each = length(unique(pops$ID)))

colnames(pops_melt)<-c("sp_lpi.ID","year", "obs_abundance")

abun_pred_obs<-merge(pops_melt, melt_df)


plot(abun_pred_obs$obs_abundance, abun_pred_obs$value)

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
  
  return(ial)
}

gam_lpi_r<-apply(popsm,  1, gam_lpi)
gam_r<-do.call( "rbind", gam_lpi_r)

gam_r<-gam_r[gam_r$Year <=2005,]

fill<-data.frame(rep(pops$ID, each=length(1950:2005)), 1950:2005)
colnames(fill)<-c("ID", "Year")

all_year_ab<-join(fill, gam_r, type="right")

colnames(melt_lambda_short)[1]<-"ID"

all_year_ab$ID<-as.numeric(as.character(all_year_ab$ID))
melt_lambda_short$ID<-as.numeric(as.character(melt_lambda_short$ID))

melt_lambda_short<-melt_lambda_short[melt_lambda_short$ID %in% all_year_ab$ID,]




ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x = Year, y= Lambdas, group=interaction(ldd, SD)), colour = "black", se=FALSE)+
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red")+
 # scale_x_discrete(breaks = melt_lambda_short$Year[c(10,20,30,40,50)])+
  #geom_line(data =  gam_r_lambda, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ID)


smooth_vals_pred = predict(loess(Lambdas~Year,melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[1],]), melt_lambda_short$year[melt_lambda_short$ID == melt_lambda_short$ID[1]])

smooth_vals_obs = predict(loess(Lambdas~Year,all_year_ab[all_year_ab$ID == all_year_ab$ID[1],]), all_year_ab$Year[all_year_ab$ID == all_year_ab$ID[1]])

####need to match the years up 
#rmse(smooth_vals_pred, smooth_vals_obs)

both_lambdas<-merge(melt_lambda_short, all_year_ab, by = c("Year", "ID"))

colnames(both_lambdas)[5:6]<-c("Lambdas_pred", "Lambdas_obs")
###sdm trends

plot(both_lambdas$Lambdas_obs, both_lambdas$Lambdas_pred)


sdm_folder<-paste(wd, "/Legion/snow_cervus_bias/Cervus_elaphus/SDM_folder", sep="")


sdm_stack<-stack(paste(sdm_folder, "/hyde_weighted_ensemble_sdm_", years,".tif", sep=""))
patch_stack<-stack(paste(sdm_folder, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years,".tif", sep=""))

pyr<-subset(lpi, Binomial ==binomial & Specific_location==1& Region=="Europe")    #record 11470 had wrong longitude - in Russia!
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
pyrs<-pyr[,c("ID","Longitude","Latitude")]
xy_lpi<-data.frame(pyrs$Longitude, pyrs$Latitude)

sdm_lpi<-raster:::extract(sdm_stack,xy_lpi, buffer = 50000, fun = mean, na.rm = TRUE)

colnames(sdm_lpi)<-1950:2005
sdm_lpi_melt<-melt(sdm_lpi)
colnames(sdm_lpi_melt)<-c("ID", "Year", "HSI")
sdm_lpi_melt$ID<-rep(pyrs$ID, length(1950:2005))


ggplot(sdm_lpi_melt, aes(x = Year, y= HSI, group = ID))+
  geom_smooth()+
  facet_grid(.~ID)


sdl<-split(sdm_lpi_melt, list(sdm_lpi_melt$ID))
sdm_smooth<-function(x){
  
  #ind<-seq(length(unique(x$rep_id)),(length(unique(x$rep_id)) * length(unique(x$year))* length(unique(x$sp_lpi.ID))), by=length(unique(x$rep_id)))
  id<-as.numeric(as.character(x$ID[1]))
  mg<-gam(log10(HSI)~s(Year, bs="cs", k = -1), data = x)
  fmg<-fitted.values(mg)
  lambdas<-diff(fmg)
  out<-cbind(id,lambdas, unique(x$Year)[-1])
  gam_out<-data.frame(out)
  colnames(gam_out)<-c("ID","HSI_Lambdas", "Year")
  return(gam_out)
}

df_out<-lapply(sdl, sdm_smooth)
sdm_lambdas_melt<-do.call( "rbind", df_out)


#with loess
ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x = Year, y= Lambdas, group=interaction(ldd, SD)), colour = "black",SE = FALSE, method = "loess")+   #demoniche
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red", fill="lightcoral",method = "loess")+    #real
  geom_smooth(data = sdm_lambdas_melt, aes(x = Year, y= HSI_Lambdas, group=ID), colour = "blue", fill="light blue",method = "loess")+   #sdm trend
  facet_grid(.~ID,labeller=label_both)

#with gam
ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x = Year, y= Lambdas, group=interaction(ldd, SD)), colour = "black", method = "gam", formula = y ~ s(x, bs = "cs", k = -1))+   #demoniche
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red", fill="lightcoral",method = "gam", formula = y ~ s(x, bs = "cs", k =-1))+    #real
  geom_smooth(data = sdm_lambdas_melt, aes(x = Year, y= HSI_Lambdas, group=ID), colour = "blue", fill="light blue",method = "gam", formula = y ~ s(x, bs = "cs", k = -1))+   #sdm trend
  facet_grid(.~ID,labeller=label_both)


#need start year and end year


library(Metrics)


vg<-expand.grid(unique(melt_lambda_short$ID), unique(melt_lambda_short$ldd))
colnames(vg)<-c("ID", "ldd")
rmse_get_sdm<-function(x){
  
  
  cnd_x<-melt_lambda_short[melt_lambda_short$ID == vg[x,1]& melt_lambda_short$ldd == vg[x,2],]
  obs_x<-all_year_ab[all_year_ab$ID == vg[x,1],]
  sdm_x<-sdm_lambdas_melt[sdm_lambdas_melt$ID == vg[x,1],]
  
  if (nrow(obs_x) >= 5){
    
    #need to get separate line per iteration of parameters for cnd
    cnd_gam = mgcv:::gam(Lambdas~s(Year, bs="cs", k = -1),data = cnd_x )
    obs_gam = gam(Lambdas~s(Year, bs="cs", k = -1),data = obs_x)
    sdm_gam = gam(HSI_Lambdas~s(Year, bs="cs", k =-1),data =sdm_x)
    
    smooth_vals_cnd = predict(cnd_gam,newdata = cnd_x)
    smooth_vals_obs = predict(obs_gam,newdata = obs_x)
    smooth_vals_sdm = predict(sdm_gam,newdata = sdm_x)    
    
    start_cnd<-which(unique(cnd_x$Year) == min(obs_x$Year))
    end_cnd<-which(unique(cnd_x$Year) == max(obs_x$Year))
    
    start_obs<-which(unique(obs_x$Year) == min(obs_x$Year))
    end_obs<-which(unique(obs_x$Year) == max(obs_x$Year))
    
    start_sdm<-which(unique(sdm_x$Year) == min(obs_x$Year))
    end_sdm<-which(unique(sdm_x$Year) == max(obs_x$Year))
    
    cnd_rmse<-rmse(smooth_vals_cnd[start_cnd:end_cnd], smooth_vals_obs[start_obs:end_obs])
    sdm_rmse<-rmse(smooth_vals_sdm[start_sdm:end_sdm], smooth_vals_obs[start_obs:end_obs])
    
    rmse_out<-data.frame(cnd_rmse,sdm_rmse)
    colnames(rmse_out)<-c("cnd", "sdm")
    
    return(rmse_out)
  }
}
rmse_scores<-sapply(1:nrow(vg), rmse_get_sdm)
rmse_out<-do.call( "rbind", rmse_scores)

rmse_out$diff<-rmse_out$cnd - rmse_out$sdm




ccf_get<-function(x){
  
  cnd_gam = gam(value~s(year, bs="cs"),data = melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[x],])
  obs_gam = gam(Lambdas~s(Year, bs="cs"),data = all_year_ab[all_year_ab$ID == all_year_ab$ID[x],])
  sdm_gam = gam(HSI~s(Year, bs="cs"),data = sdm_lambdas_melt[sdm_lambdas_melt$ID == sdm_lambdas_melt$ID[x],])
  
  smooth_vals_cnd = predict(cnd_gam,newdata =melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[x] & melt_lambda_short$rep_id==1,] )
  smooth_vals_obs = predict(obs_gam,newdata =all_year_ab[all_year_ab$ID == all_year_ab$ID[x],] )
  smooth_vals_sdm = predict(sdm_gam,newdata =sdm_lambdas_melt[sdm_lambdas_melt$ID == sdm_lambdas_melt$ID[x],] )    
  # #need to check this is okay on another dataset
  # smooth_vals_cnd = predict(cnd_gam,newdata =melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[x] & melt_lambda_short$rep_id==1,] )    
  # smooth_vals_obs = predict(loess(Lambdas~Year,all_year_ab[all_year_ab$ID == all_year_ab$ID[x],]), all_year_ab$Year[all_year_ab$ID == all_year_ab$ID[x]])
  # smooth_vals_sdm = predict(loess(HSI~Year,sdm_lambdas_melt[sdm_lambdas_melt$ID == sdm_lambdas_melt$ID[x],]), sdm_lambdas_melt$Year[sdm_lambdas_melt$ID == sdm_lambdas_melt$ID[x]])
  # 
  obs_x<-all_year_ab[all_year_ab$ID == all_year_ab$ID[x],]
  cnd_x<-melt_lambda_short[melt_lambda_short$ID==melt_lambda_short$ID[x],]
  sdm_x<-sdm_lambdas_melt[sdm_lambdas_melt$ID==sdm_lambdas_melt$ID[x],]
  
  start_obs<-which(unique(all_year_ab$Year) == min(obs_x$Year))
  end_obs<-which(unique(all_year_ab$Year) == max(obs_x$Year))
  
  start_cnd<-which(unique(melt_lambda_short$year) == min(obs_x$Year))
  end_cnd<-which(unique(melt_lambda_short$year) == max(obs_x$Year))
  
  start_sdm<-which(unique(sdm_lambdas_melt$Year) == min(obs_x$Year))
  end_sdm<-which(unique(sdm_lambdas_melt$Year) == max(obs_x$Year))
  
  cnd_ccf<-ccf(smooth_vals_cnd[start_cnd:end_cnd], smooth_vals_obs[start_obs:end_obs], type="correlation")
  sdm_ccf<-ccf(smooth_vals_sdm[start_sdm:end_sdm], smooth_vals_obs[start_obs:end_obs], type="correlation")
  
  lag_0_cnd<-cnd_ccf$acf[which(as.numeric(cnd_ccf$lag) == 0)]
  lag_0_sdm<-sdm_ccf$acf[which(as.numeric(sdm_ccf$lag) == 0)]
  
  ccf_out<-data.frame(lag_0_cnd, lag_0_sdm)
  colnames(ccf_out)<-c("cnd", "sdm")
  
  return(ccf_out)
}



ccf_scores<-sapply(1:length(unique(melt_lambda_short$ID)), ccf_get)

