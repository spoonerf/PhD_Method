

#There are some NAs in rep_id - not sure why..... find out!
#need to take out the population trends which crashed out and make a note of which ones these were
#need to talk to damaris about optimising the matrix

library(reshape2)

lpi<-read.csv("LPI_pops_20160523_edited.csv")

spin_years<-1940:1949
years<-1950:2005

binomial = "Cervus_elaphus"
demoniche_folder<-"C:/Users/Fiona/Documents/PhD/PhD_Method/Legion/cervus_output"
demoniche_folder<-"D:/Fiona/Git_Method/Git_Method/Legion/snow_cervus_bias/output_new"


l<-list.files(demoniche_folder)
nf<-length(list.files(paste(demoniche_folder, l[1], sep="/")))

highfoldernames<-list.files(demoniche_folder)
lowfoldernames<-rep(1:nf, each=length(l))

foldernames<-paste(highfoldernames, lowfoldernames, sep="/")

sp_lpi<-lpi[lpi$Binomial == binomial & lpi$Specific_location ==1 & lpi$Region == "Europe",]

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
  rasterDF <- stack(pop_out)
  trends<-raster:::extract(rasterDF,xy)
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
dfm<-as.matrix(df)

lambda<-function(x){
  
  l10<-diff(log1p(as.numeric(x[10:length(x)])))
  #l10<-10^diff(log1p(as.numeric(x[20:length(x)])))
  
}

dft<-t(apply(dfm,1,lambda))

df_lambda<-data.frame(dfm[,1:8],dft)

colnames(df_lambda)[9:ncol(df_lambda)]<-colnames(dfm)[9:ncol(df_lambda)]

melt_df<-melt(df, id=1:8)
melt_df$year<-as.numeric(gsub("Year_", "", melt_df$variable))

melt_lambda<-melt(df_lambda, id=1:8)
melt_lambda$year<-as.numeric(gsub("Year_", "", melt_lambda$variable))

melt_short<-melt_df[melt_df$year>spin_years[length(spin_years)] ,]
melt_short$sp_lpi.ID<-as.factor(melt_short$sp_lpi.ID)

melt_lambda_short<-melt_lambda[melt_lambda$year>years[1] ,]
melt_lambda_short$sp_lpi.ID<-as.factor(melt_lambda_short$sp_lpi.ID)


library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=rep_id, colour= sp_lpi.ID))+
  geom_line()+
  facet_grid(med_disp~ sp_lpi.ID)

ggplot(melt_lambda_short, aes(x= year, y=value, group=interaction(rep_id, med_disp), colour= sp_lpi.ID))+
  geom_line()+
  geom_smooth()+
  facet_grid(med_disp~ sp_lpi.ID)


ggplot(melt_lambda_short, aes(x= year, y=value, group=interaction(rep_id, med_disp), colour= sp_lpi.ID))+
  geom_line()+
  facet_grid(sdd ~ sp_lpi.ID)



library(taRifx)
library(plyr)
library(mgcv)

pops<-sp_lpi[,c(1,65:120)]

pop_counts<-sp_lpi[,c(65:120)]
pop_counts <- (pop_counts !="NULL")
points_per_pop1950_2005 = rowSums(pop_counts)
if(sum(points_per_pop1950_2005<5)>=1){
  #pyr<-pyr[-which(points_per_pop1950_2005<5),]
  pops<-pops[-which(points_per_pop1950_2005<5),]
}

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
  geom_smooth(data = melt_lambda_short, aes(x = year, y= value, group=ID), colour = "black")+
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red")+
  #geom_line(data =  gam_r_lambda, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ID)


smooth_vals_pred = predict(loess(value~year,melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[1],]), melt_lambda_short$year[melt_lambda_short$ID == melt_lambda_short$ID[1]])

smooth_vals_obs = predict(loess(Lambdas~Year,all_year_ab[all_year_ab$ID == all_year_ab$ID[1],]), all_year_ab$Year[all_year_ab$ID == all_year_ab$ID[1]])

####need to match the years up 
#rmse(smooth_vals_pred, smooth_vals_obs)



###sdm trends

sdm_stack<-stack(paste(sdm_folder, "/hyde_weighted_ensemble_sdm_", years,".tif", sep=""))
patch_stack<-stack(paste(sdm_folder, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years,".tif", sep=""))

pyr<-subset(lpi, Binomial ==binomial & Specific_location==1& Region=="Europe")    #record 11470 had wrong longitude - in Russia!
pops<-pyr[,c(65:120)]
pop_counts <- (pops !="NULL")
points_per_pop1950_2005 = rowSums(pop_counts)
if(sum(points_per_pop1950_2005<5)>=1){
  pyr<-pyr[-which(points_per_pop1950_2005<5),]
  pops<-pops[-which(points_per_pop1950_2005<5),]
}
pops<-data.frame(pyr$ID, pops)
colnames(pops)<-c("ID", 1950:2005)


#formatting the lpi data for use in demoniche
pyrs<-pyr[,c("ID","Longitude","Latitude")]
xy_lpi<-data.frame(pyrs$Longitude, pyrs$Latitude)

sdm_lpi<-extract(sdm_stack, xy_lpi)
colnames(sdm_lpi)<-1950:2005
sdm_lpi_melt<-melt(sdm_lpi)
colnames(sdm_lpi_melt)<-c("ID", "Year", "HSI")
sdm_lpi_melt$ID<-rep(pyrs$ID, length(1950:2005))

lambda<-function(x){
  
  lambdas<-diff(log10(as.numeric(x)))
  
}

sdm_lambdas<-t(apply(sdm_lpi,1,lambda))

colnames(sdm_lambdas)<-1950:2004
sdm_lambdas_melt<-melt(sdm_lambdas)
colnames(sdm_lambdas_melt)<-c("ID", "Year", "HSI")
sdm_lambdas_melt$ID<-rep(pyrs$ID, length(1950:2004))


#with loess
ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x = year, y= value, group=ID), colour = "black", fill="grey", method = "loess")+   #demoniche
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red", fill="lightcoral",method = "loess")+    #real
  geom_smooth(data = sdm_lambdas_melt, aes(x = Year, y= HSI, group=ID), colour = "blue", fill="light blue",method = "loess")+   #sdm trend
  facet_grid(.~ID,labeller=label_both)

#with gam
ggplot()+
  geom_smooth(data = melt_lambda_short, aes(x = year, y= value, group=ID), colour = "black", fill="grey", method = "gam", formula = y ~ s(x, bs = "cs"))+   #demoniche
  geom_smooth(data = all_year_ab, aes(x = Year, y= Lambdas, group=ID), colour="red", fill="lightcoral",method = "gam", formula = y ~ s(x, bs = "cs"))+    #real
  geom_smooth(data = sdm_lambdas_melt, aes(x = Year, y= HSI, group=ID), colour = "blue", fill="light blue",method = "gam", formula = y ~ s(x, bs = "cs"))+   #sdm trend
  facet_grid(.~ID,labeller=label_both)


#need start year and end year


rmse_get<-function(x){

  cnd_gam = gam(value~s(year, bs="cs"),data = melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[x],])
  obs_gam = gam(Lambdas~s(Year, bs="cs"),data = all_year_ab[all_year_ab$ID == all_year_ab$ID[x],])
  sdm_gam = gam(HSI~s(Year, bs="cs"),data = sdm_lambdas_melt[sdm_lambdas_melt$ID == sdm_lambdas_melt$ID[x],])
  
  smooth_vals_cnd = predict(cnd_gam,newdata =melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[x] & melt_lambda_short$rep_id==1,] )
  smooth_vals_obs = predict(obs_gam,newdata =all_year_ab[all_year_ab$ID == all_year_ab$ID[x],] )
  smooth_vals_sdm = predict(sdm_gam,newdata =sdm_lambdas_melt[sdm_lambdas_melt$ID == sdm_lambdas_melt$ID[x],] )    
  # smooth_vals_cnd = predict(loess(value~year,melt_lambda_short[melt_lambda_short$ID == melt_lambda_short$ID[x],]), unique(melt_lambda_short$year[melt_lambda_short$ID == melt_lambda_short$ID[x]]))
  # smooth_vals_obs = predict(loess(Lambdas~Year,all_year_ab[all_year_ab$ID == all_year_ab$ID[x],]), all_year_ab$Year[all_year_ab$ID == all_year_ab$ID[x]])
  # smooth_vals_sdm = predict(loess(HSI~Year,sdm_lambdas_melt[sdm_lambdas_melt$ID == sdm_lambdas_melt$ID[x],]), sdm_lambdas_melt$Year[sdm_lambdas_melt$ID == sdm_lambdas_melt$ID[x]])

  obs_x<-all_year_ab[all_year_ab$ID == all_year_ab$ID[x],]
  cnd_x<-melt_lambda_short[melt_lambda_short$ID==melt_lambda_short$ID[x],]
  sdm_x<-sdm_lambdas_melt[sdm_lambdas_melt$ID==sdm_lambdas_melt$ID[x],]

  start_obs<-which(unique(all_year_ab$Year) == min(obs_x$Year))
  end_obs<-which(unique(all_year_ab$Year) == max(obs_x$Year))

  start_cnd<-which(unique(melt_lambda_short$year) == min(obs_x$Year))
  end_cnd<-which(unique(melt_lambda_short$year) == max(obs_x$Year))

  start_sdm<-which(unique(sdm_lambdas_melt$Year) == min(obs_x$Year))
  end_sdm<-which(unique(sdm_lambdas_melt$Year) == max(obs_x$Year))

  cnd_rmse<-rmse(smooth_vals_cnd[start_cnd:end_cnd], smooth_vals_obs[start_obs:end_obs])
  sdm_rmse<-rmse(smooth_vals_sdm[start_sdm:end_sdm], smooth_vals_obs[start_obs:end_obs])

  rmse_out<-data.frame(cnd_rmse, sdm_rmse)
  colnames(rmse_out)<-c("cnd", "sdm")

return(rmse_out)
}

rmse_scores<-t(sapply(1:length(unique(melt_lambda_short$ID)), rmse_get))

boxplot(rmse_scores, use.cols = TRUE)


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

