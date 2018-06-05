library(reshape2)

lpi<-read.csv("LPI_pops_20160523_edited.csv")

spin_years<-1850:1949
years<-1950:2005

binomial = "Capra_ibex"
demoniche_folder<-"C:/Users/Fiona/Desktop/Grace_Output/Ibex/Output"

l<-list.files(demoniche_folder)
nf<-length(list.files(paste(demoniche_folder, l[1], sep="/")))

highfoldernames<-list.files(demoniche_folder)
lowfoldernames<-rep(1:nf, each=length(l))

foldernames<-paste(highfoldernames, lowfoldernames, sep="/")

sp_lpi<-lpi[lpi$Binomial == binomial & lpi$Specific_location ==1,]

xy<-data.frame(sp_lpi$Longitude, sp_lpi$Latitude)
coordinates(xy)<-c("sp_lpi.Longitude", "sp_lpi.Latitude")
proj4string(xy) <- CRS("+init=epsg:4326") # WGS 84
#CRS.new <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
#xy <- spTransform(xy, CRS.new)

convert_pop_out<-function(foldername){
  
  pop_out<-read.csv(paste(demoniche_folder ,foldername, "Reference_matrix_pop_output.csv", sep="/"), header = TRUE)
  pop_out<-pop_out[,-1]
  coordinates(pop_out) <- ~ X + Y
  gridded(pop_out) <- TRUE
  rasterDF <- stack(pop_out)
  trends<-raster:::extract(rasterDF,xy)
  med_disp<-strsplit(foldername, "[/_]")[[1]][6]
  sdd<-strsplit(foldername, "[/_]")[[1]][7]
  ldd<-strsplit(foldername, "[/_]")[[1]][8]
  kern<-strsplit(foldername, "[/_]")[[1]][9]
  SD<-strsplit(foldername, "[/_]")[[1]][10]
  dens<-strsplit(foldername, "[/_]")[[1]][11]
  link<-strsplit(foldername, "[/_]")[[1]][13]
  rep_id<-strsplit(foldername, "[/_]")[[1]][14]
  trends_df<-data.frame(sp_lpi$ID,med_disp,sdd,ldd,kern,SD,dens,link,rep_id,trends)
}

demoniche_pop_out<-lapply(foldernames, convert_pop_out)
df <- do.call("rbind", demoniche_pop_out)
dfm<-as.matrix(df)

lambda<-function(x){
  
  l10<-10^diff(log10(as.numeric(x[10:length(x)])))
  
}

dft<-t(apply(dfm,1,lambda))

df_lambda<-data.frame(dfm[,1:9],dft)

colnames(df_lambda)[10:ncol(df_lambda)]<-colnames(dfm)[11:ncol(dfm)]

melt_df<-melt(df, id=1:9)
melt_df$year<-as.numeric(gsub("Year_", "", melt_df$variable))

melt_lambda<-melt(df_lambda, id=1:9)
melt_lambda$year<-as.numeric(gsub("Year_", "", melt_lambda$variable))

melt_short<-melt_df[melt_df$year>spin_years[length(spin_years)] ,]
melt_short$sp_lpi.ID<-as.factor(melt_short$sp_lpi.ID)

melt_lambda_short<-melt_lambda[melt_lambda$year>years[1] ,]
melt_lambda_short$sp_lpi.ID<-as.factor(melt_lambda_short$sp_lpi.ID)


library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=interaction(rep_id, med_disp), colour= sp_lpi.ID))+
  geom_line()+
  facet_grid(med_disp~ sp_lpi.ID)

ggplot(melt_lambda_short, aes(x= year, y=value, group=interaction(rep_id, med_disp), colour= sp_lpi.ID))+
  geom_line()+
  facet_grid(med_disp~ sp_lpi.ID)


ggplot(melt_short, aes(x= year, y=value, group=interaction(rep_id, med_disp), colour= sp_lpi.ID))+
  geom_line()+
  facet_grid(sdd ~ sp_lpi.ID)



library(taRifx)
library(plyr)
library(mgcv)

pops<-sp_lpi[,c(1,65:130)]
colnames(pops)[2:ncol(pops)]<-paste("Year", 1950:2015, sep="_")
pops[pops=="NULL"]<-NA
pops$rep_id<-"Observed"
pops$md_id<-"Observed"

popsm<-as.matrix(pops)

gam_lpi<-function(x){
  #subsetting the population data by each population
  spid = x[2:(length(x)-2)]                     #subsetting only the dates
  names(spid)<-1950:2015              #renaming the date column names as R doesn't like numbered column names
  spid<-as.numeric(spid)
  pop_datab <- (!is.na(spid) )
  points = sum(pop_datab)
  id<-x[1]
  Date<-1950:2015
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
  
  
  if (points >=6) {
    PopN = df$Population
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
    
    
    lambda2<-pv2$fit
    
    ial<-data.frame(id, Year,lambda2)
    
    colnames(ial)<-c("ID", "Year", "Abundance")
  }
  
  return(ial)
}

gam_lpi_r<-apply(popsm,  1, gam_lpi)
gam_r<-do.call( "rbind", gam_lpi_r)

gam_r<-gam_r[gam_r$Year <=2005,]

fill<-data.frame(rep(pops$ID, each=length(1950:2005)), 1950:2005)
colnames(fill)<-c("ID", "Year")

all_year_ab<-join(fill, gam_r, type="right")

all_year_ab$med_disp<-"Observed"
all_year_ab$sdd<-"Observed"
all_year_ab$ldd<-"Observed"
all_year_ab$kern<-"Observed"
all_year_ab$SD<-"Observed"
all_year_ab$dens<-"Observed"
all_year_ab$link<-"Observed"
all_year_ab$rep_id<-"Observed"

colnames(all_year_ab)[1:3]<-c("sp_lpi.ID", "year", "value")



mldab<-melt_short[,-10]
all_year_ab$sp_lpi.ID<-as.factor(all_year_ab$sp_lpi.ID)
all_year_ab$med_disp<-as.factor(all_year_ab$med_disp)
all_year_ab$rep_id<-as.factor(all_year_ab$rep_id)
all_year_ab$value<-as.numeric(all_year_ab$value)
all_year_ab$year<-as.numeric(all_year_ab$year)


both_df_ab<-rbind(mldab, all_year_ab)

both_df_ab[both_df_ab$sp_lpi.ID == "  539",]$sp_lpi.ID<-539


library(mgcv)

gam_demon<-function(id){
  #for (i in 1:length(unique(melt_short$sp_lpi.ID))){
  
  #id<-unique(melt_short$sp_lpi.ID)[i]
  
  df<-melt_short[melt_short$sp_lpi.ID ==id,]
  PopN = df$value
  Year = df$year
  max_disp = df$max_disp
  rep_id = df$rep_id
  SmoothParm = round(length(na.omit(PopN))/2)
  
  mg2<-mgcv:::gam(PopN ~ s(Year, k=15), fx=TRUE)
  pv2<-unique(fitted.values(mg2))
  
  ial<-data.frame(id, unique(Year),pv2)
  
  colnames(ial)<-c("sp_lpi.ID", "Year", "Abundance")
  #print(id)
  return(ial)
}


gam_demon_r<-lapply(unique(melt_short$sp_lpi.ID), gam_demon)
gam_r<-do.call( "rbind", gam_demon_r)


library(mgcv)

gam_demon<-function(id){
  df<-melt_lambda_short[melt_lambda_short$sp_lpi.ID ==id,]
  PopN = df$value
  Year = df$year
  max_disp = df$max_disp
  rep_id = df$rep_id
  SmoothParm = round(length(na.omit(PopN))/2)
  
  mg2<-mgcv:::gam(PopN ~ s(Year, k=10), fx=TRUE)
  pv2<-unique(fitted.values(mg2))
  
  ial<-data.frame(id, unique(Year),pv2)
  
  colnames(ial)<-c("sp_lpi.ID", "Year", "Abundance")
  return(ial)
}

gam_demon_r_lambda<-lapply(unique(melt_lambda_short$sp_lpi.ID), gam_demon)
gam_r_lambda<-do.call( "rbind", gam_demon_r_lambda)



#There are some NAs in rep_id - not sure why..... find out!

library(ggplot2)


ggplot(mldab, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey")+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_r, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)




