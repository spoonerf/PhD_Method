

#There are some NAs in rep_id - not sure why..... find out!
#need to take out the population trends which crashed out and make a note of which ones these were
#need to talk to damaris about optimising the matrix

library(reshape2)

lpi<-read.csv("LPI_pops_20160523_edited.csv")

spin_years<-1850:1949
years<-1950:2005

binomial = "Cervus_elaphus"
demoniche_folder<-"C:/Users/Fiona/net_docs/Red_Deer_new/Output"

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

demoniche_pop_out<-lapply(highfoldernames[1:100], convert_pop_out)
df <- do.call("rbind", demoniche_pop_out)
dfm<-as.matrix(df)

lambda<-function(x){
  
  l10<-10^diff(log1p(as.numeric(x[10:length(x)])))
  
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


#write.csv(melt_short, "cervus_elaphus_melt_short.csv")
#write.csv(melt_lambda_short, "cervus_elaphus_melt_lambda_short.csv")


melt_short<-read.csv("cervus_elaphus_melt_short.csv")
melt_lambda_short<-read.csv("cervus_elaphus_melt_lambda_short.csv")


melt_df<-melt_df

library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=(rep_id), colour= sp_lpi.ID))+
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

pops<-sp_lpi[,c(1,65:120)]
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
  } else {
    lint<-approx(df$Population, n = length(df$Population))$y
    ial<-data.frame(id, Year, lint)
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

both_df_ab<-rbind(mldab, all_year_ab)

both_df_ab[both_df_ab$sp_lpi.ID == "  539",]$sp_lpi.ID<-539

library(ggplot2)
library(dplyr)
library(mgcv)
library(purrr)

model_gam<-function(df) {
  gam(value~s(year), data=df)
}

gammy<-function(df, start_year = 1950, end_year = 2005){
  
  gam_out<-df %>%
    split(.$sp_lpi.ID) %>%
    map(.,model_gam)%>%
    map(fitted.values)%>%
    map_dfr(unique)
  
  gam_out<-melt(gam_out)
  gam_out<-data.frame(gam_out, start_year:end_year)
  colnames(gam_out)<-c("sp_lpi.ID", "Abundance", "Year")
  
  return(gam_out)
}

gam_r<-gammy(melt_short)


#remove time series of all zeros - identify which these are

ggplot(mldab, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_r, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()

mldab_no_z<-group_by(mldab, sp_lpi.ID, med_disp, sdd, ldd, kern,SD, dens, link, rep_id)%>%
  mutate(value_sum = sum(value)) %>% 
  filter(value_sum >= 0)
  
group_by(mldab, sp_lpi.ID, med_disp, sdd, ldd, kern,SD, dens, link, rep_id)%>%
  mutate(value_sum = sum(value)) %>% 
  filter(value_sum <= 0)%>%
  distinct(sp_lpi.ID, med_disp, sdd, ldd, kern,SD, dens, link, rep_id)


#plots without low populations

gam_no_z<-gammy(mldab_no_z)

ggplot(mldab_no_z, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_no_z, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()

#dispersal = 0.1 & sd = 0.1

mldab_no_z_sdd_0_1_sd_0_1<-filter(mldab_no_z, sdd == 0.1 & ldd ==0.1 & SD == 0.1)

gam_no_z_sdd_0_1_sd_0_1<-gammy(mldab_no_z_sdd_0_1_sd_0_1)

ggplot(mldab_no_z_sdd_0_1_sd_0_1, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_no_z_sdd_0_1_sd_0_1, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()


#dispersal = 0.1 & sd = 0.25

mldab_no_z_sdd_0_1_sd_0_25<-filter(mldab_no_z, sdd == 0.1 & ldd ==0.1 & SD == 0.25)

gam_no_z_sdd_0_1_sd_0_25<-gammy(mldab_no_z_sdd_0_1_sd_0_25)

ggplot(mldab_no_z_sdd_0_1_sd_0_25, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_no_z_sdd_0_1_sd_0_25, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()


#dispersal = 0.25 & sd = 0.1


mldab_no_z_sdd_0_25_sd_0_1<-filter(mldab_no_z, sdd == 0.25 & ldd ==0.25 & SD == 0.1)

gam_no_z_sdd_0_25_sd_0_1<-gammy(mldab_no_z_sdd_0_25_sd_0_1)

ggplot(mldab_no_z_sdd_0_25_sd_0_1, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_no_z_sdd_0_25_sd_0_1, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()

#dispersal = 0.25 & sd = 0.25


mldab_no_z_sdd_0_25_sd_0_25<-filter(mldab_no_z, sdd == 0.25 & ldd ==0.25 & SD == 0.25)

gam_no_z_sdd_0_25_sd_0_25<-gammy(mldab_no_z_sdd_0_25_sd_0_25)

ggplot(mldab_no_z_sdd_0_25_sd_0_25, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_no_z_sdd_0_25_sd_0_25, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()

#dispersal = 0.5 & SD = 0.1

mldab_no_z_sdd_0_5_sd_0_1<-filter(mldab_no_z, sdd == 0.5 & ldd ==0.5 & SD == 0.1)

gam_no_z_sdd_0_5_sd_0_1<-gammy(mldab_no_z_sdd_0_5_sd_0_1)

ggplot(mldab_no_z_sdd_0_5_sd_0_1, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_no_z_sdd_0_5_sd_0_1, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()

#dispersal = 0.5 & SD = 0.25
#this one not working cos some pops have no pops that didn't crash out

mldab_no_z_sdd_0_5_sd_0_25<-filter(mldab_no_z, sdd == 0.5 & ldd ==0.5 & SD == 0.25)

gam_no_z_sdd_0_5_sd_0_25<-gammy(mldab_no_z_sdd_0_5_sd_0_25)

ggplot(mldab_no_z_sdd_0_5_sd_0_25, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_no_z_sdd_0_5_sd_0_25, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()



######too many NAs

#melt_lambda_short[melt_lambda_short$sp_lpi.ID == "  539",]$sp_lpi.ID<-539




####plotting in lambda


pops<-sp_lpi[,c(1,65:120)]
colnames(pops)[2:ncol(pops)]<-paste("Year", 1950:2005, sep="_")
pops[pops=="NULL"]<-NA
pops$rep_id<-"Observed"
pops$med_disp<-"Observed"
pops$sdd<-"Observed"
pops$ldd<-"Observed"
pops$kern<-"Observed"
pops$SD<-"Observed"
pops$dens<-"Observed"
pops$link<-"Observed"


library(taRifx)
popsm<-as.matrix(pops)

#change this to purrr?

gam_lpi<-function(x){
  #subsetting the population data by each population
  spid = x[2:(length(x)-8)]                     #subsetting only the dates
  names(spid)<-1950:2005              #renaming the date column names as R doesn't like numbered column names
  spid<-as.numeric(spid)
  pop_datab <- (!is.na(spid) )
  points = sum(pop_datab)
  id<-x[1]
  id<-as.numeric(id)
  Date<-1950:2005
  spidt<-destring(t(spid))
  time<-length(min(which(!is.na(spidt))):max(which(!is.na(spidt))))
  missing<-time-points
  
  Year<-Date[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  Population<-spidt[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  Population[Population == 0] <- mean(Population, na.rm=TRUE)*0.01 #if a population is zero one year thhis is replaced with 1% of the average population estimate - because you can log zeros
  
  df<-data.frame(Year,Population)
  

  if (points >=6) {
    PopN = log1p(df$Population)
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
    
    ial<-data.frame(id, Year[-length(Year)], exp(lambda2))
    
    colnames(ial)<-c("ID", "Year", "R")
  }
  
  return(ial)
}

gam_lpi_r<-apply(popsm,  1, gam_lpi)
gam_r<-do.call( "rbind", gam_lpi_r)

fill<-data.frame(rep(pops$ID, each=length(1950:2005)), 1950:2005)
colnames(fill)<-c("ID", "Year")

all_year_r<-join(fill, gam_r, type="right")

all_year_r$med_disp<-"Observed"
all_year_r$sdd<-"Observed"
all_year_r$ldd<-"Observed"
all_year_r$kern<-"Observed"
all_year_r$SD<-"Observed"
all_year_r$dens<-"Observed"
all_year_r$link<-"Observed"
all_year_r$rep_id<-"Observed"

colnames(all_year_r)[1:3]<-c("sp_lpi.ID", "year", "value")

mld<-melt_lambda_short[,-10]
all_year_r$sp_lpi.ID<-as.factor(all_year_r$sp_lpi.ID)
all_year_r$med_disp<-as.factor(all_year_r$med_disp)
all_year_r$rep_id<-as.factor(all_year_r$rep_id)
all_year_r$value<-as.numeric(all_year_r$value)
all_year_r$year<-as.numeric(all_year_r$year)

both_df<-rbind(mld, all_year_r)




#should consider splitting both_df anf gam_r_lambda into several plots based on each of the variables e.g one facet plot for LDD = 0.5 - yes do this - copy as above


gam_r_lambda<-gammy(melt_lambda_short, 1951, 2005)

both_df[both_df$sp_lpi.ID == "  539",]$sp_lpi.ID<-539
gam_r_lambda$sp_lpi.ID<-as.numeric(as.character(gam_r_lambda$sp_lpi.ID))

ggplot(both_df, aes(x= year, y=value, group=interaction(rep_id, med_disp,sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.2)+
  geom_line(data=both_df[both_df$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_r_lambda, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)



mldab_no_z_sdd_0_1_sd_0_1<-filter(mldab_no_z, sdd == 0.1 & ldd ==0.1 & SD == 0.1)

gam_no_z_sdd_0_1_sd_0_1<-gammy(mldab_no_z_sdd_0_1_sd_0_1)

ggplot(mldab_no_z_sdd_0_1_sd_0_1, aes(x= year, y=value, group=interaction(rep_id, med_disp, sp_lpi.ID), colour= sp_lpi.ID))+
  geom_line(colour="grey", alpha = 0.4)+
  geom_line(data=both_df_ab[both_df_ab$med_disp=="Observed",], aes(x=year, y=value), colour="red")+
  geom_line(data =  gam_no_z_sdd_0_1_sd_0_1, aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )+
  facet_grid(.~ sp_lpi.ID)+
  theme_bw()







# 
# for (i in 1:length(unique(both_df$sp_lpi.ID))){
#   pp<-ggplot(both_df[both_df$sp_lpi.ID == unique(both_df$sp_lpi.ID)[i],], aes(x= year, y=value, group=interaction(rep_id, max_disp)))+
#     geom_line(colour="grey")+
#     geom_line(data=both_df[both_df$max_disp=="Observed" & both_df$sp_lpi.ID == both_df$sp_lpi.ID[i],], aes(x=year, y=value), colour="red")+
#     geom_line(data =  gam_r_lambda[gam_r_lambda$sp_lpi.ID ==gam_r_lambda$sp_lpi.ID[i], ], aes(x =Year, y = Abundance, group=sp_lpi.ID), colour = "blue" )
#   plot(pp)
# }
# 
# 
