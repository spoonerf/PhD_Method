library(dplyr)
library(raster)
library(reshape2)

####pop trend plots

spin_years = 10
######################
####Cervus elaphus####
######################

df<-read.csv("cervus_elaphus_results_all.csv")
df<-df[,-1]

dfm<-as.matrix(df)

melt_df<-melt(df, id=1:8)
melt_df$year<-as.numeric(gsub("Year_", "", melt_df$variable))

melt_short<-melt_df[melt_df$year>spin_years[length(spin_years)] ,]
melt_short$sp_lpi.ID<-as.factor(melt_short$sp_lpi.ID)

melt_short<-melt_short %>%
  filter(sp_lpi.ID == "543" | sp_lpi.ID == "6555" | sp_lpi.ID == "11176" | sp_lpi.ID == "11180")


melt_short$sp_lpi.ID <- factor(melt_short$sp_lpi.ID, levels = c("543", "6555", "11176", "11178", "11180", "11489", "11494"), 
                               labels = c("n", "o", "p", " ", "q", " "," "))

melt_short_summary<-melt_short %>%
  group_by(sp_lpi.ID, year) %>%
  mutate(mean_abund = mean(value), median_abund = median(value))


library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=interaction(ldd, SD, rep_id)))+
  geom_line(colour= alpha("grey", 0.2))+
  geom_line(data = melt_short_summary, aes(x = year, y = mean_abund ,colour = sp_lpi.ID))+
  facet_grid(~ sp_lpi.ID)+
  #facet_grid(~ sp_lpi.ID, labeller = supp.labs)+
  xlab("Year")+
  ylab("Total Predicted Abundance")+theme(legend.position = "none") 

#########Maps


library(raster)
library(viridis)

wd<-getwd()

lf<-list.files(paste(wd, "/Legion/snow_cervus_bias_new/", sep=""))
files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]
cervus<-stack(paste(wd, "/Legion/snow_cervus_bias_new/", files, sep = ""))
cervus_mean<-calc(cervus, fun = function(x){sum(diff(x))})
plot(cervus_mean, col = viridis(100))

flat<-cervus_mean
flat[!is.na(flat[])] <- 1 



demoniche_folder<-paste(wd, "/Legion/snow_cervus_bias_faster/output", sep="")
highfoldernames<-list.files(demoniche_folder, pattern ="*hyde_new_patch")


for (i in 1545:length(highfoldernames)){
  pop_out<-read.csv(paste(demoniche_folder ,highfoldernames[[i]], "Reference_matrix_pop_output.csv", sep="/"), header = TRUE)
  pop_out<-pop_out[,-1]
  coordinates(pop_out) <- ~ X + Y
  gridded(pop_out) <- TRUE
  rasterDF <- raster:::stack(pop_out)
  rasterDF<-rasterDF[[11:66]]
  cnd_ras<-calc(rasterDF, fun = function(x){sum(diff(x))})
  
  writeRaster(cnd_ras, paste0(demoniche_folder, "/stacked_raster_out/net_abundance_", highfoldernames[[i]], ".tif"), overwrite = TRUE)
  print(i)
}


all_diffs<-stack(paste0(demoniche_folder, "/stacked_raster_out/net_abundance_", highfoldernames, ".tif"))
cnd_ras<-calc(all_diffs, mean)

writeRaster(cnd_ras, "mean_cervus_abudance_change.tif")

plot(cnd_ras*flat)


#diff_ras<-scale(cervus_mean) - scale(cnd_ras)

#plot(scale(cervus_mean))

#plot(diff_ras)


par(mfrow=c(1,2)) 

plot(cervus_mean, col = viridis(100), main = "Net Change in Habitat Suitability\nRed Deer (1950 - 2005)")
plot(cnd_ras*flat, col= viridis(100), main = "Net Change in Predicted Abundance\nRed Deer (1950 - 2005)")



#Presence/Absence

files_pa<-lf[grepl("^hyde_pres_abs_sss_.*.tif$", lf)]

cervus_pa<-stack(paste(wd, "/Legion/snow_cervus_bias_new/", files_pa, sep = ""))
cervus_sum<-calc(cervus_pa,sum)
plot(cervus_sum, col = viridis(100))

files_pa<-lf[grepl("^hyde_pres_abs_sss_.*.tif$", lf)]

cervus_pa<-stack(paste(wd, "/Legion/snow_cervus_bias_new/", files_pa, sep = ""))
cervus_sum<-calc(cervus_pa,sum)
plot(cervus_sum, col = viridis(100), main = "Total Years Red Deer Predicted Presence")
plot(cnd_ras*flat, col= viridis(100), main = "Net Change in Predicted Abundance\nRed Deer (1950 - 2005)")


xy<-cbind(0,48)
fr<-extract(rasterDF, xy)
points(xy)



################
###Capra ibex###
################

df<-read.csv("capra_ibex_results_all.csv")
df<-df[,-1]

dfm<-as.matrix(df)

melt_df<-melt(df, id=1:8)
melt_df$year<-as.numeric(gsub("Year_", "", melt_df$variable))

melt_short<-melt_df[melt_df$year>spin_years[length(spin_years)] ,]
melt_short$sp_lpi.ID<-as.factor(melt_short$sp_lpi.ID)

# melt_short<-melt_short %>%
#   filter(sp_lpi.ID == "543" | sp_lpi.ID == "6555" | sp_lpi.ID == "11176" | sp_lpi.ID == "11180")
# 
# 
# melt_short$sp_lpi.ID <- factor(melt_short$sp_lpi.ID, levels = c("543", "6555", "11176", "11178", "11180", "11489", "11494"), 
#                                labels = c("n", "o", "p", " ", "q", " "," "))

melt_short_summary<-melt_short %>%
  group_by(sp_lpi.ID, year) %>%
  mutate(mean_abund = mean(value), median_abund = median(value))


library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=interaction(ldd, SD, rep_id)))+
  geom_line(colour= alpha("grey", 0.2))+
  geom_line(data = melt_short_summary, aes(x = year, y = mean_abund ,colour = sp_lpi.ID))+
  facet_grid(~ sp_lpi.ID)+
  #facet_grid(~ sp_lpi.ID, labeller = supp.labs)+
  xlab("Year")+
  ylab("Total Predicted Abundance")+theme(legend.position = "none") 

#########Maps


library(raster)
library(viridis)

wd<-getwd()

lf<-list.files(paste(wd, "/Legion/snow_capra_bias/", sep=""))
files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]
cervus<-stack(paste(wd, "/Legion/snow_capra_bias/", files, sep = ""))
cervus_mean<-calc(cervus, fun = function(x){sum(diff(x))})
plot(cervus_mean, col = viridis(100))

flat<-cervus_mean
flat[!is.na(flat[])] <- 1 



demoniche_folder<-paste(wd, "/Legion/snow_capra_bias/output", sep="")
highfoldernames<-list.files(demoniche_folder)

for (i in 1:length(highfoldernames)){
  pop_out<-read.csv(paste(demoniche_folder ,highfoldernames[[i]], "Reference_matrix_pop_output.csv", sep="/"), header = TRUE)
  pop_out<-pop_out[,-1]
  coordinates(pop_out) <- ~ X + Y
  gridded(pop_out) <- TRUE
  rasterDF <- raster:::stack(pop_out)
  rasterDF<-rasterDF[[11:66]]
  cnd_ras<-calc(rasterDF, fun = function(x){sum(diff(x))})
  
  writeRaster(cnd_ras, paste0(demoniche_folder, "/stacked_raster_out/net_abundance_", highfoldernames[[i]], ".tif"), overwrite = TRUE)
  print(i)
}


all_diffs<-stack(paste0(demoniche_folder, "/stacked_raster_out/net_abundance_", highfoldernames), ".tif")

plot(cnd_ras*flat)


diff_ras<-scale(cervus_mean) - scale(cnd_ras)

plot(scale(cervus_mean))

plot(diff_ras)


par(mfrow=c(1,2)) 

plot(cervus_mean, col = viridis(100), main = "Net Change in Habitat Suitability\nAlpine ibex (1950 - 2005)")
plot(cnd_ras*flat, col= viridis(100), main = "Net Change in Predicted Abundance\nAlpine ibex (1950 - 2005)")



#Presence/Absence

files_pa<-lf[grepl("^hyde_pres_abs_sss_.*.tif$", lf)]

cervus_pa<-stack(paste(wd, "/Legion/snow_capra_bias/", files_pa, sep = ""))
cervus_sum<-calc(cervus_pa,sum)
plot(cervus_sum, col = viridis(100))

files_pa<-lf[grepl("^hyde_pres_abs_sss_.*.tif$", lf)]

cervus_pa<-stack(paste(wd, "/Legion/snow_capra_bias/", files_pa, sep = ""))
cervus_sum<-calc(cervus_pa,sum)
plot(cervus_sum, col = viridis(100), main = "Total Years Alpine ibex Predicted Presence")
plot(cnd_ras*flat, col= viridis(100), main = "Net Change in Predicted Abundance\nAlpine ibex (1950 - 2005)")


xy<-cbind(0,48)
fr<-extract(rasterDF, xy)
points(xy)


################
###Ursus arctos###
################

df<-read.csv("capra_ibex_results_all.csv")
df<-df[,-1]

dfm<-as.matrix(df)

melt_df<-melt(df, id=1:8)
melt_df$year<-as.numeric(gsub("Year_", "", melt_df$variable))

melt_short<-melt_df[melt_df$year>spin_years[length(spin_years)] ,]
melt_short$sp_lpi.ID<-as.factor(melt_short$sp_lpi.ID)

# melt_short<-melt_short %>%
#   filter(sp_lpi.ID == "543" | sp_lpi.ID == "6555" | sp_lpi.ID == "11176" | sp_lpi.ID == "11180")
# 
# 
# melt_short$sp_lpi.ID <- factor(melt_short$sp_lpi.ID, levels = c("543", "6555", "11176", "11178", "11180", "11489", "11494"), 
#                                labels = c("n", "o", "p", " ", "q", " "," "))

melt_short_summary<-melt_short %>%
  group_by(sp_lpi.ID, year) %>%
  mutate(mean_abund = mean(value), median_abund = median(value))


library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=interaction(ldd, SD, rep_id)))+
  geom_line(colour= alpha("grey", 0.2))+
  geom_line(data = melt_short_summary, aes(x = year, y = mean_abund ,colour = sp_lpi.ID))+
  facet_grid(~ sp_lpi.ID)+
  #facet_grid(~ sp_lpi.ID, labeller = supp.labs)+
  xlab("Year")+
  ylab("Total Predicted Abundance")+theme(legend.position = "none") 

#########Maps


library(raster)
library(viridis)

wd<-getwd()

lf<-list.files(paste(wd, "/Legion/snow_bear_bias_new/", sep=""))
files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]
cervus<-stack(paste(wd, "/Legion/snow_bear_bias_new/", files, sep = ""))
cervus_mean<-calc(cervus, fun = function(x){sum(diff(x))})
plot(cervus_mean, col = viridis(100))

flat<-cervus_mean
flat[!is.na(flat[])] <- 1 



demoniche_folder<-paste(wd, "/Legion/snow_bear_bias_new/output_smooth", sep="")
highfoldernames<-list.files(demoniche_folder, pattern ="*hyde_new_patch")

for (i in 1:length(highfoldernames)){
  pop_out<-read.csv(paste(demoniche_folder ,highfoldernames[[i]], "Reference_matrix_pop_output.csv", sep="/"), header = TRUE)
  pop_out<-pop_out[,-1]
  coordinates(pop_out) <- ~ X + Y
  gridded(pop_out) <- TRUE
  rasterDF <- raster:::stack(pop_out)
  rasterDF<-rasterDF[[11:66]]
  cnd_ras<-calc(rasterDF, fun = function(x){sum(diff(x))})
  
  writeRaster(cnd_ras, paste0(demoniche_folder, "/stacked_raster_out/net_abundance_", highfoldernames[[i]], ".tif"), overwrite = TRUE)
  print(i)
}


all_diffs<-stack(paste0(demoniche_folder, "/stacked_raster_out/net_abundance_", highfoldernames, ".tif"))

plot(cnd_ras*flat)


diff_ras<-scale(cervus_mean) - scale(cnd_ras)

plot(scale(cervus_mean))

plot(diff_ras)


par(mfrow=c(1,2)) 

plot(cervus_mean, col = viridis(100), main = "Net Change in Habitat Suitability\nBrown Bear (1950 - 2005)")
plot(cnd_ras*flat, col= viridis(100), main = "Net Change in Predicted Abundance\nBrown Bear (1950 - 2005)")



#Presence/Absence

files_pa<-lf[grepl("^hyde_pres_abs_sss_.*.tif$", lf)]

cervus_pa<-stack(paste(wd, "/Legion/snow_bear_bias_new/", files_pa, sep = ""))
cervus_sum<-calc(cervus_pa,sum)
plot(cervus_sum, col = viridis(100))

files_pa<-lf[grepl("^hyde_pres_abs_sss_.*.tif$", lf)]

cervus_pa<-stack(paste(wd, "/Legion/snow_bear_bias_new/", files_pa, sep = ""))
cervus_sum<-calc(cervus_pa,sum)
plot(cervus_sum, col = viridis(100), main = "Total Years Brown Bear Predicted Presence")
plot(cnd_ras*flat, col= viridis(100), main = "Net Change in Predicted Abundance\nBrown Bear (1950 - 2005)")


xy<-cbind(0,48)
fr<-extract(rasterDF, xy)
points(xy)



