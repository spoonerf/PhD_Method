
####pop trend plots

df<-read.csv("cervus_elaphus_results_all.csv")
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

melt_short<-melt_short %>%
  filter(sp_lpi.ID == "543" | sp_lpi.ID == "6555" | sp_lpi.ID == "11176" | sp_lpi.ID == "11180")


melt_short_summary<-melt_short %>%
  group_by(sp_lpi.ID, year) %>%
  mutate(mean_abund = mean(value))

supp.labs <- c("n","o", "p", "q")
names(supp.labs) <- c("543", "6555", "11176", "11180")


library(ggplot2)
ggplot(melt_short, aes(x= year, y=value, group=interaction(ldd, SD, rep_id)))+
  geom_line(colour= alpha("grey", 0.5))+
  geom_line(data = melt_short_summary, aes(x = year, y = mean_abund ,colour = sp_lpi.ID))+
  #geom_smooth()+
  facet_grid(~ sp_lpi.ID, labeller = supp.labs)+
  xlab("Year")+
  ylab("Total Predicted Abundance")

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
highfoldernames<-list.files(demoniche_folder)

pop_out<-read.csv(paste(demoniche_folder ,foldername, "Reference_matrix_pop_output.csv", sep="/"), header = TRUE)
pop_out<-pop_out[,-1]
coordinates(pop_out) <- ~ X + Y
gridded(pop_out) <- TRUE
rasterDF <- raster:::stack(pop_out)

rasterDF<-rasterDF[[11:66]]


cnd_ras<-calc(rasterDF, fun = function(x){sum(diff(x))})
plot(cnd_ras*flat)


diff_ras<-scale(cervus_mean) - scale(cnd_ras)

plot(scale(cervus_mean))

plot(diff_ras)


par(mfrow=c(1,2)) 

plot(cervus_mean, col = viridis(100), main = "Net Change in Habitat Suitability\nRed deer (1950 - 2005)")
plot(cnd_ras*flat, col= viridis(100), main = "Net Change in Predicted Abundance\nRed deer (1950 - 2005)")



#Presence/Absence

files_pa<-lf[grepl("^hyde_pres_abs_sss_.*.tif$", lf)]

cervus_pa<-stack(paste(wd, "/Legion/snow_cervus_bias_new/", files_pa, sep = ""))
cervus_sum<-calc(cervus_pa,sum)
plot(cervus_sum, col = viridis(100))

xy<-cbind(0,48)
fr<-extract(rasterDF, xy)
points(xy)


####Alpine Ibex

library(raster)
library(viridis)

wd<-getwd()

lf<-list.files(paste(wd, "/Legion/snow_capra_bias/", sep=""))
files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]
capra<-stack(paste(wd, "/Legion/snow_capra_bias/", files, sep = ""))
capra_mean<-calc(cervus, fun = function(x){sum(diff(x))})
plot(capra_mean, col = viridis(100))

flat<-capra_mean
flat[!is.na(flat[])] <- 1 



demoniche_folder<-paste(wd, "/Legion/snow_capra_bias/output", sep="")
highfoldernames<-list.files(demoniche_folder)

pop_out<-read.csv(paste(demoniche_folder ,highfoldernames[[1]], "Reference_matrix_pop_output.csv", sep="/"), header = TRUE)
pop_out<-pop_out[,-1]
coordinates(pop_out) <- ~ X + Y
gridded(pop_out) <- TRUE
rasterDF <- raster:::stack(pop_out)

rasterDF<-rasterDF[[11:66]]


cnd_ras<-calc(rasterDF, fun = function(x){sum(diff(x))})
plot(cnd_ras*flat)


diff_ras<-scale(cervus_mean) - scale(cnd_ras)

plot(scale(cervus_mean))

plot(diff_ras)


par(mfrow=c(1,2)) 

plot(cervus_mean, col = viridis(100), main = "Net Change in Habitat Suitability\nAlpine ibex (1950 - 2005)")
plot(cnd_ras*flat, col= viridis(100), main = "Net Change in Predicted Abundance\nAlpine ibex (1950 - 2005)")



#Presence/Absense

files_pa<-lf[grepl("^hyde_pres_abs_sss_.*.tif$", lf)]

cervus_pa<-stack(paste(wd, "/Legion/snow_cervus_bias_new/", files_pa, sep = ""))
cervus_sum<-calc(cervus_pa,sum)
plot(cervus_sum, col = viridis(100))

xy<-cbind(0,48)
fr<-extract(rasterDF, xy)
points(xy)





