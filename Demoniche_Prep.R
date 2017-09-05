#Chapter 3 workflow

install.packages('climates',,'http://www.rforge.net/')
install.packages("dismo")

library(zoo)
library(climates)
library(dismo)

######Creating Bioclim varibles for each year 1950-2016 for Europe - this only needs to be done once

rr<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/rr_0.25deg_reg_v15.0.nc/rr_0.25deg_reg_v15.0.nc") #precipitation
tg<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/tg_0.25deg_reg_v15.0.nc/tg_0.25deg_reg_v15.0.nc") #mean temp
tn<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/tn_0.25deg_reg_v15.0.nc/tn_0.25deg_reg_v15.0.nc") #min temp
tx<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/tx_0.25deg_reg_v15.0.nc/tx_0.25deg_reg_v15.0.nc") #max temp
#dates<-as.Date((gsub("X", "",names(rr))), format="%Y.%m.%d")


rr_mon<-zApply(rr, by=as.yearmon, fun = mean)
tg_mon<-zApply(tg, by=as.yearmon, fun = mean)
tn_mon<-zApply(tn, by=as.yearmon, fun = mean)
tx_mon<-zApply(tx, by=as.yearmon, fun = mean)

jans<-seq(1,804, by=12)
years<-as.character(1950:2016)

for (i in 1:length(years)){
  jan_sel<-jans[i]
  dec_sel<-jan_sel+11
  bio_vars_all<-biovars(rr_mon[[jan_sel:dec_sel]], tn_mon[[jan_sel:dec_sel]], tx_mon[[jan_sel:dec_sel]])  
  writeRaster(bio_vars_all, paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Bioclim/", years[i],"_bioclim_variable_stack.tif", sep=""), overwrite=T)
  print(years[i])
  }


#################################Dismo SDMs

library(dismo)
library(demoniche)
library(maptools)

wd<- getwd()
capra <-gbif("capra", "ibex", geo=FALSE)
capgeo <- subset(capra, !is.na(lon) & !is.na(lat) & (lon!="NA" & lat !="NA") | year!="NA") 
dups <- duplicated(capgeo[, c("lon", "lat")])
capg <-capgeo[!dups, ]
capg2 <- capg[capg$lon > 0 & capg$lon<25 & capg$lat > 43 & capg$lat < 47.9 & capg$year>=1985, ] 
capg2<-data.frame(capg2$lon,capg2$lat)
colnames(capg2)<-c("Longitude", "Latitude")

species<-"Capra_ibex"
df2<-read.csv("LPI_pops_20160523_edited.csv")
pyr<-subset(df2, Binomial ==species & Specific_location==1)    #record 11470 had wrong longitude - in Russia!
pyrs<-pyr[,c("Longitude","Latitude")]
capg2<-rbind(capg2, pyrs)
capg2<-na.omit(capg2)
capg2$presence<-rep(1)
capg2$ID<-1:nrow(capg2)

capc<-as.matrix(capg2[ , c( "Longitude","Latitude", "presence")])
capc<-matrix(capc[complete.cases(capc)], ncol=3)
xy<-as.matrix(capc[,c(1,2)])
df<-data.frame(capc[,3])
sp<-SpatialPointsDataFrame(coords=xy, data=df)
x <- circles(sp, d=50000, lonlat=TRUE)
pol <- polygons(x)
plot(pol)
points(sp)

e<-extent(sp)+4

lf<-list.files(paste(wd, "/Bioclim/", sep=""))

first<-which(grepl("1985_bioclim", lf))
last<-which(grepl("2016_bioclim", lf))
all_years<-stack(paste(wd, "/Bioclim/",lf[first:last], sep=""))

bios<-seq(1,nlayers(all_years), by=19)

#creating a 1985-2016 average of each bioclim variable
for (i in 1:19){
  layers<-bios
  bio_layer<-mean(all_years[[layers]])
  writeRaster(bio_layer, paste(wd, "/Bioclim/Bio_",i,"_1985_2016_average.tif",sep=""), overwrite=TRUE)
  bios<-bios+1
  print(layers)
}

bio_layer_pred<-c(1,5,6,13,15,18,19)  #picking out the bioclim layers we want to use in the model
pred_nf<-stack(paste(wd, "/Bioclim/Bio_", bio_layer_pred,"_1985_2016_average.tif",sep="" ))
sp_rast<-rasterize(sp, pred_nf[[1]])
sp_rast<-sp_rast[[2]]
sp_rast[is.na(sp_rast)]<-0
sp_rast<-crop(sp_rast, e)
plot(sp_rast)

k<-4
group_pres<-kfold(sp, k)
pres_train<-sp[group_pres!=1,]
pres_test<-sp[group_pres==1,]

set.seed(10)
backg <- randomPoints(pred_nf, n=1000, ext=e, extf = 1.25)
colnames(backg) = c('lon', 'lat')
group_back <- kfold(backg, k)
backg_train <- backg[group_back != 1, ]
backg_test <- backg[group_back == 1, ]

pres_train<-data.frame(pres_train)[,c(2,3)]
colnames(pres_train)<-c("lon", "lat")
train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(pred_nf, train)
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )

#Bioclim
bc <- bioclim(pred_nf, sp)
ev <- dismo:::evaluate(pres_test, backg_test, bc, pred_nf)

group_back <- kfold(backg, k)

#group_back <- kfold(pb_train, k)
evl<- list()

for (i in 1:k){
  pres_train<-sp[group_pres!=i,]
  pres_test<-sp[group_pres==i,]
  #backg_train<-sp[group_back!=i,]
  backg_test<-backg[group_back==i,]
  bc <- bioclim(pred_nf,pres_train)
  evl[[i]] <- dismo:::evaluate(pres_test, backg_test, bc, pred_nf)
  print(i)
}

auc <- sapply( evl, function(x){slot(x, "auc")} )

bioclim_auc<-mean(auc)


library(mgcv)
gm1<-gam(pa~ s(Bio_1_1985_2016_average)+s(Bio_5_1985_2016_average)+ s(Bio_6_1985_2016_average)+ 
           s(Bio_13_1985_2016_average)+ s(Bio_15_1985_2016_average)+ s(Bio_18_1985_2016_average)+ 
           s(Bio_19_1985_2016_average), data=envtrain)

env_pres<-envtrain[envtrain$pa==1,]
env_back<-envtrain[envtrain$pa==0,]
group_pres<-  kfold(env_pres, k)
group_back<- kfold(env_back, k)

evl<- list()

for (i in 1:k){
  pres_train<-env_pres[group_pres!=i ,]
  pres_test<-env_pres[(group_pres==i) ,]
  back_test<-env_back[(group_back==i),]
  back_train<-env_back[(group_back!=i),]
  envtrain<-rbind(pres_train, back_train)
  gm1<-gam(pa~ s(Bio_1_1985_2016_average)+ s(Bio_5_1985_2016_average)+ s(Bio_6_1985_2016_average)+ s(Bio_13_1985_2016_average)+ s(Bio_15_1985_2016_average)+ s(Bio_18_1985_2016_average)+ s(Bio_19_1985_2016_average), data=envtrain)
  evl[[i]] <- dismo:::evaluate(pres_test, back_test, gm1)
  print(i)
}

auc <- sapply( evl, function(x){slot(x, "auc")} )

gam_auc<-mean(auc)

library(randomForest)

model<-pa~  Bio_1_1985_2016_average+Bio_5_1985_2016_average+ Bio_6_1985_2016_average+ Bio_13_1985_2016_average+ Bio_15_1985_2016_average+ Bio_18_1985_2016_average+ Bio_19_1985_2016_average

evl<- list()

for (i in 1:k){
  pres_train<-env_pres[group_pres!=i,]
  pres_test<-env_pres[(group_pres==i) ,]
  back_test<-env_back[(group_back==i),]
  back_train<-env_back[group_back !=i,]
  envtrain<-rbind(pres_train, back_train)
  rf1 <- randomForest(model, data=envtrain)
  evl[[i]] <- dismo:::evaluate(pres_test, back_test, rf1)
}

auc <- sapply( evl, function(x){slot(x, "auc")} )

rf_auc<-mean(auc)


#using the models to predict

pb <- predict(pred_nf, bc, ext=e, progress='') #bioclim predict
pg <- predict(pred_nf, gm1, ext=e) #gam predict
pr <- predict(pred_nf, rf1, ext=e) #random predict

models <- stack(pb, pg, pr)
names(models) <- c("bioclim", "gam", "random forest")
plot(models)

auc<-c(bioclim_auc, gam_auc, rf_auc)
w <- (auc-0.5)^2
wm <- weighted.mean( models[[c("bioclim", "gam", "random.forest")]], w)
plot(wm)

#envtrain
# bc <- bioclim(pred_nf, sp)
# gm1<-gam(pa~ s(Bio_1_1985_2016_average)+ s(Bio_5_1985_2016_average)+ s(Bio_6_1985_2016_average)+ s(Bio_13_1985_2016_average)+ s(Bio_15_1985_2016_average)+s(Bio_18_1985_2016_average)+ s(Bio_19_1985_2016_average), data=envtrain)
# rf1 <- randomForest(model, data=envtrain)

k<-4
group_pres<-kfold(sp, k)
pres_train<-sp[group_pres!=1,]
pres_test<-sp[group_pres==1,]

set.seed(10)
backg <- randomPoints(pred_nf, n=1000, ext=e, extf = 1.25)
colnames(backg) = c('lon', 'lat')
group_back <- kfold(backg, k)
backg_train <- backg[group_back != 1, ]
backg_test <- backg[group_back == 1, ]

ev <- dismo:::evaluate(pres_test, backg_test, bc, pred_nf)
tr_bc <- threshold(ev, 'spec_sens')
tr_bc_k <- threshold(ev, 'kappa')
tr_bc_p <- threshold(ev, 'prevalence')
tr_bc_e<-threshold(ev, 'equal_sens_spec')

gam_ev<-dismo:::evaluate(pres_test, backg_test, gm1, pred_nf)
tr_gam<-threshold(gam_ev, 'spec_sens')
tr_gam_k<-threshold(gam_ev, 'kappa')
tr_gam_p<-threshold(gam_ev, 'prevalence')
tr_gam_e<-threshold(gam_ev, 'equal_sens_spec')

rf_ev<-dismo:::evaluate(pres_test, backg_test, rf1, pred_nf)
tr_rf<-threshold(rf_ev, 'spec_sens')
tr_rf_k<-threshold(rf_ev, 'kappa')
tr_rf_p<-threshold(rf_ev, 'prevalence')
tr_rf_e<-threshold(rf_ev, 'equal_sens_spec')

thresh<-mean(c(tr_bc, tr_gam, tr_rf))
thresh_weighted<-weighted.mean(c(tr_bc, tr_gam, tr_rf), w)
thresh_weighted_k<-weighted.mean(c(tr_bc_k, tr_gam_k, tr_rf_k), w)
thresh_weighted_p<-weighted.mean(c(tr_bc_p, tr_gam_p, tr_rf_p), w)
thresh_weighted_e<-weighted.mean(c(tr_bc_e, tr_gam_e, tr_rf_e), w)

tholds<-c(thresh, thresh_weighted, thresh_weighted_k, thresh_weighted_p, thresh_weighted_e)
names<-c("non_weighted", "spec_sens", "kappa", "prevalence", "equal_sens_spec")

cbind(names, tholds)

years<-1950:2016
library(raster)
for (i in 1:length(years)){
  
  pred_nf<-stack(paste(wd, "/Bioclim/", years[i], "_bioclim_variable_stack.tif", sep="" ))  
  pred_nf<-pred_nf[[bio_layer_pred]]
  names(pred_nf)<-c("Bio_1_1985_2016_average","Bio_5_1985_2016_average", "Bio_6_1985_2016_average", "Bio_13_1985_2016_average" ,"Bio_15_1985_2016_average","Bio_18_1985_2016_average", "Bio_19_1985_2016_average")
  
  pb <- predict(pred_nf, bc, ext=e, progress='')
  pg <- predict(pred_nf, gm1, ext=e)
  pr <- predict(pred_nf, rf1, ext=e)
  
  models <- stack(pb, pg, pr)
  names(models) <- c("bioclim", "gam", "random forest")
  wm <- weighted.mean( models[[c("bioclim", "gam", "random.forest")]], w)
  
  pa<-wm>thresh_weighted
  pa_k<-wm>thresh_weighted_k
  pa_p<-wm>thresh_weighted_p
  pa_e<-wm>thresh_weighted_e
  
  writeRaster(wm , paste(wd, "/Alp_SDMs/Ensembles/weighted_ensemble_sdm_", years[i], ".tif", sep=""), overwrite=TRUE)
  writeRaster(pa , paste(wd, "/Alp_SDMs/Ensembles/pres_abs_weighted_ensemble_sdm_", years[i], ".tif", sep=""), overwrite=TRUE)
  writeRaster(pa_k , paste(wd, "/Alp_SDMs/Ensembles/pres_abs_kappa_weighted_ensemble_sdm_", years[i], ".tif", sep=""), overwrite=TRUE)
  writeRaster(pa_p , paste(wd, "/Alp_SDMs/Ensembles/pres_abs_prevalence_weighted_ensemble_sdm_", years[i], ".tif", sep=""), overwrite=TRUE)
  writeRaster(pa_e , paste(wd, "/Alp_SDMs/Ensembles/pres_abs_equal_sens_spec_weighted_ensemble_sdm_", years[i], ".tif", sep=""), overwrite=TRUE)
  
  print(years[i])
  print(cellStats(pa, "max"))
  print(cellStats(pa_k, "max"))
  print(cellStats(pa_p, "max"))
  print(cellStats(pa_e, "max"))
  #plot(wm, main=years[i])
}

############################# Demoniche

library(demoniche)
library(rgdal)
library(raster)
species<-"Capra_ibex"

df2<-read.csv("LPI_pops_20160523_edited.csv")
pyr<-subset(df2, Binomial ==species & Specific_location==1)    #record 11470 had wrong longitude - in Russia!
tl<-pyr[,c(65:130)]
tl[tl=="NULL"]<-NA
tlm<-as.matrix(tl)
tlmn<-as.numeric(tlm)
tl<-matrix(tlmn, nrow=10, ncol=66)
ann_sum<-colSums(tl, na.rm=T)
plot(1950:2015, ann_sum)

plot(pyr$Longitude, pyr$Latitude)

#formatting the data for use in demoniche
pyrs<-pyr[,c("ID","Longitude","Latitude")]

id<-pyrs$ID*100
lam<-rep(1,length(id))    #not sure what the value here pertains to - think it sets starting population so should use values from LPI?
pyrxy<-SpatialPoints(pyr[,c("Longitude","Latitude")])

wd<-getwd()
sdm<-raster(paste(wd,"/Alp_SDMs/Ensembles/pres_abs_weighted_ensemble_sdm_1950.tif", sep=""))
e2<-extent(sdm)

r<-raster(e2, resolution=res(sdm))

rz<-rasterize(pyrxy,r,lam )
rid<-rasterize(pyrxy,r,id)
plot(rz)
plot(rid)

rz_spdf<-xyFromCell(rz, 1:ncell(rid))

rzm<-as.vector(rz)
ridm<-as.vector(rid)

df<-data.frame(ridm,rz_spdf,rzm)
colnames(df)<-c( "PatchID","X","Y","area")

Populations<-data.frame(na.omit(df)) 


sdm_df<-data.frame(ID=1:ncell(rid))

#formatting data for demoniche
for (i in 1:length(years)){
  
  #a selection of different threshold techniques for presence absence, as well as a suitability surface
  sdm<-raster(paste(wd,"/Alp_SDMs/Ensembles/weighted_ensemble_sdm_", years[i],".tif", sep=""))   #14.6
  #sdm<-raster(paste(wd,"/Alp_SDMs/Ensembles/pres_abs_weighted_ensemble_sdm_", years[i],".tif", sep="")) #15.4%  
  #sdm<-raster(paste(wd,"/Alp_SDMs/Ensembles/pres_abs_kappa_weighted_ensemble_sdm_", years[i],".tif", sep="")) #21.1%  
  #sdm<-raster(paste(wd,"/Alp_SDMs/Ensembles/pres_abs_prevalence_weighted_ensemble_sdm_", years[i],".tif", sep="")) #14,3%  
  #sdm<-raster(paste(wd,"/Alp_SDMs/Ensembles/pres_abs_equal_sens_spec_weighted_ensemble_sdm_", years[i],".tif", sep=""))  #13.1%
  if (i ==1){
    vec<-as.data.frame(sdm, xy = TRUE)      
  } else{
    vec<-as.data.frame(sdm)  
  }
  
  sdm_df<-cbind(sdm_df, vec)
  print(i)
}

sdm_df$ID[which(!is.na(df$PatchID))]<-df$PatchID[!is.na(df$PatchID)]

niche_map_mine<-sdm_df
colnames(niche_map_mine)[1:3]<-c("gridID", "X", "Y")

col_years<-paste("Year_", 1950:2016, sep="")
colnames(niche_map_mine)[4:length(colnames(niche_map_mine))]<-col_years

niche_formulas <- as.formula(paste(paste(colnames(niche_map_mine)[-c(1:3)], collapse="+"),"X+Y",sep="~"))

print(levelplot(niche_formulas, niche_map_mine, col.regions=rev(heat.colors(100)), main = "Niche Values"))

no_yrs_mine<-1 #number of years each time period represents - could be 1 when I do it for reals


####matrix set up


load(paste(wd, "COMADRE_v.2.0.1.RData", sep="/"))

#load(paste(wd, "COMADRE_v.1.0.0.RData", sep="/"))

species<-"Capra_ibex"

tempMetadata<-subset(comadre$metadata, SpeciesAccepted==species)

keep<-as.numeric(rownames(tempMetadata))

tempMat<-comadre$mat[keep]   #MatA is matrix pop model, can be split into U, F and/or C


MatList<-list(tempMat[[1]][[1]])  #varies depending on number of matrices - need to find a way to code this better - now have five matrices available so need to sort this
AllMat<-unlist(MatList)
matrices<-matrix(AllMat, ncol=length(MatList))
colnames(matrices)<- c("Reference_matrix")

#tempMat<-(tempMat[[2]][[1]] + tempMat[[3]][[1]] + tempMat[[4]][[1]]+ tempMat[[5]][[1]])/4
#MatList<-list(tempMat[[2]][[1]], tempMat[[3]][[1]] ,tempMat[[4]][[1]], tempMat[[5]][[1]])
# AllMat<-unlist(MatList)
# matrices<-matrix(AllMat, ncol=length(MatList))
# #matrices<-matrix(tempMat, ncol=1)
# colnames(matrices)<- c("Reference_matrix", "Matrix_1", "Matrix_2", "Matrix_3")

prob_scenario<-c(0.5,0.5)    #need to check this

noise<-0.90     #need to check this

stages<-comadre$matrixClass[keep][[1]]$MatrixClassAuthor
#stages<-comadre$matrixClass[keep][[2]]$MatrixClassAuthor
#stagesf<-stages[1:3]

list_names_matrices<-colnames(matrices)

sumweight<-rep(1, length(stages)) #weight of stages  - should be equal for all mine just 
#in plants seed not included in calculating population sizes - or if you wanted to just 
#calculate the female population it would be c(1,1,1,0,0,0)
#sumweightf<-c(1,1,1)

transition_affected_niche<-"all"    #which parts of the matrix are affected by the 
#niche values

transition_affected_env <- "all"

transition_affected_demogr <- "all"

env_stochas_type<-"normal"   #can also be lognormal

matrices_var <- matrix(0.01, ncol = 1, nrow = nrow(matrices), dimnames = list(NULL, "sd")) 
#standard deviation of matrices

proportion_initial<- rep(1/length(stages), length(stages)) #proportion of population in 
#each stage - no idea what this should be and will likely have a big impact on results! 
#- just doing eqaul splits for now
#proportion_initialf<- c(1/3,1/3,1/3)

density_individuals <- 4292.32  #4292.32 to 16096.2 based on density being between 8 and 30 per 100 ha and the area of each cell being 53654 ha 

K<-10000   #carrying capacity

K_weight<-c(rep(1, length(stages)))  #the weight with which carrying capacity affects each stage was FALSE

dispersal_constants_mine <-c(0.7,0.7,0.1,2)

niche_map_mine[is.na(niche_map_mine)]<-0

demoniche_setup(modelname = "Capra_ibex",Populations = Populations, Nichemap = niche_map_mine,
                matrices = matrices,matrices_var = matrices_var, prob_scenario = prob_scenario,
                stages = stages, proportion_initial = proportion_initial,
                density_individuals = density_individuals,
                fraction_LDD = 0.1, fraction_SDD = 0.2,
                dispersal_constants = dispersal_constants_mine,
                transition_affected_niche = transition_affected_niche,
                transition_affected_demogr = transition_affected_demogr,
                transition_affected_env=transition_affected_env,
                env_stochas_type = env_stochas_type,
                no_yrs = no_yrs_mine, K=K, Kweight = K_weight, 
                sumweight =sumweight)

RPyran_disp_niche <- demoniche_model(modelname = "Capra_ibex", Niche = TRUE, 
                                     Dispersal = TRUE, repetitions = 1,
                                     foldername = "RPyran_disp_niche")

RPyran_disp_niche[,"Meanpop","Reference_matrix"]

years<-1950:2016
#plot(RPyran_disp_niche[,"Meanpop","Reference_matrix"])
bleh<-(RPyran_disp_niche[,"Meanpop","Reference_matrix"])

plot(years, bleh, type="l")

# RPyran_niche <- demoniche_model(modelname = "Capra_ibex", Niche = TRUE, 
#                                 Dispersal = FALSE, repetitions = 1,
#                                 foldername = "RPyran_minimal")
# 
# RPyran_niche[,"Meanpop","Reference_matrix"]
# 
# years<-1950:2016
# bleh<-(RPyran_niche[,"Meanpop","Reference_matrix"])
# max<-(RPyran_niche[,"Max","Reference_matrix"])
# min<-(RPyran_niche[,"Min","Reference_matrix"])
# plot(bleh, type="l")
# lines(max, col="red")
# lines(min, col="red")
# 
# plot(years,log10(bleh), type="l")




niche_disp<-cbind(years,RPyran_disp_niche[,"Meanpop","Reference_matrix"], RPyran_disp_niche[,"Max","Reference_matrix"])

#lines(niche_disp, type="l", col="red")
plot(niche_disp, type="l", col="red")
mean(diff(log10(niche_disp[,2])))
























