#install.packages("dismo")
library(dismo)
library(demoniche)
library(maptools)
data("wrld_simpl")

wd<- getwd()
capra <-gbif("capra", "ibex", geo=FALSE)
capgeo <- subset(capra, !is.na(lon) & !is.na(lat) & (lon!="NA" & lat !="NA") | year!="NA") 
dups <- duplicated(capgeo[, c("lon", "lat")])
capg <-capgeo[!dups, ]
capg2 <- capg[capg$lon > 0 & capg$lon<25 & capg$lat > 43 & capg$lat < 47.9 & capg$year>=1985, ] 

capg2<-data.frame(capg2$lon,capg2$lat)
colnames(capg2)<-c("Longitude", "Latitude")

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

set.seed(0)
k<-4
group_pres<-kfold(sp, k)
pres_train<-sp[group_pres!=1,]
pres_test<-sp[group_pres==1,]

e<-extent(sp)+4

lf<-list.files(paste(wd, "/Bioclim/", sep=""))

first<-which(grepl("1985_bioclim", lf))
last<-which(grepl("2016_bioclim", lf))
all_years<-stack(paste(wd, "/Bioclim/",lf[first:last], sep=""))

bios<-seq(1,nlayers(all_years), by=19)

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

set.seed(10)
backg <- randomPoints(pred_nf, n=1000, ext=e, extf = 1.25)
colnames(backg) = c('lon', 'lat')
group_back <- kfold(backg, k)
backg_train <- backg[group_back != 1, ]
backg_test <- backg[group_back == 1, ]

r <- raster(pred_nf, 1)
plot(!is.na(r), col=c('white', 'light grey'), legend=FALSE, ext=e)
plot(e, add=TRUE, col='red', lwd=2)
points(backg_train, pch='-', cex=0.5, col='yellow')
points(backg_test, pch='-',  cex=0.5, col='black')
points(pres_train, pch= '+', col='green')
points(pres_test, pch='+', col='blue')


bc <- bioclim(pred_nf, sp)
#plot(bc, a=1, b=2, p=0.85)
#plot(bc)

ev <- dismo:::evaluate(pres_test, backg_test, bc, pred_nf)
plot(ev, 'ROC')

evl<- list()

group_back <- kfold(backg, k)

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

tr <- threshold(ev, 'spec_sens')
pb <- predict(pred_nf, bc, ext=e, progress='')
pb

par(mfrow=c(1,2))
plot(pb, main='Bioclim, raw values')
plot(pb > tr, main='presence/absence')
points(pres_train, pch='+')


#######Regression

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

group <- kfold(envtrain, k)
#testpres <- data.frame( extract(pred_nf, pres_test) )
#testbackg <- data.frame( extract(pred_nf, backg_test) )

library(mgcv)


gm1<-gam(pa~ s(Bio_1_1985_2016_average)+s(Bio_5_1985_2016_average)+ s(Bio_6_1985_2016_average)+ 
           s(Bio_13_1985_2016_average)+ s(Bio_15_1985_2016_average)+ s(Bio_18_1985_2016_average)+ 
           s(Bio_19_1985_2016_average), data=envtrain)

#gam_ev<-dismo:::evaluate(testpres, testbackg, gm1)
#plot(gam_ev, "ROC")

evl<- list()

for (i in 1:k){
  pres_train<-envtrain[group!=i ,]
  pres_test<-envtrain[(group==i & envtrain$pa ==1),]
  back_test<-envtrain[(group==i & envtrain$pa ==0),]
  gm1<-gam(pa~ s(Bio_1_1985_2016_average)+ s(Bio_5_1985_2016_average)+ s(Bio_6_1985_2016_average)+ s(Bio_13_1985_2016_average)+ s(Bio_15_1985_2016_average)+ s(Bio_18_1985_2016_average)+ s(Bio_19_1985_2016_average), data=pres_train)
  evl[[i]] <- dismo:::evaluate(pres_test, back_test, gm1)
  print(i)
  }

auc <- sapply( evl, function(x){slot(x, "auc")} )

gam_auc<-mean(auc)



pg <- predict(pred_nf, gm1, ext=e)
par(mfrow=c(1,2))
plot(pg, main='GAM, raw values')
tr <- threshold(gam_ev, 'spec_sens')
plot(pg > tr, main='presence/absence')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


#RandomForest

library(randomForest)

model<-pa~  Bio_1_1985_2016_average+Bio_5_1985_2016_average+ Bio_6_1985_2016_average+ Bio_13_1985_2016_average+ Bio_15_1985_2016_average+ Bio_18_1985_2016_average+ Bio_19_1985_2016_average

rf1 <- randomForest(model, data=envtrain)
rf2 <- randomForest(model, data=envtrain)

evl<- list()


for (i in 1:k){
  pres_train<-envtrain[group!=i,]
  pres_test<-envtrain[(group==i & envtrain$pa ==1),]
  back_test<-envtrain[(group==i & envtrain$pa ==0),]
  rf1 <- randomForest(model, data=pres_train)
  evl[[i]] <- dismo:::evaluate(pres_test, back_test, rf1)
}

auc <- sapply( evl, function(x){slot(x, "auc")} )

rf_auc<-mean(auc)



#erf <- dismo:::evaluate(testpres, testbackg, rf1)
#plot(erf, "ROC")


pr <- predict(pred_nf, rf1, ext=e)

par(mfrow=c(1,2))
plot(pr, main='Random Forest, regression')

tr <- threshold(erf, 'spec_sens')
plot(pr > tr, main='presence/absence')

points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


par(mfrow=c(3,1))
models <- stack(pb, pg, pr)
names(models) <- c("bioclim", "gam", "random forest")
plot(models)

#auc <- sapply(list(ev, gam_ev, erf), function(x) x@auc)
auc<-c(bioclim_auc, gam_auc, rf_auc)
w <- (auc-0.5)^2
wm <- weighted.mean( models[[c("bioclim", "gam", "random.forest")]], w)


par(mfrow=c(1,1))
plot(wm, main='weighted mean of bioclim, gam and random forest')
points(sp)

years<-1950:2016
predict_stack<-stack(paste(wd, "/Bioclim/",years,"_bioclim_variable_stack.tif", sep=""))
#predict_alps<-crop(predict_stack, e)
bio_layer_pred<-c(1,5,6,13,15,18,19)

# for (i in 1:length(years)){
#   predict_stack_year<-predict_stack[[bio_layer_pred]]
#   predict_alps_year<-crop(predict_stack_year, e)
#   names(predict_alps_year)<-c("bio5", "bio6", "bio7", "bio13", "bio19")
#   pbc<-predict(predict_alps_year, bc)
#   writeRaster(pbc, paste("D:/Fiona/Git_Method/Git_Method/Alp_SDMs/Bioclim/bioclim_",years[i],"capra_ibex.tif", sep=""), overwrite=T)
#   bio_layer_pred<-bio_layer_pred+19
#   print(i)
# }
# 


#bc, gm1,rf1


#for each year predict with each model then create ensemble prediction
plot(models)

#auc <- sapply(list(ev, gam_ev, erf), function(x) x@auc)
auc<-c(bioclim_auc, gam_auc, rf_auc)
w <- (auc-0.5)^2

#years<-1950:2016
bio_layer_pred<-c(1,5,6,13,15,18,19)


model<-pa~ Bio_1_1985_2016_average+Bio_5_1985_2016_average+ Bio_6_1985_2016_average+ Bio_13_1985_2016_average+ Bio_15_1985_2016_average+ Bio_18_1985_2016_average+ Bio_19_1985_2016_average

set.seed(0)
k<-4
group_pres<-kfold(sp, k)
pres_train<-sp[group_pres!=1,]
pres_test<-sp[group_pres==1,]

pres_train<-data.frame(sp)[,c(2,3)]   #not sure if this should be sp or pres train? sp has more points
colnames(pres_train)<-c("lon", "lat")
train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(pred_nf, train)
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )

#envtrain
bc <- bioclim(pred_nf, sp)
gm1<-gam(pa~ s(Bio_1_1985_2016_average)+ s(Bio_5_1985_2016_average)+ s(Bio_6_1985_2016_average)+ s(Bio_13_1985_2016_average)+ s(Bio_15_1985_2016_average)+s(Bio_18_1985_2016_average)+ s(Bio_19_1985_2016_average), data=envtrain)
rf1 <- randomForest(model, data=envtrain)

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


ms<-stack(paste(wd, "/Alp_SDMs/Ensembles/weighted_ensemble_sdm_", years, ".tif", sep=""))
pa_s<-stack(paste(wd, "/Alp_SDMs/Ensembles/pres_abs_weighted_ensemble_sdm_", years, ".tif", sep=""))

mean_hab<-cellStats(ms, stat="mean")
plot(years, mean_hab, type="l")

mean_pa<-cellStats(pa_s, stat="mean")
plot(years, mean_pa, type="l")

###try with presence/absence threshold - check damaris SM



values(sp_rast)
values(pb)

mat<-as.matrix(values(pb), values(sp_rast))

truth<-as.factor(values(sp_rast))
library(caret)

thresh<-seq(0.01, 1, by=0.01)

pred<-values(pb)

spec_l<-numeric()
sens_l<-numeric()
tss_l<-numeric()
for (i in 1:length(thresh)){
  pred_i<-as.factor((pred>thresh[i])*1)
  spec<-caret:::specificity(pred_i, truth)
  sens<-caret:::sensitivity(pred_i, truth)
  spec_l<-rbind(spec, spec_l)
  sens_l<-rbind(sens, sens_l)
  tss<-(spec + sens)-1
  tss_l<-rbind(tss, tss_l)
  print(i)
}


