install.packages('climates',,'http://www.rforge.net/')
install.packages("dismo")

library(zoo)
library(climates)
library(dismo)

######Creating Bioclim varibles for each year 1950-2016 for Europe

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

#########Creating a biomod2 SDM
library(biomod2)
library(dismo)

#getting location points for capra ibex from gbif and then cleaning this data
capra <-gbif("capra", "ibex", geo=FALSE)
capgeo <- subset(capra, !is.na(lon) & !is.na(lat)) 
dups <- duplicated(capgeo[, c("lon", "lat")])
capg <-capgeo[!dups, ]
capg2 <- capg[capg$lon > 0 & capg$lon<25 & capg$lat > 43 , ] 
capg2$presence<-rep(1)
capg2$ID<-1:nrow(capg2)
capc<-as.matrix(capg2[ , c( "lon","lat", "presence")])
xy<-as.matrix(capc[,c(1,2)])
df<-data.frame(capg2$presence)
sp<-SpatialPointsDataFrame(coords=xy, data=df)

# eval<-capg2[sample(1:nrow(capg2),sampleSize),]
# test<-capg2[!capg2$ID %in% eval$ID,]
# 
# capc<-as.matrix(test[ , c( "lon","lat", "presence")])
# xy<-as.matrix(capc[,c(1,2)])
# df<-data.frame(test$presence)
# 
# cap_eval<-as.matrix(eval[ , c( "lon","lat", "presence")])
# xy_eval<-as.matrix(cap_eval[,c(1,2)])
# df_eval<-data.frame(eval$presence)
# 
# 
# sp_test<-SpatialPointsDataFrame(coords=xy, data=df)
# sp_eval<-SpatialPointsDataFrame(coords=xy_eval, data=df_eval)

bio_1950<-stack("C:/Users/Fiona/Documents/PhD/PhD_Method/Bioclim/1950_bioclim_variable_stack.tif")
bio_1951<-stack("C:/Users/Fiona/Documents/PhD/PhD_Method/Bioclim/1951_bioclim_variable_stack.tif")

e<-extent(sp)+4

alps<-crop(bio_1950,e)
plot(alps[[1]])
alps<-stack(alps)

species<-"Capra ibex"


my_biomod_data<-BIOMOD_FormatingData(resp.var = sp,
                     expl.var =  alps,
                     #eval.resp.var = sp_eval,
                     #eval.expl.var = bio_1950,
                     resp.name = "Capra_ibex",
                     PA.nb.rep = 100,
                     PA.nb.absences = 100,
                     PA.strategy = 'random',
                     na.rm=TRUE)

#plot(my_biomod_data)

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  my_biomod_data,
  models = c('SRE','CTA','RF','MARS','FDA'),
  models.options = myBiomodOption,
  NbRunEval=3,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(species,"FirstModeling",sep=""))

myBiomodModelEval <- get_evaluations(myBiomodModelOut)

dimnames(myBiomodModelEval)

myBiomodModelEval["TSS","Testing.data","RF",,]

myBiomodModelEval["ROC","Testing.data",,,]

get_variables_importance(myBiomodModelOut)


myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.5),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )


myBiomodEM 

get_evaluations(myBiomodEM)

myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = alps,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')


myBiomodProj


list.files("Capra.ibex/proj_current/")


plot(myBiomodProj, str.grep="MARS")

myCurrentProj <- get_predictions(myBiomodProj)

ave_proj<-mean(myCurrentProj)

plot(ave_proj)





