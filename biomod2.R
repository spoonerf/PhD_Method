install.packages("biomod2",repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))
install.packages("slam",repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))
install.packages("spocc", dependencies = T,repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))
install.packages("stringi")
install.packages("ggplot2")
install.packages("zoo")
install.packages("SparseM")
library(biomod2)
library(spocc)
library(rebird)
sp<-"Tyto alba"

sloth<-occ("acridotheres tristis", from="gbif", limit=10000)
df<-occ2df(sloth)

myRespName = "Bradypusvariegatus"

df<-occ2df(sloth)

myRespXY<-data.frame(df$longitude, df$latitude)
myRespXY<-unique(na.omit(myRespXY))

myRespName = sp

env <- getData("worldclim", var="bio", res=10)

myExpl<-env[[c(3,4,7,11,12)]]
#myExpl<-env

presence.absence.raster <- function (mask.raster,species.data,raster.label="") {
  require(raster)
  
  # set the background cells in the raster to 0
  mask.raster[!is.na(mask.raster)] <- 0
  
  #set the cells that contain points to 1
  speciesRaster <- rasterize(species.data,mask.raster,field=1)
  speciesRaster <- merge(speciesRaster,mask.raster)
  
  #label the raster
  names(speciesRaster) <- raster.label
  return(speciesRaster)
}



myRaster <- env[[1]]
#myRaster<-aggregate(myRaster, fact=6)
myRaster

# create presence absence raster for foxes
myResp.ras <- presence.absence.raster(mask.raster=myRaster, species.data=myRespXY, raster.label=myRespName)
plot(myResp.ras, main=names(myResp.ras))


myRespXY <- xyFromCell(object=myResp.ras,
                       cell=which(myResp.ras[]>0))

myResp <- extract(x=myResp.ras, y=myRespXY)

# 
# myResp<-data.frame(rep(1, length(myRespXY$df.longitude)))
# 
# myResp<-data.frame(myRespXY, myResp)


myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 2,
                                     PA.nb.absences = 2000,
                                     PA.strategy = 'random')

myBiomodData


plot(myBiomodData)


# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()


myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
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
  modeling.id = paste(myRespName,"FirstModeling",sep=""))

myBiomodModelOut



# get all models evaluation
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
# print the dimnames of this object
dimnames(myBiomodModelEval)

# let's print the TSS scores of Random Forest
myBiomodModelEval["TSS","Testing.data","RF",,]


# let's print the ROC scores of all selected models
myBiomodModelEval["ROC","Testing.data",,,]

# print variable importances
#getModelsVarImport(myBiomodModelOut)

get_variables_importance(myBiomodModelOut)

myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.65),  #was 0.85, should probs be higher
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

# print summary
myBiomodEM
  
# get evaluation scores
get_evaluations(myBiomodEM)

# projection over the globe under current conditions
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

myBiomodProj

list.files("Tyto.alba/proj_current/")

# make some plots sub-selected by str.grep argument
plot(myBiomodProj, str.grep = 'SRE')

plot(myBiomodProj, str.grep = 'CTA')

plot(myBiomodProj, str.grep = 'RF')

plot(myBiomodProj, str.grep = 'MARS')

plot(myBiomodProj, str.grep = 'FDA')

# if you want to make custom plots, you can also get the projected map
myCurrentProj <- get_predictions(myBiomodProj)
plot(myCurrentProj[[9]], main=names(myCurrentProj[[9]]))


myExpl2050<- getData('CMIP5', var='bio', res=10, rcp=85, model='HG', year=50)

myExpl2050<-myExpl2050[[c(3,4,7,11,12)]]
names(myExpl2050)<-names(myExpl)

myBiomodProjFuture <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl2050,
  proj.name = 'future',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')


# make some plots, sub-selected by str.grep argument
plot(myBiomomodProj2050, str.grep = 'MARS')


my2050Proj <- get_predictions(myBiomomodProj2050)

names(my2050Proj)

now<-myCurrentProj[[9]]
fut<-my2050Proj[[9]]

diff<-(fut-now)
mean(na.omit(values(diff)))
mean(na.omit(values(fut)))
mean(na.omit(values(now)))



myBiomodEF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProj,
  output.format = '.img',
  do.stack = FALSE)

myBiomodEF

# reduce layer names for plotting convegences
plot(myBiomodEF)


proj_t2050_Myocastor_PA1_Full_AllAlgos_EMbyTSS <- stack("Myocastor/proj_t2050/proj_t2050_Myocastor_PA1_Full_AllAlgos_EMbyTSS.grd")
proj_t2050_Myocastor_PA1_Full_AllAlgos_EMbyTSS





smp_size<-round(nrow(myResp)*0.7)
train_ind <- sample(seq_len(nrow(myResp)), size = smp_size)

train <- myResp[train_ind, ]
test <- myResp[-train_ind, ]

xy<-data.frame(train[,c(1:2)])
dt<-data.frame(train[,3])

myResp<-SpatialPointsDataFrame(xy,dt)

xye<-data.frame(test[,c(1:2)])
dte<-data.frame(test[,3])

evalResp<-SpatialPointsDataFrame(xye,dte)



myBiomodData<-BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, 
                     resp.xy = myRespXY, resp.name= myRespName,
                     PA.nb.rep = 1,
                     PA.nb.absences = 200, PA.strategy = "sre")
                     
#                     eval.resp.var = dte, eval.resp.xy = xye)

plot(myBiomodData)

