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
bio_2016<-stack("C:/Users/Fiona/Documents/PhD/PhD_Method/Bioclim/2016_bioclim_variable_stack.tif")

e<-extent(sp)+4

alps<-crop(bio_1950,e)
alps_now<-crop(bio_2016,e)
plot(alps[[1]])
plot(alps_now[[1]])
alps<-stack(alps)
alps_now<-stack(alps_now)

species<-"Capra ibex"


my_biomod_data<-BIOMOD_FormatingData(resp.var = sp,
                     expl.var =  alps_now,  #alps for 1950, alps now for 2016
                     #eval.resp.var = sp_eval,
                     #eval.expl.var = bio_1950,
                     resp.name = "Capra_ibex",
                     PA.nb.rep = 10,    #switch these three off to run without pseudo absences
                     PA.nb.absences = 10,
                     PA.strategy = 'random',
                     na.rm=TRUE)

#plot(my_biomod_data)

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  my_biomod_data,
  models = c("GLM", "GAM", "RF"),
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
  eval.metric.quality.threshold = c(0.7),
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

#pres.only.eval <- BIOMOD_presenceonly(myBiomodModelOut, myBiomodEM)

myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = alps_now,    #alps for 1950 and alps now for 2016
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')


myBiomodProj


list.files("Capra.ibex/proj_current/")


plot(myBiomodProj, str.grep="RF")

myCurrentProj <- get_predictions(myBiomodProj)

ave_proj<-mean(myCurrentProj)

old<-raster("Average_SDM_1950.tif")

plot(ave_proj)

plot(myCurrentProj[[11]])

test<-myCurrentProj[[11]]/10

writeRaster(ave_proj, "Average_SDM_2016.tif") #glm, gam, rf


writeRaster(test, "Example_biomod_SDM.tif", overwrite=T)#3
writeRaster(test, "Example_biomod_SDM_loss.tif", overwrite=T)#11

test[test<50]<-0
test[test>50]<-1
plot(test)


writeRaster(test, "Example_presence_absence_biomod_SDM.tif", overwrite=T)
writeRaster(test, "Example_presence_absence_biomod_loss_SDM.tif", overwrite=T)

############## Demoniche
install.packages("popbio")
install.packages("demoniche", repos="http://R-Forge.R-project.org")
library(demoniche)

sdm<-raster("Average_SDM_1950.tif")
sdm_loss<-raster("Average_SDM_2016.tif")
sdm_table<-as.data.frame(sdm, xy = TRUE)
sdm_loss_table<-as.data.frame(sdm_loss, xy = TRUE)

dir<-getwd()
load(paste(dir, "COMADRE_v.1.0.0.RData", sep="/"))

species<-"Capra_ibex"

tempMetadata<-subset(comadre$metadata, SpeciesAccepted==species)

keep<-as.numeric(rownames(tempMetadata))

tempMat<-comadre$mat[keep]   #MatA is matrix pop model, can be split into U, F and/or C

MatList<-list(tempMat[[1]][[1]])  #varies depending on number of matrices - need to find a way to code this better
AllMat<-unlist(MatList)
matrices<-matrix(AllMat, ncol=length(MatList))
colnames(matrices)<- c("Reference_matrix")

#setting up

pyr<-subset(df2, Binomial ==species)    #record 11470 had wrong longitude - in Russia!

plot(pyr$Longitude, pyr$Latitude)

pyrs<-pyr[,c("ID","Longitude","Latitude","lambda_mean")]

id<-pyrs$ID
#lam<-as.numeric(pyrs$lambda_mean)
lam<-rep(200,length(id))    #not sure what the value here pertains to - think it sets starting population so should use values from LPI?


library(rgdal)

pyrxy<-SpatialPoints(pyr[,c("Longitude","Latitude")])

library(raster)

# e<-extent(pyrxy)
# xmn<-floor(e[1])
# xmx<-ceiling(e[2])
# ymn<-floor(e[3])
# ymx<-ceiling(e[4])

#e2<-extent(xmn,xmx,ymn,ymx)

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

#linking up the sdm map cells with the population cells
sdm[is.na(sdm)] <- 0
sdm_loss[is.na(sdm_loss)] <- 0

sdm_table<-as.data.frame(sdm, xy = TRUE)
sdm_loss_table<-as.data.frame(sdm_loss, xy = TRUE)

sdm_table$ID<-1:length(sdm_table$ID)
sdm_table$ID[which(!is.na(df$PatchID))]<-df$PatchID[!is.na(df$PatchID)]

niche_map_mine<-cbind(sdm_table[,c(4,1,2,3)],sdm_loss_table$Average_SDM_2016)
colnames(niche_map_mine)<-c("gridID", "X", "Y", "Year_1950", "Year_2016")

niche_formulas <- as.formula(paste(paste(colnames(niche_map_mine)[-c(1:3)], collapse="+"),"X+Y",sep="~"))

print(levelplot(niche_formulas, niche_map_mine, col.regions=rev(heat.colors(100)), main = "Niche Values"))

no_yrs_mine<-33 #number of years each time period represents - could be 1 when I do it for reals

####matrix set up
prob_scenario<-c(0.5,0.5)    #need to check this

noise<-0.95     #need to check this

stages<-comadre$matrixClass[keep][[1]]$MatrixClassAuthor
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

density_individuals <- 1   #also compulsory not sure what best value would be

K<-1000   #carrying capacity
 
K_weight<-c(0, 1.5,1,1,1,1,1,1)  #the weight with which carrying capacity affects each stage was FALSE

fraction_SDD <- 0.05  #short distance dispersal

dispersal_constants_mine <-c(0.7,0.7,0.1,3)

#fraction_LDD_mine <- 0.05 #long distance dispersal
##Populations - prioritise sorting this out
  
niche_map_mine2<-c("period_A", "period_B")  


demoniche_setup(modelname = "Capra_ibex",Populations = Populations, Nichemap = niche_map_mine,
                matrices = matrices,matrices_var = matrices_var, prob_scenario = prob_scenario,
                stages = stages, proportion_initial = proportion_initial,
                density_individuals = density_individuals,
                fraction_LDD = 0.5, fraction_SDD = 0.5,
                dispersal_constants = dispersal_constants_mine,
                transition_affected_niche = transition_affected_niche,
                transition_affected_demogr = transition_affected_demogr,
                transition_affected_env=transition_affected_env,
                env_stochas_type = env_stochas_type,
                no_yrs = no_yrs_mine, K=K, Kweight = K_weight, 
                sumweight =sumweight)


#works and is simple - exponential increase
#demoniche_setup(modelname = "Capra_ibex",Populations = Populations, 
#                matrices_var = matrices_var,matrices = matrices,
#                 stages = stages, proportion_initial = proportion_initial,
#                density_individuals = density_individuals,
#                no_yrs = 80, sumweight =sumweight)


#important to include sumweight, I think the default is FALSE but that 
#causes the population to be 0 in all years

RPyran_min_run <- demoniche_model(modelname = "Capra_ibex", Niche = TRUE, 
                                  Dispersal = FALSE, repetitions = 1,
                                  foldername = "RPyran_minimal")

RPyran_min_run[,"Meanpop","Reference_matrix"]

years<-1950:2016
plot(RPyran_min_run[,"Meanpop","Reference_matrix"])

bleh<-cbind(years,RPyran_min_run[,"Meanpop","Reference_matrix"])
plot(bleh)



























