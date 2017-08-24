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
library(dismo)
library(demoniche)
library(biomod2)

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

years<-as.character(1950:2016)
e<-extent(sp)+4
species<-"Capra ibex"

myBiomodOption <- BIOMOD_ModelingOptions()


for (i in 1:length(years)){
  
  bio<-stack(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Bioclim/",years[i],"_bioclim_variable_stack.tif", sep=""))  
  alps<-crop(bio,e)  
  alps<-stack(alps)
  my_biomod_data<-BIOMOD_FormatingData(resp.var = sp,
                                       expl.var =  alps,  #alps for 1950, alps now for 2016
                                       #eval.resp.var = sp_eval,
                                       #eval.expl.var = bio_1950,
                                       resp.name = "Capra_ibex",
                                       PA.nb.rep = 10,    #switch these three off to run without pseudo absences
                                       PA.nb.absences = 10,
                                       PA.strategy = 'random',
                                       na.rm=TRUE)
  
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
  
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = alps,    
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd')
  
  myCurrentProj <- get_predictions(myBiomodProj)
  
  ave_proj<-mean(myCurrentProj)
  writeRaster(ave_proj, paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Alp_SDMs/", years[i],"_SDM.tif", sep=""), overwrite=TRUE)
  print(years[i])
}

#############################

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
plot(ann_sum)

plot(pyr$Longitude, pyr$Latitude)

pyrs<-pyr[,c("ID","Longitude","Latitude")]

id<-pyrs$ID*100
lam<-rep(200,length(id))    #not sure what the value here pertains to - think it sets starting population so should use values from LPI?
pyrxy<-SpatialPoints(pyr[,c("Longitude","Latitude")])

sdm<-raster("C:/Users/Fiona/Documents/PhD/PhD_Method/Alp_SDMs/1950_SDM.tif")
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


sdm_df<-data.frame(ID=1:2030)

#formatting data for demoniche
for (i in 1:length(years)){
  
  sdm<-raster(paste("C:/Users/Fiona/Documents/PhD/PhD_Method/Alp_SDMs/", years[i],"_SDM.tif", sep=""))  
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



prob_scenario<-c(0.5,0.5)    #need to check this

noise<-0.90     #need to check this

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

density_individuals <- 4292.32   #4292.32 to 16096.2 based on density being between 8 and 30 per 100 ha and the area of each cell being 53654 ha 

K<-NULL   #carrying capacity

K_weight<-c(1, 1,1,1,1,1,1,1)  #the weight with which carrying capacity affects each stage was FALSE

fraction_SDD <- 0.05  #short distance dispersal

dispersal_constants_mine <-c(0.7,0.7,0.1,3)

niche_map_mine[is.na(niche_map_mine)]<-0

demoniche_setup(modelname = "Capra_ibex",Populations = Populations, Nichemap = niche_map_mine,
                matrices = matrices,matrices_var = matrices_var, prob_scenario = prob_scenario,
                stages = stages, proportion_initial = proportion_initial,
                density_individuals = density_individuals,
                fraction_LDD = 0.2, fraction_SDD = 0.2,
                dispersal_constants = dispersal_constants_mine,
                transition_affected_niche = transition_affected_niche,
                transition_affected_demogr = transition_affected_demogr,
                transition_affected_env=transition_affected_env,
                env_stochas_type = env_stochas_type,
                no_yrs = no_yrs_mine, K=K, Kweight = K_weight, 
                sumweight =sumweight)

RPyran_niche <- demoniche_model(modelname = "Capra_ibex", Niche = TRUE, 
                                Dispersal = FALSE, repetitions = 1,
                                foldername = "RPyran_minimal")

RPyran_niche[,"Meanpop","Reference_matrix"]

years<-1950:2016
plot(RPyran_niche[,"Meanpop","Reference_matrix"])

bleh<-cbind(years,RPyran_niche[,"Meanpop","Reference_matrix"])
plot(bleh, type="l")

RPyran_disp_niche <- demoniche_model(modelname = "Capra_ibex", Niche = TRUE, 
                                     Dispersal = TRUE, repetitions = 1,
                                     foldername = "RPyran_disp_niche")

RPyran_disp_niche[,"Meanpop","Reference_matrix"]

years<-1950:2016
#plot(RPyran_disp_niche[,"Meanpop","Reference_matrix"])

bleh<-cbind(years,RPyran_disp_niche[,"Meanpop","Reference_matrix"], RPyran_min_run[,"Max","Reference_matrix"])
plot(bleh, type="l", col="red")


























