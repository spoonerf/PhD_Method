wd<-getwd()

.libPaths(c(wd,.libPaths()))

library(raster)
library(demoniche)
library(doParallel)
library(sp)
library(popbio)

source("demoniche_setup_me.R")
source("demoniche_model_me.R")
source("demoniche_dispersal_me.R")
source("demoniche_population_function.R")

wd<-getwd()
genus<-"Cervus"
species<-"elaphus"
Europe_only<-TRUE
binomial<-paste(genus, species, sep="_")
min_lat<- 35    #one degree more/less than the iucn range extent
max_lat<- 68
min_lon<- -11
max_lon<- 33
bioclim_layers<-c(2,4,5,6,7,10,11)   
no_background_points<-1000
bioclim_names<-c("Bio_2_2006_2016_average", "Bio_4_2006_2016_average","Bio_5_2006_2016_average","Bio_6_2006_2016_average","Bio_7_2006_2016_average","Bio_10_2006_2016_average","Bio_11_2006_2016_average")


#Comadre
transition_affected_niche<-"all"  #which parts of the matrix are affected by the c(1,2) is juveniles
#niche values
transition_affected_env <- "all"
transition_affected_demogr <- "all"


SDD_seq<-c(0.1, 0.25, 0.5)
LDD_seq<-c(0.1, 0.25, 0.5)

SD<-c(0.1,0.25, 0.5)

kern_seq<-list(c(49.7,284.2),c(13.7, 89.4), c(1000,1000))   #first values is median derived from median home range and second value is maximum dispersal distance

density_mine<- c(500, 1000,1564)  
carry_k<-c(1000, 1617, 2000)

var_grid<-expand.grid(SD,SDD_seq, LDD_seq, 1:length(kern_seq), density_mine)
colnames(var_grid)<-c("SD","SDD" ,"LDD", "Kern", "Density")


binomial<-paste(genus, species, sep="_")

years<-1950:2005 #can only go up to 2005 with hyde
spin_years<-1750:1949
lpi<-read.csv("LPI_pops_20160523_edited.csv")
species_directory<-paste(wd, binomial, sep="/")
dir.create(species_directory)
sdm_folder<-paste(species_directory, "SDM_folder", sep = "/")
dir.create(sdm_folder)
demoniche_folder<-paste(species_directory, "Demoniche_Output", sep = "/")
dir.create(demoniche_folder)

no_yrs_mine<-1
prob_scenario<-c(0.5,0.5)    #need to check this
noise<-0.90

env_stochas_type<-"normal"  


pyr<-subset(lpi, Binomial ==binomial & Specific_location==1)    #record 11470 had wrong longitude - in Russia!

#formatting the lpi data for use in demoniche
pyrs<-pyr[,c("ID","Longitude","Latitude")]
pyrs$ID<-pyrs$ID * 100
xy_lpi<-data.frame(pyrs$Longitude, pyrs$Latitude)

xy<-read.csv(paste(wd, "/", binomial, "_locs.csv", sep=""))

#formatting the gbif data for use in demoniche

gbif_xy<-data.frame(xy$V1, xy$V2)
coordinates(gbif_xy)<-c("xy.V1","xy.V2")
proj4string(gbif_xy)<- CRS("+init=epsg:4326")
patch<-raster(paste(wd, "/hyde_pres_abs_sss_weighted_ensemble_sdm_1950.tif", sep=""))

gbif_xy<-raster:::rasterize(gbif_xy, patch)
gbif_xy<-as.data.frame(gbif_xy, xy=TRUE)
gbif_xy<-na.omit(gbif_xy)
gbif_xy$ID<-(1:nrow(gbif_xy))+ncell(patch)*2
gbif_xy<-gbif_xy[,c("ID", "x", "y")]
colnames(gbif_xy)<-c("ID", "Longitude", "Latitude")

# coordinates(xy_lpi)<-c("pyrs.Longitude", "pyrs.Latitude")
# proj4string(xy_lpi) <- CRS("+init=epsg:4326") # WGS 84
# xy_lpi<-data.frame(xy_lpi)


Populations_lpi<-data.frame(pyrs$ID, xy_lpi$pyrs.Longitude, xy_lpi$pyrs.Latitude)
colnames(Populations_lpi)<-c("PatchID", "X", "Y")

Populations_gbif<-data.frame(gbif_xy$ID, gbif_xy$Longitude, gbif_xy$Latitude)
colnames(Populations_gbif)<-c("PatchID", "X", "Y")

Populations<-rbind(Populations_lpi, Populations_gbif)

Populations$area<-1
pop_years<-pyr[,65:130]
pop_years[pop_years == "NULL"]<-NA

first_year<-function(x){
  pop_first<-min(which(!is.na(x)))
  pop_value<-x[pop_first]
  return(pop_value)
}
pop_values<-apply(pop_years, 1, first_year)

density_mine<-as.numeric(pop_values)


id<-pyrs$ID*100
lam<-rep(1,length(id))    #not sure what the value here pertains to - think it sets starting population so should use values from LPI?
pyrxy<-SpatialPoints(pyr[,c("Longitude","Latitude")])
sdm<-raster(paste(sdm_folder, "/pres_abs_sss_weighted_ensemble_sdm_1950.tif", sep=""))
e2<-extent(sdm)
r<-raster(e2, resolution=res(sdm))
rz<-rasterize(pyrxy,r,lam )
crs(rz)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
rid<-rasterize(pyrxy,r,id)
crs(rid)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
rz_spdf<-xyFromCell(rz, 1:ncell(rid))
rzm<-as.vector(rz)
ridm<-as.vector(rid)
df<-data.frame(ridm,rz_spdf,rzm)
colnames(df)<-c( "PatchID","X","Y","area")
#df<-data.frame(na.omit(df))
#Populations<-data.frame(na.omit(df))


###formatting environmental data
patch<-raster(paste(sdm_folder,"/hyde_pres_abs_sss_weighted_ensemble_sdm_1950.tif", sep=""))

sdm_patch_df<-data.frame(ID=1:ncell(patch))
sdm_df<-data.frame(ID=1:ncell(patch))

for (i in 1:length(years)){
  sdm<-raster(paste(sdm_folder,"/hyde_weighted_ensemble_sdm_", years[i],".tif", sep=""))   #14.6
  patch<-raster(paste(sdm_folder,"/hyde_pres_abs_sss_weighted_ensemble_sdm_", years[i],".tif", sep=""))
  
  if (i ==1){
    vec<-as.data.frame(sdm, xy = TRUE)
    vec_pat<-as.data.frame(patch, xy=TRUE)
  } else{
    vec<-as.data.frame(sdm)
    
    vec_pat<-as.data.frame(patch)
  }
  
  sdm_df<-cbind(sdm_df, vec)
  sdm_patch_df<-cbind(sdm_patch_df,vec_pat)
}

sdm_df$ID[which(!is.na(df$PatchID))]<-df$PatchID[!is.na(df$PatchID)]
sdm_patch_df$ID[which(!is.na(df$PatchID))]<-df$PatchID[!is.na(df$PatchID)]

niche_map_mine<-sdm_df
colnames(niche_map_mine)[1:3]<-c("gridID", "X", "Y")

patch_map_mine<-sdm_patch_df
colnames(patch_map_mine)[1:3]<-c("gridID", "X", "Y")

niche_spin_up<-matrix(rep(niche_map_mine[,4], length(spin_years)), nrow=nrow(niche_map_mine))
niche_spin_up<-cbind(niche_map_mine[,1:3],niche_spin_up, niche_map_mine[,4:ncol(niche_map_mine)])

patch_spin_up<-matrix(rep(patch_map_mine[,4], length(spin_years)), nrow=nrow(patch_map_mine))
patch_spin_up<-cbind(patch_map_mine[,1:3],patch_spin_up, patch_map_mine[,4:ncol(patch_map_mine)]) #last patch used to be niche

col_years_short<-paste("Year_", years, sep="")
col_years<-paste("Year_", min(spin_years):max(years), sep="")

colnames(niche_map_mine)[4:length(colnames(niche_map_mine))]<-col_years_short
colnames(patch_map_mine)[4:length(colnames(patch_map_mine))]<-col_years_short

colnames(niche_spin_up)[4:length(colnames(niche_spin_up))]<-col_years
colnames(patch_spin_up)[4:length(colnames(patch_spin_up))]<-col_years

opt_patch_spin_up<-patch_spin_up
opt_patch_spin_up[,4:704]<-1

patch_map_mine[is.na(patch_map_mine)]<-0

niche_spin_up[is.na(niche_spin_up)]<-0
patch_spin_up[is.na(patch_spin_up)]<-0

###COMADRE

load(paste(wd, "COMADRE_v.2.0.1.RData", sep="/"))

tempMetadata<-subset(comadre$metadata, SpeciesAccepted==binomial)

keep<-as.numeric(rownames(tempMetadata))

tempMat<-comadre$mat[keep]   #MatA is matrix pop model, can be split into U, F and/or C

MatList<-list(tempMat[[2]][[1]])  #varies depending on number of matrices - need to find a way to code this better - now have five matrices available so need to sort this
AllMat<-unlist(MatList)
matrices<-matrix(AllMat, ncol=length(MatList))
colnames(matrices)<- c("Reference_matrix")


MatListA<-list(tempMat[[3]][[1]])  #varies depending on number of matrices - need to find a way to code this better - now have five matrices available so need to sort this
AllMatA<-unlist(MatListA)
matricesA<-matrix(AllMatA, ncol=length(MatListA))
colnames(matricesA)<- c("Matrix A")

MatListB<-list(tempMat[[4]][[1]])  #varies depending on number of matrices - need to find a way to code this better - now have five matrices available so need to sort this
AllMatB<-unlist(MatListB)
matricesB<-matrix(AllMatB, ncol=length(MatListB))
colnames(matricesB)<- c("Matrix B")


matrices<-cbind(matrices, matricesA, matricesB)


#not yet active in demoniche

stages<-comadre$matrixClass[keep][[2]]$MatrixClassAuthor
#stages<-comadre$matrixClass[keep][[2]]$MatrixClassAuthor
#stagesf<-stages[1:3]
proportion_initial<- rep(1/length(stages), length(stages)) #Think spin up sorts this
sumweight<-rep(1, length(stages))#weight of stages  - should be equal for all mine just
#in plants seed not included in calculating population sizes
list_names_matrices<-colnames(matrices)
K_weight<-c(rep(1, length(stages)))  #the weight with which carrying capacity affects each stage was FALSE


####Scaling Carrying Capacity

lin<-function(x, carry_k){
  x*carry_k
}

lf<-list.files(sdm_folder)

files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]

sdms<-stack(paste(files, sep="/"))

hsi<-raster:::extract(sdms, Populations[,c(2,3)])

link1<-lin(hsi, carry_k[1])
link2<-lin(hsi, carry_k[2])
link3<-lin(hsi, carry_k[3])

spin1<-replicate(length(spin_years),link1[,1])
link_spin1<-cbind(spin1, link1)
colnames(link_spin1)[1:length(spin_years)]<-paste("hyde_weighted_ensemble_sdm_", spin_years, sep="")

spin2<-replicate(length(spin_years),link2[,1])
link_spin2<-cbind(spin2, link2)
colnames(link_spin2)[1:length(spin_years)]<-paste("hyde_weighted_ensemble_sdm_", spin_years, sep="")

spin3<-replicate(length(spin_years),link3[,1])
link_spin3<-cbind(spin3, link3)
colnames(link_spin3)[1:length(spin_years)]<-paste("hyde_weighted_ensemble_sdm_", spin_years, sep="")

link_spin<-list(link_spin1, link_spin2, link_spin3)



###Running Demoniche Model

reps<-12

for (s in 1:nrow(var_grid)){
  
  print(paste (s, " out of ", nrow(var_grid) ), sep="")
  
  SDD<-var_grid[s,"SDD"]
  LDD<-var_grid[s,"LDD"]
  kern<-var_grid[s, "Kern"]
  SD<-var_grid[s, "SD"]
  SD<-0.1
  
  K<-var_grid[s, "K"]
  dens<-var_grid[s, "Density"]
  link_id<-var_grid[s, "K_scale"]
  
  matrices_var<-matrix(SD, ncol = 1, nrow = nrow(matrices), dimnames = list(NULL, "sd"))
  
  start.time <- Sys.time()
  
  med_disp<-as.character(kern_seq[[kern]][1])
  dir.create(paste(demoniche_folder,"/hyde_new_patch_disp_test_",med_disp,"_",SDD,"_",LDD,"_",kern,"_",SD,"_",K,"_",dens,"_",link_id, "/",sep=""),showWarnings = TRUE)
  
  rep_demoniche<-function(i){
    wd<-getwd()
    
    .libPaths(c(wd,.libPaths()))
    library(demoniche)
    library(doParallel)
    source("demoniche_setup_me.R")
    source("demoniche_model_me.R")
    source("demoniche_dispersal_me.R")
    demoniche_setup_me(modelname = binomial ,Populations = Populations, Nichemap = patch_spin_up,
                       matrices = matrices,matrices_var = matrices_var, prob_scenario = prob_scenario,
                       stages = stages, proportion_initial = proportion_initial,
                       density_individuals = dens,
                       fraction_LDD = LDD, fraction_SDD = SDD,
                       dispersal_constants = kern_seq[[kern]],
                       transition_affected_niche = transition_affected_niche,
                       transition_affected_demogr = transition_affected_demogr,
                       transition_affected_env=transition_affected_env,
                       env_stochas_type = env_stochas_type,
                       no_yrs = no_yrs_mine, K=link_spin[link_id], Kweight = K_weight, Ktype="ceiling",
                       sumweight =sumweight)
    
    
    c_ibex_k_16000 <- demoniche_model_me(modelname = binomial, Niche = TRUE,
                                         Dispersal = TRUE, repetitions = 1,
                                         foldername = paste(binomial,"/Demoniche_Output/hyde_new_patch_disp_test_",med_disp,"_",SDD,"_",LDD,"_",kern,"_",SD,"_",K,"_",dens,"_",link_id,"/",i, sep = ""))
  }
  
  if (Sys.info()["nodename"] == "FIONA-PC"){
    cl <- makeCluster(4)
  } else {
    cl <- makeCluster(16)
  }
  
  registerDoParallel(cl)
  foreach(i=(1:reps)) %dopar% rep_demoniche(i)
  stopCluster(cl)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}





