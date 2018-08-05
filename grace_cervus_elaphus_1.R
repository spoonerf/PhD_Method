
#running with only one matrix

wd<-getwd()

.libPaths(c(wd,.libPaths()))

library(raster)
#library(demoniche)
library(doParallel)
library(sp)
library(popbio)

source("demoniche_setup_csv.R")
source("demoniche_model_csv.R")
source("demoniche_dispersal_csv.R")
source("demoniche_population_function.R")
source("randomLHS.R")

genus<-"Cervus"
species<-"elaphus"
Europe_only<-TRUE
binomial<-paste(genus, species, sep="_")
min_lat<- 35    #one degree more/less than the iucn range extent
max_lat<- 68
min_lon<- -11
max_lon<- 33

e_iucn<-extent(min_lon, max_lon, min_lat, max_lat)
bioclim_layers<-c(2,4,5,6,7,10,11)
no_background_points<-1000
bioclim_names<-paste("Bio_", bioclim_layers, "_2006_2016_average", sep="")

transition_affected_niche<-"all"  #which parts of the matrix are affected by the c(1,2) is juveniles
#niche values
transition_affected_env <- "all"
transition_affected_demogr <- "all"

env_stochas_type<-"normal" 

tetra_dens<-read.csv("TetraDENSITY.csv")
tetra_dens$binomial<-paste(tetra_dens$Genus, tetra_dens$Species, sep="_")

bin_dens<-tetra_dens[tetra_dens$binomial == binomial,]
qexp_rate<-1/mean(bin_dens$Density)

med_disp_l<-13.7
med_disp_h<-49.7

max_disp_l<-89.4
max_disp_h<-284.2

carry_k_l<-25
carry_k_h<-33

SDD<-0

n<-6

var_grid <- randomLHS(round(15.6*n),n)

var_grid[,1]<-qunif(var_grid[,1],0, 0.25) #SD
var_grid[,2]<-qunif(var_grid[,2],0, 0.25) #LDD
var_grid[,3]<-qexp(var_grid[,3], qexp_rate)#dens
var_grid[,4]<-qunif(var_grid[,4],carry_k_l, carry_k_h)#carrying capacity
var_grid[,5]<-qunif(var_grid[,5],med_disp_l, med_disp_h)#med disp
var_grid[,6]<-qunif(var_grid[,6],max_disp_l, max_disp_h)#max disp

var_grid<-cbind(var_grid, SDD)

var_grid<-round(var_grid, 2)

colnames(var_grid)<-c("SD","LDD", "Density", "Carry_K", "Med_Disp", "Max_Disp","SDD" )

var_grid<-read.csv("var_grid.csv")

var_grid<-var_grid[,-1]

colnames(var_grid)<-c("SD","LDD", "Density", "Carry_K", "Med_Disp", "Max_Disp","SDD" )

binomial<-paste(genus, species, sep="_")

years<-1950:2005 #can only go up to 2005 with hyde
spin_years<-1850:1949
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


env_stochas_type<-"normal"   #can also be lognormal

###############

xy<-read.csv(paste(binomial, "_locs.csv", sep=""))
patch<-stack(paste(wd, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years[1],".tif", sep=""))
patch_stack<-stack(paste(wd, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years,".tif", sep=""))

patch<-crop(patch, e_iucn)
patch_stack<-crop(patch_stack, e_iucn)

patch_sum<-raster::calc(patch_stack, sum)
values(patch_sum)[(values(patch_sum)>=1)]<-1

test_patch<-patch
values(test_patch)[values(test_patch==0) & !is.na(values(test_patch))]<-1
values(test_patch)[is.na(values(test_patch))]<-0

pyr<-subset(lpi, Binomial ==binomial & Specific_location==1)    #record 11470 had wrong longitude - in Russia!

#formatting the lpi data for use in demoniche
pyrs<-pyr[,c("ID","Longitude","Latitude")]
pyrs$ID<-pyrs$ID * 100
xy_lpi<-data.frame(pyrs$Longitude, pyrs$Latitude)

xy<-read.csv(paste(binomial, "_locs.csv", sep=""))
patch<-stack(paste(wd, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years[1],".tif", sep=""))
patch_stack<-stack(paste(wd, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years,".tif", sep=""))

patch<-crop(patch, e_iucn)
patch_stack<-crop(patch_stack, e_iucn)

patch_sum<-raster::calc(patch_stack, sum)
values(patch_sum)[(values(patch_sum)>=1)]<-1

start_pops<-which(values(!is.na(patch_sum)))
start_pops<-sample(start_pops, 1000, replace = FALSE)
sdm_pops<-xyFromCell(patch_sum, start_pops)

ID<-1:length(sdm_pops[,1])*100
sdm_pops<-cbind(ID, sdm_pops)

#formatting the lpi data for use in demoniche

pyrs$ID<-pyrs$ID * 100
pyrs<-pyrs[pyrs$Longitude > e_iucn[1] & pyrs$Longitude < e_iucn[2] & pyrs$Latitude > e_iucn[3] & pyrs$Latitude < e_iucn[4],]
xy_lpi<-data.frame(pyrs$Longitude, pyrs$Latitude)

coordinates(xy_lpi)<-c("pyrs.Longitude", "pyrs.Latitude")
proj4string(xy_lpi) <- CRS("+init=epsg:4326") # WGS 84
lpi_xy<-raster:::rasterize(xy_lpi, patch)
values(lpi_xy)[values(!is.na(lpi_xy))]<-1
values(lpi_xy)[values(is.na(lpi_xy))]<-0

lpi_pops<-which(values(lpi_xy) ==1)
lpi_pops<-xyFromCell(lpi_xy, lpi_pops)

lpi_pops<-cbind(pyrs$ID, lpi_pops)

overlap<-which(sdm_pops[,2] %in% lpi_pops[,2] & sdm_pops[,3] %in% lpi_pops[,3])

if (length(overlap)>0){
  
  sdm_pops<-sdm_pops[-overlap,]
}

Populations<-rbind(sdm_pops, lpi_pops)
Populations<-data.frame(Populations)

colnames(Populations)<-c("PatchID", "X", "Y")
Populations<-Populations[Populations$X > e_iucn[1] &  Populations$X < e_iucn[2] & Populations$Y > e_iucn[3]& Populations$Y < e_iucn[4],]

pxy<-cbind(Populations$X, Populations$Y)
a<-area(patch)

Populations$area<-extract(a, pxy)




###formatting environmental data
patch<-raster(paste("hyde_pres_abs_sss_weighted_ensemble_sdm_1950.tif", sep=""))
patch<-crop(patch, e_iucn)

sdm_patch_df<-data.frame(ID=1:ncell(patch))
sdm_df<-data.frame(ID=1:ncell(patch))

for (i in 1:length(years)){
  sdm<-raster(paste("hyde_weighted_ensemble_sdm_", years[i],".tif", sep=""))   #14.6
  patch<-raster(paste("hyde_pres_abs_sss_weighted_ensemble_sdm_", years[i],".tif", sep=""))
  
  sdm<-crop(sdm, e_iucn)
  patch<-crop(patch, e_iucn)
  
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

niche_map_mine<-sdm_df
colnames(niche_map_mine)[1:3]<-c("gridID", "X", "Y")

patch_map_mine<-sdm_patch_df
colnames(patch_map_mine)[1:3]<-c("gridID", "X", "Y")

niche_spin_up<-matrix(rep(rowMeans(niche_map_mine[,4:ncol(niche_map_mine)]), length(spin_years)), nrow=nrow(niche_map_mine))
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
opt_niche_spin_up<-niche_spin_up

opt_patch_spin_up[is.na(opt_patch_spin_up)]<--1
opt_niche_spin_up[is.na(opt_niche_spin_up)]<--1

opt_patch_spin_up[,4:(length(spin_years)+3)][opt_patch_spin_up[,4:(length(spin_years)+3)]>=0]<-1
opt_niche_spin_up[,4:(length(spin_years)+3)][opt_niche_spin_up[,4:(length(spin_years)+3)]>=0]<-1

av_hab<-mean(na.omit(niche_map_mine[,4:ncol(niche_map_mine)][niche_map_mine[,4:ncol(niche_map_mine)] * patch_map_mine[,4:ncol(niche_map_mine)]  > 0]))
max_hab<-max(na.omit(niche_map_mine[,4:ncol(niche_map_mine)][niche_map_mine[,4:ncol(niche_map_mine)] * patch_map_mine[,4:ncol(niche_map_mine)]  > 0]))

patch_spin_up[is.na(patch_spin_up)]<--1
niche_spin_up[is.na(niche_spin_up)]<--1

###COMADRE
load(paste(wd, "COMADRE_v.2.0.1.RData", sep="/"))

tempMetadata<-subset(comadre$metadata, SpeciesAccepted==binomial)

keep<-as.numeric(rownames(tempMetadata))

tempMat<-comadre$mat[keep] 

A<-tempMat[[1]][[1]]
AS<-A*1/av_hab
SM<-AS   ##capping survival at 1
# SM[SM>1]<-1/max_hab
# AS[-1, ] <- SM[-1, ]

MatList_S<-list(AS)  #varies depending on number of matrices - need to find a way to code this better - now have five matrices available so need to sort this
AllMat_S<-unlist(MatList_S)
matrices_scale<-matrix(AllMat_S, ncol=length(MatList_S))
colnames(matrices_scale)<- c("Reference_matrix")

matrices_scale<-cbind(matrices_scale)

MatList<-list(tempMat[[2]][[1]])  #varies depending on number of matrices - need to find a way to code this better - now have five matrices available so need to sort this
AllMat<-unlist(MatList)
matrices<-matrix(AllMat, ncol=length(MatList))
colnames(matrices)<- c("Reference_matrix")

matrices<-cbind(matrices)

#not yet active in demoniche

stages<-comadre$matrixClass[keep][[1]]$MatrixClassAuthor
proportion_initial<- rep(1/length(stages), length(stages)) #Think spin up sorts this
sumweight<-rep(1, length(stages))#weight of stages  - should be equal for all mine just
#in plants seed not included in calculating population sizes
list_names_matrices<-colnames(matrices)
K_weight<-c(rep(1, length(stages)))  #the weight with which carrying capacity affects each stage was FALSE
#######################



#######################

####Scaling Carrying Capacity


lin<-function(x, carry_k){
  x*carry_k
}

#################

lf<-list.files(wd)

files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]

sdms<-stack(paste(wd, files, sep="/"))
sdms<-crop(sdms, e_iucn)

hsi<-raster:::extract(sdms, Populations[,c(2,3)])

###Running Demoniche Model
dispersal_probabilities<-paste(wd, "/disp_probs.csv", sep="")

#dist_populations<-nc_open("dist_pops_all.nc")

reps<-100

#for (s in 1:nrow(var_grid)){
for (s in 11:20){
  
  print(paste (s, " out of ", nrow(var_grid) ), sep="")
  
  SDD<-var_grid[s,"SDD"]
  LDD<-var_grid[s,"LDD"]
  SD<-var_grid[s, "SD"]
  dens<-var_grid[s, "Density"]
  K<-var_grid[s, "Carry_K"]
  med_disp<-as.numeric(var_grid[s,"Med_Disp"])
  max_disp<-as.numeric(var_grid[s,"Max_Disp"])
  
  link<-lin(hsi, K)
  link<-link*Populations$area
  spin_k<-lin(1, K)*Populations$area
  spin<-replicate(length(spin_years),spin_k)
  link_spin<-cbind(spin, link)
  colnames(link_spin)[1:length(spin_years)]<-paste("hyde_weighted_ensemble_sdm_", spin_years, sep="")
  
  matrices_var<-matrix(SD, ncol = 1, nrow = nrow(matrices), dimnames = list(NULL, "sd"))
  
  start.time <- Sys.time()
  
  dir.create(paste(demoniche_folder,"/hyde_new_patch_disp_test_",SD,"_",SDD,"_",LDD,"_",dens,"_",K,"_",med_disp,"_",max_disp, "/",sep=""),showWarnings = TRUE)
  
  rep_demoniche<-function(i){
    wd<-getwd()
    
    .libPaths(c(wd,.libPaths()))
    library(demoniche)
    library(doParallel)
    source("demoniche_setup_csv.R")
    source("demoniche_model_csv.R")
    source("demoniche_dispersal_csv.R")
    source("demoniche_population_function.R")
    demoniche_setup_csv(modelname = binomial ,Populations = Populations, Nichemap = opt_niche_spin_up,
                       matrices = matrices_scale,matrices_var = matrices_var, prob_scenario = prob_scenario,
                       stages = stages, proportion_initial = proportion_initial,
                       density_individuals = dens,
                       fraction_LDD = LDD, fraction_SDD = SDD,
                       dispersal_constants = c(med_disp, max_disp),
                       transition_affected_niche = transition_affected_niche,
                       transition_affected_demogr = transition_affected_demogr,
                       transition_affected_env=transition_affected_env,
                       env_stochas_type = env_stochas_type,
                       no_yrs = no_yrs_mine, K=link_spin, Kweight = K_weight, Ktype="ceiling",
                       sumweight =sumweight, dispersal_probabilities = dispersal_probabilities)
    
    
    c_ibex_k_16000 <- demoniche_model_csv(modelname = binomial, Niche = TRUE,
                                         Dispersal = TRUE, repetitions = 1,
                                         foldername = paste(binomial,"/Demoniche_Output/hyde_new_patch_disp_test_",SD,"_",SDD,"_",LDD,"_",dens,"_",K,"_",med_disp,"_",max_disp, "/",i,sep=""))
    
    
    files<-lf[grepl("*.rda", lf)]
    file.remove(files)
    

    
  }
  
  if (Sys.info()["nodename"] == "FIONA-PC"){
    cl <- makeCluster(4)
  } else {
    cl <- makeCluster(100)
  }
  
  registerDoParallel(cl)
  foreach(i=(1:reps)) %dopar% rep_demoniche(i)
  stopCluster(cl)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}


