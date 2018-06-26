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
source("C:/Users/Fiona/Documents/PhD/PhD_Method/demoniche_population_function_optim.R")

genus<-"Capra"
species<-"ibex"
Europe_only<-TRUE
binomial<-paste(genus, species, sep="_")
min_lat<- 43
max_lat<- 47.9
min_lon<- 0
max_lon<- 16
bioclim_layers<-c(1,5,6,13,15,18,19)
no_background_points<-1000
bioclim_names<-paste("Bio_",bioclim_layers,"_2006_2016_average", sep="")

transition_affected_niche<-"all"  #which parts of the matrix are affected by the c(2,11,20) is  juveniles
transition_affected_env <- "all"
transition_affected_demogr <- "all"

env_stochas_type<-"normal"   #can also be lognormal

density_mine <- c(8,30) #per km2  #4292.32 to 16096.2 based on density being between 8 and 30 per 100 ha and the area of each cell being 53654 ha -
carry_k<-c(30) #per km2

#Things to vary
SDD_seq<-c(0.1, 0.25, 0.5)
LDD_seq<-c(0.1, 0.25, 0.5)

SD<-c(0.1,0.25, 0.5)

med_disp<-seq(13.5975, 19.775, length.out = 3)

kern_seq<-list(c(med_disp[1], 113),c(med_disp[2],113.000),c(med_disp[3],113.000))
var_grid<-expand.grid(SDD_seq, LDD_seq, 1:length(kern_seq), SD, carry_k, density_mine, 1:length(carry_k))
colnames(var_grid)<-c("SDD", "LDD", "Kern", "SD", "K", "Density", "K_scale")

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

xy<-read.csv(paste(binomial, "_locs.csv", sep=""))
patch<-raster(paste(wd, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years[1],".tif", sep=""))
patch_stack<-stack(paste(wd, "/hyde_pres_abs_sss_weighted_ensemble_sdm_", years,".tif", sep=""))
patch_sum<-raster::calc(patch_stack, sum)
values(patch_sum)[(values(patch_sum)>=1)]<-1

test_patch<-patch
values(test_patch)[values(test_patch==0) & !is.na(values(test_patch))]<-1
values(test_patch)[is.na(values(test_patch))]<-0

pyr<-subset(lpi, Binomial ==binomial & Specific_location==1)  

#formatting the lpi data for use in demoniche
pyrs<-pyr[,c("ID","Longitude","Latitude")]
pyrs$ID<-pyrs$ID * 100
xy_lpi<-data.frame(pyrs$Longitude, pyrs$Latitude)


#formatting the gbif data for use in demoniche

gbif_xy<-data.frame(xy$V1, xy$V2)
coordinates(gbif_xy)<-c("xy.V1","xy.V2")
proj4string(gbif_xy)<- CRS("+init=epsg:4326")
gbif_xy<-raster:::rasterize(gbif_xy, patch)
values(gbif_xy)[values(!is.na(gbif_xy))]<-1

coordinates(xy_lpi)<-c("pyrs.Longitude", "pyrs.Latitude")
proj4string(xy_lpi) <- CRS("+init=epsg:4326") # WGS 84
lpi_xy<-raster:::rasterize(xy_lpi, patch)
values(lpi_xy)[values(!is.na(lpi_xy))]<-1
values(lpi_xy)[values(is.na(lpi_xy))]<-0


gbif_xy2<-gbif_xy - lpi_xy
#values(gbif_xy2)[values(gbif_xy2==0)]<-NA
odd<-(test_patch +gbif_xy2)
values(odd)[values(odd) <2 | is.na(values(odd))]<-NA
values(odd)[values(odd) ==2]<-1

#removing populations from any cells which are unsuitable in all years (0/1)
pop_suit<-odd+patch_sum
values(pop_suit)[values(pop_suit) ==1]<-NA
values(pop_suit)[values(pop_suit) ==2]<-1

gbif_xy<-as.data.frame(pop_suit, xy=TRUE)
gbif_xy<-na.omit(gbif_xy)
gbif_xy$ID<-(1:nrow(gbif_xy))+ncell(patch)*2
Populations<-gbif_xy[,c("ID", "x", "y")]
colnames(Populations)<-c("PatchID", "X", "Y")

pxy<-cbind(Populations$X, Populations$Y)

a<-area(patch)

Populations$area<-extract(a, pxy)

id<-pyrs$ID
lam<-rep(1,length(id))    #not sure what the value here pertains to - think it sets starting population so should use values from LPI?
pyrxy<-SpatialPoints(pyr[,c("Longitude","Latitude")])
sdm<-raster(paste(wd, "/hyde_pres_abs_sss_weighted_ensemble_sdm_1950.tif", sep=""))
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
df<-data.frame(na.omit(df))

df<-rbind(df, Populations)

dfxy<-cbind(df$X, df$Y)
df$area<-extract(a, dfxy)



# buff<-0.75
# 
# e<-extent(df_test$X - buff, df_test$X + buff, df_test$Y - buff, df_test$Y + buff)
# 
# df<-df[df$X > e@xmin & df$X < e@xmax & df$Y < e@ymax & df$Y > e@ymin,]
#  

###formatting environmental data
patch<-raster(paste(wd,"/hyde_pres_abs_sss_weighted_ensemble_sdm_1950.tif", sep=""))
#patch<-crop(patch, e)

sdm_patch_df<-data.frame(ID=1:ncell(patch))
sdm_df<-data.frame(ID=1:ncell(patch))

for (i in 1:length(years)){
  sdm<-raster(paste(wd,"/hyde_weighted_ensemble_sdm_", years[i],".tif", sep=""))   #14.6
  patch<-raster(paste(wd,"/hyde_pres_abs_sss_weighted_ensemble_sdm_", years[i],".tif", sep=""))
  #sdm<-crop(sdm, e)
  #patch<-crop(sdm, e)
  
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
opt_patch_spin_up[,4:length(spin_years)+4]<-1

patch_map_mine[is.na(patch_map_mine)]<-0

niche_spin_up[is.na(niche_spin_up)]<-0
patch_spin_up[is.na(patch_spin_up)]<-0

opt_patch_spin_up

###COMADRE
load(paste(wd, "COMADRE_v.2.0.1.RData", sep="/"))

tempMetadata<-subset(comadre$metadata, SpeciesAccepted==binomial)

keep<-as.numeric(rownames(tempMetadata))

tempMat<-comadre$mat[keep]   #MatA is matrix pop model, can be split into U, F and/or C

MatList<-list(tempMat[[1]][[1]])  #varies depending on number of matrices - need to find a way to code this better - now have five matrices available so need to sort this
AllMat<-unlist(MatList)
matrices<-matrix(AllMat, ncol=length(MatList))
colnames(matrices)<- c("Reference_matrix")

stages<-comadre$matrixClass[keep][[1]]$MatrixClassAuthor
proportion_initial<- rep(1/length(stages), length(stages)) #Think spin up sorts this
sumweight<-rep(1, length(stages))#weight of stages  - should be equal for all mine just
#in plants seed not included in calculating population sizes
list_names_matrices<-colnames(matrices)
K_weight<-c(rep(1, length(stages)))  #the weight with which carrying capacity affects each stage was FALSE


####Scaling Carrying Capacity

lin<-function(x, carry_k){
  x*carry_k
}

lf<-list.files(wd)

files<-lf[grepl("^hyde_weighted_ensemble_sdm_.*.tif$", lf)]

sdms<-stack(paste(wd,files, sep="/"))

hsi<-raster:::extract(sdms, Populations[,c(2,3)])

link1<-lin(hsi, carry_k[1])
link1<-link1*Populations$area

spin1<-replicate(length(spin_years),link1[,1])
link_spin1<-cbind(spin1, link1)
colnames(link_spin1)[1:length(spin_years)]<-paste("hyde_weighted_ensemble_sdm_", spin_years, sep="")

link_spin<-list(link_spin1)

s=1

SDD<-var_grid[s,"SDD"]
LDD<-var_grid[s,"LDD"]
kern<-var_grid[s, "Kern"]
SD<-var_grid[s, "SD"]

K<-var_grid[s, "K"]
dens<-var_grid[s, "Density"]
link_id<-var_grid[s, "K_scale"]

matrices_var<-matrix(SD, ncol = 1, nrow = nrow(matrices), dimnames = list(NULL, "sd"))

start.time <- Sys.time()

med_disp<-as.character(kern_seq[[kern]][1])

link_k<-matrix(unlist(link_spin[link_id]), nrow=nrow(link_spin1), ncol=ncol(link_spin1))
colnames(link_k)<-colnames(link_spin1)

lower_est = c(0.5,0.5,0.5,0.08,0.5,0.08,0.5,0.08,0.5,0.08,0.5,0.08,0.5)
upper_est = c(1,1,1,1.5,1,1.5,1,1.5,1,1.5,1,1.5,1)

csv_read<-function(file){
  matf<-read.csv(file)
  return(matf)
}


file.remove(list.files(pattern = "*matrix_spin_up"))

cal_demoniche=function(x) {
  matrices[c(2,11,20,25,29,33,38,41,47,49,56,57,64)]=x
  print(x)
  
  demoniche_setup_me(modelname = binomial ,Populations = df, Nichemap = opt_patch_spin_up,
                     matrices = matrices,matrices_var = matrices_var, prob_scenario = prob_scenario,
                     stages = stages, proportion_initial = proportion_initial,
                     density_individuals = dens,
                     fraction_LDD = LDD, fraction_SDD = SDD,
                     dispersal_constants = kern_seq[[kern]],
                     transition_affected_niche = transition_affected_niche,
                     transition_affected_demogr = F,
                     transition_affected_env=F,
                     env_stochas_type = env_stochas_type,
                     no_yrs = no_yrs_mine, K=5000, Kweight = K_weight, Ktype="ceiling",
                     sumweight =sumweight)
  
  demoVE_model=demoniche_model_me(binomial,Niche=T,Dispersal=T,repetitions=1,foldername=binomial)	

  
  files<-paste("matrix_spin_up_",89:99, ".csv", sep="")
  file_out<-lapply(files, csv_read)
  mats<do.call("cbind", files_out)
  
  sum(abs(rowMeans(mats)[which(rowMeans(mats)>0)]-lower_est))
  
  #sum(abs(rowMeans(spin_mat)[which(rowMeans(spin_mat)>0)]-upper_est))
  #file.copy("matrix_spin_up.csv", "matrix_spin_up_copy.csv")
  #file.remove("matrix_spin_up.csv")
  
}

cal_mat=optim(c(0.93,0.93,0.93,0.28,0.93,0.28,0.93,0.28,0.93,0.28,0.93,0.28,0.93),cal_demoniche,lower=lower_est,upper=upper_est,method='L-BFGS-B', control = list(fnscale = -2)) # box constraint


save(cal_mat, "matrix_optim.Rda")








