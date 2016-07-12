body<-read.csv("bird_and_mammal_traits2.csv")
luc<-read.csv("LUC_average_annual_change_nat_025.csv")
LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")
Realm<-read.csv("selected_pops_Ecoregion.csv")
Realm<-Realm[,c("ID", "WWF_REALM2")]

pop<-read.csv("Global_Population_Trends_Rsq_Lambda_16_03_18.csv")
#EurHil<-read.csv("Europe_HILDA_5_year_pops.csv")  # data from Euro-centric analysis

temp<-temp[,c("ID", "Estimate")]

LPI<-LPI[,c("ID","Binomial","Common_name", "Order", "Protected_status", "Country","Region", "System", "Class","Specific_location", "Longitude", "Latitude", "Primary_threat", "Secondary_threat", "Tertiary_threat", "Migratory")]

df<-merge(merge(temp,luc, by="ID", all=TRUE), merge(LPI, pop, by="ID", all=TRUE),by="ID", all=TRUE)

df<-merge(df, Realm, by="ID")

df<-merge(df, body[,c(3:5)], by="ID", all=TRUE)     #41 pops bodysizes missing for birds

df2<-subset(df, !is.na(Estimate)&r_sq >= 0.5  & !is.na(Nat_change)&length_time >=10 & System!="Marine" &Specific_location == 1 & (Class=="Mammalia") & !is.na(Bodymass) )

nrow(df2)

df2$Nat_loss<-0-df2$Nat_change 

df2[is.na(df2$lambda_mean),]$lambda_mean<-0

pyr<-subset(df2, Common_name =="Pyrenean chamois")    #record 11470 had wrong longitude - in Russia!

plot(pyr$Longitude, pyr$Latitude)

pyrs<-pyr[,c("ID","Longitude","Latitude","lambda_mean")]

id<-pyrs$ID
#lam<-as.numeric(pyrs$lambda_mean)
lam<-rep(200,length(id))    #not sure what the value here pertains to

library(rgdal)

pyrxy<-SpatialPoints(pyr[,c("Longitude","Latitude")])

library(raster)

e<-extent(pyrxy)
xmn<-floor(e[1])
xmx<-ceiling(e[2])
ymn<-floor(e[3])
ymx<-ceiling(e[4])

e2<-extent(xmn,xmx,ymn,ymx)

r<-raster(e2, resolution=0.25)

rz<-rasterize(pyrxy,r,lam )
rid<-rasterize(pyrxy,r,id)
plot(rz)
plot(rid)

rz_spdf<-xyFromCell(rz, 1:ncell(rid))

rzm<-as.vector(rz)
ridm<-as.vector(rid)

df<-data.frame(ridm,rz_spdf,rzm)
colnames(df)<-c( "PatchID","X","Y","area")

Orig_Populations<-na.omit(df)   #df without squares with missing pops


write.csv(Orig_Populations, "Rupicapra_pyrenaica_populations.csv")

#############geographic data for populations#######
# lnd<-read.csv("LUH2_Land_Use_All_LPI.csv")
# 
# lnd_rp<-merge(df,lnd, by="ID")
# 
# lnd_rp2<-unique(lnd_rp[,c("ID","patchID", "x", "y")])
# 
# library(reshape2)
# 
# lnd_rp2<-lnd_rp[,c("ID","Year","natural")]
# 
# lnd_rp_cast<-dcast(lnd_rp2, ID ~ Year )
# 
######geographic data in grids###############

library(raster)
library(ncdf4)

yr<-as.character(850:2015)
date<-as.Date(yr, format="%Y" )

primf<-brick("states.nc", varname="primf")
primf<-setZ(primf, date, name="year")

primn<-brick("states.nc", varname="primn")
primn<-setZ(primn, date, name="year")

secdf<-brick("states.nc", varname="secdf")
secdf<-setZ(secdf, date, name="year")

secdn<-brick("states.nc", varname="secdn")
secdn<-setZ(secdn, date, name="year")

primf_cr<-primf[[1091:1166]]
primn_cr<-primn[[1091:1166]]
secdf_cr<-secdf[[1091:1166]]
secdn_cr<-secdn[[1091:1166]]


# writeRaster(primf_cr, "primf_1940.tif")
# writeRaster(primn_cr, "primn_1940.tif")
# writeRaster(secdf_cr, "secdf_1940.tif")
# writeRaster(secdn_cr, "secdn_1940.tif")

primf<-brick("primf_1940.tif")
primn<-brick("primn_1940.tif")
secdf<-brick("secdf_1940.tif")
secdn<-brick("secdn_1940.tif")

n40<-primf[[1]]+ primn[[1]]+secdf[[1]]+secdn[[1]]
n40c<-crop(n40, r)

n50<-primf[[11]]+ primn[[11]]+secdf[[11]]+secdn[[11]]
n50c<-crop(n50,r)

n60<-primf[[21]]+ primn[[21]]+secdf[[21]]+secdn[[21]]
n60c<-crop(n60,r)

n70<-primf[[31]]+ primn[[31]]+secdf[[31]]+secdn[[31]]
n70c<-crop(n70,r)

n80<-primf[[41]]+ primn[[41]]+secdf[[41]]+secdn[[41]]
n80c<-crop(n80,r)

n90<-primf[[51]]+ primn[[51]]+secdf[[51]]+secdn[[51]]
n90c<-crop(n90,r)

n2000<-primf[[61]]+ primn[[61]]+secdf[[61]]+secdn[[61]]
n2000c<-crop(n2000,r)

n2010<-primf[[71]]+ primn[[71]]+secdf[[71]]+secdn[[71]]
n2010c<-crop(n2010,r)


n40m<-as.vector(n40c)
n50m<-as.vector(n50c)
n60m<-as.vector(n60c)
n70m<-as.vector(n70c)
n80m<-as.vector(n80c)
n90m<-as.vector(n90c)
n2000m<-as.vector(n2000c)
n2010m<-as.vector(n2010c)

Niche_values<-cbind(n40m,n50m,n60m,n70m,n80m,n90m,n2000m,n2010m)
colnames(Niche_values)<-c("period_1940", "period_1950", "period_1960", "period_1970", "period_1980", "period_1990", "period_2000", "period_2010")

Niche_ID<-1:length(n40m) + 5000
land_use_map<-cbind(Niche_ID,df[,c(2,3,1)]) #percentage cover of natural land use (primary+secondary) for each decade 1940-2010

land_use_map[is.na(land_use_map)]<-0
colnames(land_use_map)[3]<-"y"

years_projections<-colnames(Niche_values)

no_yrs<-10
######## heat map of niche values (natural land use cover 1940-2010)

nichemap<-cbind(land_use_map[,c(1:3)],Niche_values )

niche_formulas <- as.formula(paste(paste(colnames(nichemap)[-c(1:3)],collapse="+"),"x+y",sep="~"))

print(levelplot(niche_formulas, nichemap, col.regions=rev(heat.colors(100)), main = "Niche Values"))


###### demographic information
#install.packages("demoniche", repos="http://R-Forge.R-project.org")

library(popbio)
library(demoniche)

dir<-getwd()
load(paste(dir, "COMADRE_v.1.0.0.RData", sep="/"))

tempMetadata<-subset(comadre$metadata, SpeciesAccepted=="Ovis_aries")

keep<-as.numeric(rownames(tempMetadata))

tempMat<-comadre$mat[keep]   #MatA is matrix pop model, can be split into U, F and/or C

#first matrix in list is the mean of others?

#tempMatmean<-(tempMat[[1]][[1]]+ tempMat[[2]][[1]]+tempMat[[3]][[1]]+tempMat[[4]][[1]]+tempMat[[5]][[1]]+tempMat[[6]][[1]]+tempMat[[7]][[1]])/length(keep)

MatList<-list(tempMat[[1]][[1]], tempMat[[2]][[1]], tempMat[[3]][[1]],tempMat[[4]][[1]],tempMat[[5]][[1]],tempMat[[6]][[1]],tempMat[[7]][[1]])
AllMat<-unlist(MatList)
matrices<-matrix(AllMat, ncol=7)
colnames(matrices)<- c("Reference_matrix", "Matrix_1", "Matrix_2", "Matrix_3", "Matrix_4", "Matrix_5", "Matrix_6")


prob_scenario<-c(0.5,0.5)    #need to check this

noise<-0.95     #need to check this

stages<-comadre$matrixClass[keep][[1]]$MatrixClassAuthor

list_names_matrices<-colnames(matrices)

sumweight<-c(1,1,1,1,1,1)  #weight of stages  - should be equal for all mine just in plants seed not included in calculating population sizes - or if you wanted to just calculate the female population it would be c(1,1,1,0,0,0)


transition_affected_niche<-"all"    #which parts of the matrix are affected by the niche values

transition_affected_env <- "all"

transition_affected_demogr <- "all"

env_stochas_type<-"normal"   #can also be lognormal

matrices_var <- matrix(0.01, ncol = 1, nrow = nrow(matrices), dimnames = list(NULL, "sd")) #standard deviation of matrices

proportion_initial<- c(0.9818098089, 0.0006907668, 0.0069076675,0.0036840893, 0.0057563896, 0.0011512779)  #proportion of population in each stage - no idea what this should be and will likely have a big impact on results!

density_individuals <- 20000   #also compulsory not sure what best value would be

K<-NULL

K_weight<-FALSE

fraction_SDD <- 0.05  #short distance dispersal

#fraction_LDD_mine <- 0.05 #long distance dispersal


###Minimal Setup #Remember that the demographic info is for different species - ovis aries
##doesn't work with ovis aries because the matrix is split into male and female so used the default from popbio

library(popbio)
data(hudvrs)
data(hudsonia)
matrices_mine <- cbind(meanmatrix = as.vector(hudmxdef(hudvrs$mean)),sapply(hudsonia, unlist))
head(matrices_mine)
colnames(matrices_mine) <-c("Reference_matrix", "Matrix_1", "Matrix_2", "Matrix_3", "Matrix_4")

hudsonia_stages<-c("seeds", "seedlings", "tiny", "small", "medium", "large")
####################################

demoniche_setup(modelname = "RPyran_minimal",Populations = Orig_Populations, matrices_var = matrices_var,matrices = matrices_mine,
                stages = hudsonia_stages, proportion_initial = proportion_initial,density_individuals = density_individuals,
                no_yrs = 10)

RPyran_min_run <- demoniche_model(modelname = "RPyran_minimal", Niche = FALSE, Dispersal = FALSE, repetitions = 1,foldername = "RPyran_minimal")
dimnames(RPyran_min_run)

RPyran_min_run[,"Meanpop","Reference_matrix"]



(hudsonia)














