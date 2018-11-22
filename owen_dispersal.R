

library(raster)

#each cell, each age class


dispersal_probabilities<-read.table("D:/Fiona/Git_Method/Git_Method/Legion/snow_bear_bias_faster/disp_probs_hc_max.csv")
patch<-raster("D:/Fiona/Git_Method/Git_Method/Legion/snow_bear_bias/hyde_pres_abs_sss_weighted_ensemble_sdm_1950.tif")
values(patch)<-as.numeric(dispersal_probabilities[1000,])

cumsum_disp <- function(x) {
  test <-rbind(1:ncol(dispersal_probabilities), dispersal_probabilities[x,])
  empty <- which(test[2,] == 0) #identifying which cells are sea
  if (length(empty) > 0) {
    test_em <- test[,-empty]
    test_em <- as.matrix(test_em, nrow = 2, ncol = 1)
    test_em[2,] <- cumsum(as.numeric(test_em[2,]))
    test_em[1,] <- as.integer(test_em[1,])
    
  } else {
    test_em <- matrix(nrow = 2, ncol = 1)
    test_em[1, 1] <- 0
    test_em[2, 1] <- 0
  }
  print(x)
  return(test_em)
  
}

list_disp<-lapply(1:nrow(dispersal_probabilities), cumsum_disp)

save(list_disp, file = "cumulative_dispersal_probabilities.Rdata")


load("D:/Fiona/Git_Method/Git_Method/Legion/snow_bear_bias_faster/cumulative_dispersal_probabilities.Rdata")

#same as source_patches


stages<-c( "1 yr","2 yrs","3 yrs","4 yrs","5 yrs","6 yrs","7 yrs","8 yrs","9 yrs","10 yrs","11 yrs","12 yrs","13 yrs","14 yrs","15 yrs","16 yrs","17 yrs","18 yrs","19 yrs")
seeds_per_population<-matrix(100,ncol = length(list_disp), nrow = 19) #setting to 100 individuals in each stage in each cell 
# sea<-as.numeric(which(colSums(dispersal_probabilities) == 0)) #getting rid of any individuals in the sea
# seeds_per_population[,sea]<-0

seeds_per_population_new_LDD <- matrix(0, nrow = nrow(seeds_per_population), ncol = ncol(seeds_per_population)) #each row is an age class, each column is a cell

fraction_LDD<-0.3 #dispersal rate from the cell - this varies with each model run

seeds_per_population_migrate_LDD <- round(seeds_per_population * fraction_LDD) #number of individuals to disperse each year from each cell and age class

source_patches_ldd<-which(colSums(seeds_per_population_migrate_LDD)>0)  #which cells have individuals in - all of them do initially.

which_ap<-function(probs, cumsm){
  
  if(length(cumsm) != 0 &length(cumsm)!=0 ){
    out<-min(which(cumsm[2,] >= probs)) #the cell with cumulative prob just over the randomly selected probability
    cell<-cumsm[1,out]
    return(cell)
  }
}



cell_ind<-function(ind, stage){
  print(ind)
  num<-seeds_per_population_migrate_LDD[stage,ind]
  cumsm<-list_disp[[ind]]
  probs<-runif(num,0,1)
  if(length(probs)>0){
    mig_pops<-data.frame(table(sapply(probs, which_ap, cumsm = cumsm)))
    mig_pops$Var1<-as.numeric(as.character(mig_pops$Var1))
    mig_pops$Freq<-as.numeric(as.character(mig_pops$Freq))
  } else{
    mig_pops<-data.frame(Var1 = 1,Freq = 0)
  }
  seeds_per_population_new_LDD[stage,mig_pops$Var1]<-seeds_per_population_new_LDD[stage,mig_pops$Var1]+mig_pops$Freq
  return(seeds_per_population_new_LDD[stage,])
}


# seeds_per_population_new_LDD<-matrix(0,nrow=nrow(seeds_per_population), ncol = ncol(seeds_per_population))
for (i in 1:nrow(seeds_per_population)) {
  seeds_per_population_new_LDD[i, ] <- rowSums(sapply(source_patches_ldd, cell_ind, stage =  i))
  print(i)
}

values(patch)<-seeds_per_population_new_LDD[1,]
plot(patch)



