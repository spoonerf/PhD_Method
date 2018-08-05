#seeds_per_population is Projection[yx, 1, , tx], which is from the population function - need to find out if it is multiplied by niche values at any point

# demoniche_dispersal_me<-function (seeds_per_population, fraction_LDD, fraction_SDD, dispersal_probabilities, 
#           dist_latlong, neigh_index) 
demoniche_dispersal_csv<-function (seeds_per_population, fraction_LDD, fraction_SDD, dispersal_probabilities, 
                                    dist_latlong, neigh_index, niche_values) 
{
  seeds_per_population_migrate_LDD <- round(seeds_per_population *    #number of seeds to ldd - can this number be larger than 1? added round
    fraction_LDD)
  seeds_per_population_migrate_SDD <- round(seeds_per_population * 
    fraction_SDD)
  #  seeds_per_population_new_SDD <- seeds_per_population_new_LDD <- rep(0, 
  #                                                                      length(seeds_per_population))
  seeds_per_population_new_SDD <- seeds_per_population_new_LDD <- matrix(0, nrow = nrow(seeds_per_population), ncol = ncol(seeds_per_population))
  
  if (fraction_SDD > 0) {
    # source_patches <- which(seeds_per_population_migrate_SDD >
    #                           0)
    source_patches <- which(colSums(seeds_per_population_migrate_SDD) >
                              0)
    for (px_orig in source_patches) {
      for (pxdisp_new in 1:length(seeds_per_population_migrate_SDD[1,])) {
        if (dist_latlong[pxdisp_new, px_orig] == neigh_index[1]) {
          # seeds_per_population_new_SDD[pxdisp_new] <- seeds_per_population_new_SDD[pxdisp_new] + 
          #   (seeds_per_population_migrate_SDD[px_orig] * 
          #      0.2)
          seeds_per_population_new_SDD[,pxdisp_new] <- round(seeds_per_population_new_SDD[,pxdisp_new] + 
            (seeds_per_population_migrate_SDD[,px_orig] * 
               0.2))
        }
        if (length(neigh_index) == 2) {
          if (dist_latlong[pxdisp_new, px_orig] == neigh_index[2]) {
            # seeds_per_population_new_SDD[pxdisp_new] <- seeds_per_population_new_SDD[pxdisp_new] +
            #   (seeds_per_population_migrate_SDD[px_orig] *
            #      0.05)
            seeds_per_population_new_SDD[,pxdisp_new] <- round(seeds_per_population_new_SDD[,pxdisp_new] +
              (seeds_per_population_migrate_SDD[,px_orig] *
                 0.05))
          }
        }
      }
    }
  }

  # install.packages('ncdf4', repos="http://cran.r-project.org")
  # library(ncdf4)

 
   if (fraction_LDD > 0) {
    
    # for (i in 1:nrow(seeds_per_population_migrate_LDD)){
    #   
    #  seeds_per_population_new_LDD[i,] <- round(dispersal_probabilities %*%  seeds_per_population_migrate_LDD[i,])
    #   #seeds_per_population_new_LDD[i,] <- disp_prob %*%  seeds_per_population_migrate_LDD[i,]
    #  
    # }
    
    #  source_patches_ldd<-which(colSums(seeds_per_population_migrate_LDD)>0)
    #  dp<-dispersal_probabilities$var[[1]]
    #  varsize <- dp$varsize
    #  ndims   <- dp$ndims  
    #  nt      <- varsize[1] #to go row by row?
  #install.packages("RNetCDF")
  # library(RNetCDF)
  # 
  # source_patches_ldd<-which(colSums(seeds_per_population_migrate_LDD)>0)
  # dp<-file.inq.nc(dispersal_probabilities)
  # varsize <- c(length(var.get.nc(dispersal_probabilities, 0)), length(var.get.nc(dispersal_probabilities, 0)))
  # ndims   <- dp$ndims
  # nt      <- varsize[1] #to go row by row?

  
  #library(RNetCDF)
  source_patches_ldd<-which(colSums(seeds_per_population_migrate_LDD)>0)
  

  for(i in 1:nrow(seeds_per_population_new_LDD)){ 
    for (j in source_patches_ldd){
    
      #if (seeds_per_population_migrate_LDD[i,j] > (1/fraction_LDD)){
      
      disp_prob<-read.table("disp_probs.csv", nrow = 1, skip = j - 1, sep=" ")
      disp_prob[is.na(disp_prob)]<-0
      
      if (disp_prob[1] >= 1){
        disp_prob<-disp_prob[-1]
        }
      
      if (sum(disp_prob)> 0){
      new_patches <-sample(1:length(disp_prob),seeds_per_population_migrate_LDD[i,j],prob=disp_prob, replace =TRUE)
      seeds_per_population_new_LDD[i,new_patches]<-seeds_per_population_new_LDD[i,new_patches] +1
      np_df<-as.data.frame(table(new_patches))
      np_df$new_patches<-as.numeric(as.character(np_df$new_patches))
      np_df$Freq<-as.numeric(as.character(np_df$Freq))
      } else {
        np_df<-data.frame(new_patches = 1, Freq = 0)
      }
      print(j)
      
      for(k in 1:nrow(np_df)){
      seeds_per_population_new_LDD[i,np_df$new_patches[k]]<-seeds_per_population_new_LDD[i,np_df$new_patches[k]] +np_df$Freq[k]
        }
      }
     #seeds_per_population_new_LDD <- as.vector(dispersal_probabilities %*% seeds_per_population_migrate_LDD[1,])
  }
  
   
  } 
  seeds_stay <- (seeds_per_population - seeds_per_population_migrate_SDD -   #seeds that migrate must go out of the cell and are taken off the total
                   seeds_per_population_migrate_LDD)
  
    #seeds_stay<-niche_values*seeds_stay
  #seeds_total<-niche_values*(seeds_stay + seeds_per_population_new_SDD + seeds_per_population_new_LDD)
  return(seeds_stay + seeds_per_population_new_SDD + seeds_per_population_new_LDD)
  #return(seeds_stay)

}
