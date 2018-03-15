#seeds_per_population is Projection[yx, 1, , tx], which is from the population function - need to find out if it is multiplied by niche values at any point

# demoniche_dispersal_me<-function (seeds_per_population, fraction_LDD, fraction_SDD, dispersal_probabilities, 
#           dist_latlong, neigh_index) 
demoniche_dispersal_me<-function (seeds_per_population, fraction_LDD, fraction_SDD, dispersal_probabilities, 
                                    dist_latlong, neigh_index, niche_values) 
{
  seeds_per_population_migrate_LDD <- seeds_per_population *    #number of seeds to ldd - can this number be larger than 1?
    fraction_LDD
  seeds_per_population_migrate_SDD <- seeds_per_population * 
    fraction_SDD
  #  seeds_per_population_new_SDD <- seeds_per_population_new_LDD <- rep(0, 
  #                                                                      length(seeds_per_population))
  seeds_per_population_new_SDD <- seeds_per_population_new_LDD <- matrix(0, nrow = nrow(seeds_per_population), ncol = ncol(seeds_per_population))
  if (fraction_SDD > 0) {
    # source_patches <- which(seeds_per_population_migrate_SDD >
    #                           0)
    source_patches <- which(colSums(seeds_per_population_migrate_SDD) >
                              0)
    for (px_orig in source_patches) {
      for (pxdisp_new in 1:length(seeds_per_population_migrate_SDD)) {
        if (dist_latlong[pxdisp_new, px_orig] == neigh_index[1]) {
          # seeds_per_population_new_SDD[pxdisp_new] <- seeds_per_population_new_SDD[pxdisp_new] + 
          #   (seeds_per_population_migrate_SDD[px_orig] * 
          #      0.2)
          seeds_per_population_new_SDD[,pxdisp_new] <- seeds_per_population_new_SDD[,pxdisp_new] + 
            (seeds_per_population_migrate_SDD[,px_orig] * 
               0.2)
        }
        if (length(neigh_index) == 2) {
          if (dist_latlong[pxdisp_new, px_orig] == neigh_index[2]) {
            # seeds_per_population_new_SDD[pxdisp_new] <- seeds_per_population_new_SDD[pxdisp_new] +
            #   (seeds_per_population_migrate_SDD[px_orig] *
            #      0.05)
            seeds_per_population_new_SDD[,pxdisp_new] <- seeds_per_population_new_SDD[,pxdisp_new] +
              (seeds_per_population_migrate_SDD[,px_orig] *
                 0.05)
          }
        }
      }
    }
  }
  if (fraction_LDD > 0) {
    
    for (i in 1:nrow(seeds_per_population_migrate_LDD)){
      
      seeds_per_population_new_LDD[i,] <- dispersal_probabilities %*%  seeds_per_population_migrate_LDD[i,]
      
    }
    
    # seeds_per_population_new_LDD <- as.vector(dispersal_probabilities %*% 
    #                                             seeds_per_population_migrate_LDD[1,])
  }
  seeds_stay <- (seeds_per_population - seeds_per_population_migrate_SDD -   #seeds that migrate must go out of the cell and are taken off the total
                   seeds_per_population_migrate_LDD)
  #seeds_stay<-niche_values*seeds_stay
  #seeds_total<-niche_values*(seeds_stay + seeds_per_population_new_SDD + seeds_per_population_new_LDD)
  return(seeds_stay + seeds_per_population_new_SDD + seeds_per_population_new_LDD)
  #return(seeds_stay)
>>>>>>> 6350399a964b1683bd7f56a1c7b22069732e8720
}
