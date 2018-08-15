
#seeds_per_population is Projection[yx, 1, , tx], which is from the population function - need to find out if it is multiplied by niche values at any point

# demoniche_dispersal_me<-function (seeds_per_population, fraction_LDD, fraction_SDD, dispersal_probabilities, 
#           dist_latlong, neigh_index) 
demoniche_dispersal_csv<-function (seeds_per_population, fraction_LDD, fraction_SDD, dispersal_probabilities, 
                                   dist_latlong, neigh_index, niche_values, stages) 
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
  
  if (fraction_LDD > 0) {

    source_patches_ldd<-which(colSums(seeds_per_population_migrate_LDD)>0)
    
    disp_prob<-dispersal_probabilities
    
    sample_ages<-function(x, dp){
      new_patches<-sample(1:length(dp),x,prob=dp, replace =TRUE)
      return(new_patches)
    }
    
    disp_out<-function(stg,seeds_new, stage_out){
      
      a<-data.frame(table(stage_out[[stg]])) #slowness
      a$Var1<-as.numeric(as.character(a$Var1))
      a$Freq<-as.numeric(as.character(a$Freq))
      
      if (nrow(a)>0){
        seeds_new[stg,a$Var1]<-seeds_new[stg,a$Var1] +a$Freq
      } else {
        seeds_new[stg,]<-seeds_new[stg,]
      }
      return(seeds_new[stg,])
    }
    
   
    dispersal<-function(j, seeds_ldd, disp_prob){
      
      dp<-disp_prob[j,]
      print(j)
      ages<-matrix(seeds_ldd[,j], ncol = 1)
      stage_out<-lapply(X = ages, FUN = sample_ages, dp = dp )
      st<-as.matrix(1:length(stage_out))
      
      if(length(stage_out)>0){
        seeds_new_out<-sapply( X = st,FUN  =  disp_out,  seeds_new =  seeds_new, stage_out = stage_out)
        #seeds_newt<-t(seeds_new_out)
        
        return(seeds_new_out)
        
      } else {
        
        return(seeds_new)  
        
      }
    }
    
    seeds_new<-seeds_per_population_new_LDD
    
    source_patches<-as.matrix(source_patches_ldd, ncol = 1)
    
    check_disp<-apply(X = source_patches,1, FUN = dispersal, seeds_ldd = seeds_per_population_migrate_LDD, disp_prob = disp_prob)
    
    seeds_per_population_new_LDD<-matrix(rowSums(check_disp), ncol = ncol(seeds_per_population_new_LDD) , nrow = length(stages), byrow = T)
    
    print(paste("no. source_patches ", length(source_patches_ldd), sep=""))
    print(paste("no. to migrate ",sum(colSums(seeds_per_population_migrate_LDD)), sep=""))
    print(paste("no. which migrated ", sum(colSums(seeds_per_population_new_LDD)), sep=""))
    
    }
      

  seeds_stay <- (seeds_per_population - seeds_per_population_migrate_SDD -   #seeds that migrate must go out of the cell and are taken off the total
                   seeds_per_population_migrate_LDD)
  print(sum(seeds_stay+ seeds_per_population_new_SDD + seeds_per_population_new_LDD))
  
  return(seeds_stay + seeds_per_population_new_SDD + seeds_per_population_new_LDD)
}
