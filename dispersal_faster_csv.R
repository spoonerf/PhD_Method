
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

    source_patches_ldd<-which(colSums(seeds_per_population_migrate_LDD)>0)
    
    #disp_prob<-read.table(dispersal_probabilities, sep=" ")
    disp_prob<-dispersal_probabilities
    #disp_prob[is.na(disp_prob)]<-0    
    
    sample_ages<-function(x){
      new_patches<-sample(1:length(dp),x,prob=dp, replace =TRUE)
      return(new_patches)
    }
    
    for (j in source_patches_ldd){
      
      dp<-disp_prob[j,]
      
      if (dp[1] >= 1){
        dp<-dp[-1]
      }
      
      if (sum(dp)> 0){
        
        ages<-matrix(seeds_per_population_migrate_LDD[,j], ncol = 1)
        
        stage_out<-apply(ages,1,sample_ages)
        
        # disp_stages<-function(x){
        #   a<-as.data.frame(table(stage_out))
        # }
        # 
        for(k in 1:length(stage_out)){
          a<-data.frame(table(stage_out[[k]]))
          a$Var1<-as.numeric(as.character(a$Var1))
          a$Freq<-as.numeric(as.character(a$Freq))
          
          if (nrow(a)>0){
            seeds_per_population_new_LDD[k,a$Var1]<-seeds_per_population_new_LDD[k,a$Var1] +a$Freq
            #print(seeds_per_population_new_LDD[k,a$Var1])
          } else{
            seeds_per_population_new_LDD[k,]<-seeds_per_population_new_LDD[k,]
          }
        }
      }
      #seeds_per_population_new_LDD <- as.vector(dispersal_probabilities %*% seeds_per_population_migrate_LDD[1,])
    #print(paste("cell",j,sep=" "))
      }
  }
  
  
  seeds_stay <- (seeds_per_population - seeds_per_population_migrate_SDD -   #seeds that migrate must go out of the cell and are taken off the total
                   seeds_per_population_migrate_LDD)
  print(sum(seeds_stay))
  #seeds_stay<-niche_values*seeds_stay
  #seeds_total<-niche_values*(seeds_stay + seeds_per_population_new_SDD + seeds_per_population_new_LDD)
  return(seeds_stay + seeds_per_population_new_SDD + seeds_per_population_new_LDD)
  #return(seeds_stay)
  
}
