
# transition_affected_niche <- c(2,29)
# transition_affected_env <- c(2,29)

demoniche_population_me<-function (Matrix_projection, Matrix_projection_var, n, populationmax, 
                                   K = NULL, Kweight = BEMDEM$Kweight, onepopulation_Niche, 
                                   sumweight, noise, prob_scenario, prev_mx, transition_affected_demogr, 
                                   transition_affected_niche, transition_affected_env, env_stochas_type, 
                                   yx_tx, tx, yrs_total) 
{
  prob_scenario_noise <- c(prob_scenario[prev_mx[yx_tx]] * 
                             noise, 1 - (prob_scenario[prev_mx[yx_tx]] * noise))   #this all seems strange and samples from a first or second matrix
  rand_mxs <- sample(1:2, 1, prob = prob_scenario_noise, replace = TRUE)
  
  one_mxs <- Matrix_projection[, rand_mxs]
  prev_mx[yx_tx + 1] <- rand_mxs
  if (Matrix_projection_var[1] != FALSE) {
    one_mxs_var <- one_mxs * (Matrix_projection_var[, rand_mxs])
    if (is.numeric(transition_affected_niche)) {
      one_mxs[transition_affected_niche] <- one_mxs[transition_affected_niche] * 
        onepopulation_Niche
    }
    if (is.numeric(transition_affected_env)) {
      #normal or log-normal effects of env      #sd here has already been multiplied by the corresponding values in the matrix (here 0.93)
      switch(EXPR = env_stochas_type, normal = one_mxs[transition_affected_env] <- rnorm(length(one_mxs[transition_affected_env]), 
                                                                                         mean = one_mxs[transition_affected_env], sd = one_mxs_var[transition_affected_env]), 
             lognormal = one_mxs[transition_affected_env] <- rlnorm(length(one_mxs[transition_affected_env]), 
                                                                    meanlog = one_mxs[transition_affected_env], 
                                                                    sdlog = one_mxs_var[transition_affected_env]))
    }
  }
  
  one_mxs[one_mxs < 0] <- 0   #changing any negative values to zero
  A <- matrix(one_mxs, ncol = length(n), nrow = length(n), 
              byrow = FALSE)
  Atest <- A
  Atest[1, ][-1] <- 0   #getting fertility values 
  if (sum(colSums(Atest) > 1)) {      #picking out survival scores higher than 1 
    for (zerox in which(colSums(Atest) > 1)) {
      Atest[, zerox] <- Atest[, zerox]/sum(Atest[, zerox])
    }
    A[-1, ] <- Atest[-1, ] #changing survival values
  }
  
  #yrs_total added to get matrices for last 10 years of spin up - writes for each cell which has a pop?
  
  #file.remove("matrix_spin_up.csv")
    Am<-matrix(A, ncol=1)
    if (tx >= 89 & tx <=99 ){
    if (file.exists(paste("matrix_spin_up_",tx,".csv", sep=""))){
      tmp<-read.csv(paste("matrix_spin_up_",tx,".csv", sep=""))
      new<-cbind(tmp, Am)
      write.csv(new, paste("matrix_spin_up_",tx,".csv", sep=""), row.names = FALSE)
    } else{
      write.csv(Am, paste("matrix_spin_up_",tx,".csv", sep=""), row.names=FALSE)      
    }
    }
  ##me
  n[is.na(n)]<-0
  ##
  
  
  n <- as.vector(A %*% n)    #n is the number of ibex in each stage of the matrix - a row from n0s which is all of the populations - here it is multipled by the matrix
  n <- floor(n)
  
  if (sum(n) > 0) {
    if (is.numeric(populationmax)) {
      if (sum(n * Kweight) > populationmax) {
        n <- n * (populationmax/sum(n * sumweight))   #where carrying capacity comes in - brings n back down to carrying capacity 
      }
    }
  }
  print(n[1])
  return(n)
}

