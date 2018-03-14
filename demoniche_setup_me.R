demoniche_setup_me<-function (modelname, Populations, stages, Nichemap = "oneperiod", 
          matrices, matrices_var = FALSE, prob_scenario = c(0.5, 0.5), 
          proportion_initial, density_individuals, transition_affected_niche = FALSE, 
          transition_affected_env = FALSE, transition_affected_demogr = FALSE, 
          env_stochas_type = "normal", noise = 1, fraction_SDD = FALSE, 
          fraction_LDD = FALSE, dispersal_constants = c(0.7, 0.7, 0.1, 
                                                        3), no_yrs, Ktype = "ceiling", K = NULL, Kweight = FALSE, 
          sumweight = FALSE) 
{
  require(sp)
  if (exists("BEMDEM")) 
    rm(BEMDEM, inherits = TRUE)
  if (is.vector(matrices)) {
    matrices <- matrix(matrices, ncol = 2, nrow = length(matrices))
    print("You are carrying out deterministic modelling.")
    colnames(matrices) <- c("matrixA", "matrixA")
  }
  if (length(proportion_initial) != length(stages)) 
    print("Number of stages or proportions is wrong!")
  if (nrow(matrices)%%length(stages) != 0) 
    print("Number of rows in matrix is not a multiple of stages name vector!")
  if (is.vector(Populations)) 
    print("There must be at least two populations!")
  if (sum(proportion_initial) > 1.02 | sum(proportion_initial) < 
      0.99) 
    print("Your 'proportion_initial' doesn't add to one...")
  if (is.numeric(sumweight)) {
    if (length(sumweight) != length(stages)) 
      print("Length of sumweight does not correpond to length of stages!")
  }
  list_names_matrices <- list()
  for (i in 1:ncol(matrices)) {
    M_name_one <- paste(colnames(matrices)[i], sep = "_")
    list_names_matrices <- c(list_names_matrices, list(M_name_one))
  }
  if (is.vector(Nichemap) | (Nichemap == "oneperiod")[1]) {
    min_dist <- sort(unique(dist(Populations[, 2:3])))[1]
    extent <- expand.grid(X = seq(min(Populations[, "X"]), 
                                  max(Populations[, "X"]), by = min_dist), Y = seq(min(Populations[, 
                                                                                                   "Y"]), max(Populations[, "Y"]), by = min_dist))
    Nichemap <- cbind(HScoreID = 1:nrow(extent), extent, 
                      matrix(1, ncol = length(Nichemap), nrow = nrow(extent), 
                             dimnames = list(NULL, paste(Nichemap))))
  }
  if (is.vector(Nichemap[, -c(1:3)])) {
    Nichemap <- Nichemap[Nichemap[, -c(1:3)] != 0, ]
  } else {
    #get rid of any cells which have no suitable cells in any year
    Nichemap <- Nichemap[rowSums(Nichemap[, -c(1:3)]) != 
                           0, ]
  }
  years_projections <- colnames(Nichemap)[4:ncol(Nichemap)]
  if ((ncol(Nichemap) - 3) != length(years_projections)) 
    print("Number of years of projections is not equal to the number of habitat scores!")
  colnames(Populations) <- c("PatchID", "XCOORD", "YCOORD", 
                             "area_population")
  colnames(Nichemap) <- c("HScoreID", "XCOORD", "YCOORD", years_projections)
  if (max(Nichemap[, 4:ncol(Nichemap)]) > 100) {
    Nichemap[, 4:ncol(Nichemap)] <- Nichemap[, 4:ncol(Nichemap)]/1000
  }
  if (max(Nichemap[, 4:ncol(Nichemap)]) > 10) {
    Nichemap[, 4:ncol(Nichemap)] <- Nichemap[, 4:ncol(Nichemap)]/100
  }
  Niche_ID <- data.frame(matrix(0, nrow = nrow(Nichemap), ncol = 4))
  Niche_ID[, 1:3] <- Nichemap[, 1:3]
  colnames(Niche_ID) <- c("Niche_ID", "X", "Y", "PopulationID")
  rownames(Niche_ID) <- Nichemap[, 1]
  if (length(density_individuals) == 1) {
    density_individuals <- rep(density_individuals, times = nrow(Populations))
  }
  n0_all <- matrix(0, nrow = nrow(Nichemap), ncol = length(stages))
  destination_Nicherows <- 1:nrow(Populations)
  for (pxs in 1:nrow(Populations)) {
    #distance is in kilometres as longlat = TRUE
    rows <- which(spDistsN1(as.matrix(Nichemap[, 2:3], ncol = 2), 
                            matrix(as.numeric(Populations[pxs, 2:3]), ncol = 2), 
                            longlat = FALSE) == min(spDistsN1(as.matrix(Nichemap[, 
                                                                                 2:3], ncol = 2), matrix(as.numeric(Populations[pxs, 
                                                                                                                                2:3]), ncol = 2), longlat = FALSE)))
    # rows <- which(spDistsN1(as.matrix(Nichemap[, 2:3], ncol = 2), 
    #                         matrix(as.numeric(Populations[pxs, 2:3]), ncol = 2), 
    #                         longlat = TRUE) == min(spDistsN1(as.matrix(Nichemap[, 
    #                                                                              2:3], ncol = 2), matrix(as.numeric(Populations[pxs, 
    #                                                                                                                             2:3]), ncol = 2), longlat = TRUE)))
    Niche_ID[rows, 4] <- Populations[pxs, 1]
    n0_all[rows[1], ] <- n0_all[rows[1], ] + (Populations[pxs, 
                                                          4] * proportion_initial * density_individuals[pxs])
    destination_Nicherows[pxs] <- rows[1]
  }
  Niche_values <- as.matrix(Nichemap[, 4:(length(years_projections) + 
                                            3)], ncol = length(years_projections))
  if (is.numeric(K)) {
    populationmax_all <- matrix(mean(K), ncol = length(years_projections), 
                                nrow = nrow(Nichemap))
    colnames(populationmax_all) <- years_projections
    rownames(populationmax_all) <- Niche_ID[, "Niche_ID"]
  }
  if (length(K) == 1) {
    populationmax_all <- matrix(K, ncol = length(years_projections), 
                                nrow = nrow(Nichemap))
  }
  if (length(K) == nrow(Populations)) {
    populationmax_all <- matrix(0, ncol = length(years_projections), 
                                nrow = nrow(Nichemap))
    for (rx in 1:length(destination_Nicherows)) {
      populationmax_all[destination_Nicherows[rx], ] <- populationmax_all[destination_Nicherows[rx], 
                                                                          ] + K[rx]
    }
    populationmax_all[populationmax_all == 0] <- mean(K)
  }
  if (length(K) == length(years_projections)) {
    populationmax_all[rowSums(n0_all) == 0, ] <- matrix(K, 
                                                        ncol = length(years_projections), nrow = nrow(Nichemap) - 
                                                          nrow(Populations))
    populationmax_all[rowSums(n0_all) > 0, ] <- matrix(K, 
                                                       ncol = length(years_projections), nrow = nrow(Populations), 
                                                       byrow = TRUE)
  }
  if (length(dim(K)) == 2) {
    populationmax_all[, ] <- matrix(colMeans(K), ncol = length(years_projections), 
                                    nrow = nrow(Nichemap), byrow = TRUE)
    populationmax_all[rowSums(n0_all) > 0, ] <- K
  }
  if (is.null(K)) {
    populationmax_all <- matrix("no_K", ncol = length(years_projections), 
                                nrow = nrow(Nichemap))
  }
  dist_populations <- spDists(as.matrix(Niche_ID[, 2:3]))
  dimnames(dist_populations) <- list(Niche_ID[, 1], Niche_ID[, 
                                                             1])
  # dispersal_probabilities <- dist_latlong <- neigh_index <- NA
  # if (dispersal_constants[1] != FALSE) {
  #   dispersal_probabilities <- dispersal_constants[1] * exp(-dist_populations^(dispersal_constants[3]/dispersal_constants[2]))
  #   dispersal_probabilities[dist_populations > dispersal_constants[4]] <- 0
  #   diag(dispersal_probabilities) <- 0
  # }
  # 
  #inserting a decay function with required parameters of median dispersal and maximum dispersal
  dispersal_probabilities <- dist_latlong <- neigh_index <- NA
  if (dispersal_constants[1] != FALSE) {
    half_life<-(log(2)/dispersal_constants[1]) * - 1
    dispersal_probabilities<-1 * exp(half_life * dist_populations)
    dispersal_probabilities[dist_populations > dispersal_constants[2]]<-0
    diag(dispersal_probabilities) <- 0
  }
  dist_latlong <- round(as.matrix(dist(Niche_ID[, 2:3])), 1)
  neigh_index <- sort(unique(as.numeric(dist_latlong)))[2:3] #distance two closest cells
  if (sumweight[1] == "all_stages") 
    sumweight <- rep(1, length(proportion_initial))
  if (Kweight[1] == "FALSE") 
    Kweight <- rep(1, length(proportion_initial))
  if (transition_affected_env[1] == "all") 
    transition_affected_env <- which(matrices[, 1] > 0)
  if (transition_affected_niche[1] == "all") 
    transition_affected_niche <- which(matrices[, 1] > 0)
  if (transition_affected_demogr[1] == "all") 
    transition_affected_demogr <- which(matrices[, 1] > 0)
  if (any(matrices < 0)) 
    print("There are some negative rates in the transition matrices!")
  if (any(matrices_var < 0)) 
    print("There are some negative rates in the standard deviation transition matrices!")
  if (max(transition_affected_niche) > nrow(matrices)) {
    print("Stages affected by Habitat suitability values does not comply with the size of matrix! Not that the matrix is made with 'byrow = FALSE")
  }
  if (max(transition_affected_env) > nrow(matrices)) {
    print("Stages affected by environmental stochasticity does not comply with the size of matrix! Note that the matrix is made with 'byrow = FALSE")
  }
  if (max(transition_affected_demogr) > nrow(matrices)) {
    print("Stages affected by demographic stochasticity does not comply with the size of matrix! Note that the matrix is made with 'byrow = FALSE")
  }
  BEMDEM <- list(Orig_Populations = Populations, Niche_ID = Niche_ID, 
                 Niche_values = Niche_values, years_projections = years_projections, 
                 matrices = matrices, matrices_var = matrices_var, prob_scenario = prob_scenario, 
                 noise = noise, stages = stages, proportion_initial = proportion_initial, 
                 density_individuals = density_individuals, fraction_LDD = fraction_LDD, 
                 fraction_SDD = fraction_SDD, dispersal_probabilities = dispersal_probabilities, 
                 dist_latlong = dist_latlong, neigh_index = neigh_index, 
                 no_yrs = no_yrs, K = K, Kweight = Kweight, populationmax_all = populationmax_all, 
                 n0_all = n0_all, list_names_matrices = list_names_matrices, 
                 sumweight = sumweight, transition_affected_env = transition_affected_env, 
                 transition_affected_niche = transition_affected_niche, 
                 transition_affected_demogr = transition_affected_demogr, 
                 env_stochas_type = env_stochas_type)
  assign(modelname, BEMDEM, envir = .GlobalEnv)
  eval(parse(text = paste("save(", modelname, ", file='", modelname, 
                          ".rda')", sep = "")))
}