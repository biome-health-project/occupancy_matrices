##

calcOcc <-
  function(species,   #the species name - in dataframe - that the function is to be run for
           d = d,     # the dataframe with species, site, and each date it was seen at that site - must have a columns called Species, Site and DateTime 
           all_cams = all_cams,    # a matrix with all the dates cameras, 1s for dates when the camera was working/deployed and NAs for when it wasn't
           startDate = startDate,    #start date in date format
           endDate = endDate) {
    # Make a vector of breaks 
    brks <- seq(startDate, endDate, by = "day")   #makes a sequence of all the days from start to end

        # Create an empty matrix of dim sites x time periods
    occ <-
      matrix(0, ncol= length(unique(d$Site)), nrow = length(brks)) 
    colnames(occ) <- sort(unique(d$Site))
    rownames(occ) <- strftime(brks, format = "%Y-%m-%d")
    
    for (s in unique(d$Site)) {     #this loops through each site and inserts 1's on days which there were captures  
      seen <- NA
      captures <- na.omit(d$DateTime[d$Species == species & d$Site == s]) 
      # Were animals seen at the site
      seen <- which(brks %in% captures)
      # If the species was seen, occ = 1

      col_i <- which(colnames(occ) == s)
      occ[seen,col_i] <- 1
    }
    
      occ <- occ * all_cams[, 2:ncol(all_cams)]
      print(paste0(species, " done!"))
      species_name <- gsub(" ", "", species)
      row.names(occ) <- brks
      write.csv(occ, here::here("matrices_out", paste0(species_name, "_tt_effort.csv")))
      
    return(occ)
    
    
  }



timestepper <- function(occ_in, timestep, na_mode = "include") {
  if (na_mode == "include") {
    occ_in[is.na(occ_in)] <- 0   #replacing NAs with 0s if we want to include them in analysis.
  }
  
  if (timestep > nrow(occ_in) / 2) {
    print(paste(
      "Time step is too large! Please reduce to",
      nrow(occ_in) / 2 ,
      "or less."
    ))
  } else {
    start <- seq(1, nrow(occ_in), by = timestep)
    end <- seq(timestep, nrow(occ_in), by = timestep)
    
    if (length(start) > length(end)) {
      start <- start[-length(start)]
    }
    
    timesteps <- matrix(nrow = length(start), ncol = ncol(occ_in))
    colnames(timesteps) <- colnames(occ_in)
    rownames(timesteps) <-
      paste(rownames(occ_in)[start], rownames(occ_in)[end], sep = ":")
    
    for (i in 1:length(start)) {
      timestep_out <- colSums(occ_in[start[i]:end[i],])
      timesteps[i,] <- timestep_out
      timesteps[timesteps > 0] <- 1
    }
    
    timesteps<-t(timesteps)
    return(timesteps)
    
  }
  
}


unmarker<-function(sp_matrix, species_name){
  umf.sp_matrix<-unmarkedFrameOccu(y=sp_matrix, siteCovs=covs)  
  m0_5d<-occu(~1 ~1,data=umf.sp_matrix) # null model
  m1_5d <- occu(~1 ~manag,data=umf.sp_matrix) # managemt type
  m2_5d<-occu(~1 ~DistPark_Norm,data=umf.sp_matrix) #distance to park
  
  m_list<-list(m0_5d, m1_5d, m2_5d)
  names(m_list)<-c(paste0(species_name, " m0"),
                   paste0(species_name, " m1"),
                   paste0(species_name, " m2"))
  return(m_list)
}



###function which gets predictions from a chosen model
predict_unmarked<-function(model, species_name) {
  # Predicting m2
  predPark<-predict(model,type="state",newdata=dist_seq,appendData=T)
  parkReal<- (predPark$DistPark*sdPark)+meanPark #formula=(z-value*SD of original value)+Mean of original values
  predParkReal<- cbind(species_name,parkReal,predPark) #adding the real values as a collumn to the predicted dataframe  
  return(predParkReal)
}


