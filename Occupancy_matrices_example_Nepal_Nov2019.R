# Occupancy model example for Nepal data
# Models used in the November report
library(unmarked)
library(here)
library(dplyr)

# function to aggregate survey days into longer occasions
source("occupancyfunctions.R") #all of the functions are in here so the below code does not need to be included, it just is for transparency

setwd(
  "Z:/biome_health_project_files/country_files/nepal/tagging_photos/targeted_tagging_effort"
)

###
##Functions in occupyfunctions.R - it's important that you use the same species name throughout
####  calcOcc function
calcOcc <-
  function(species,
           #the species name - in dataframe - that the function is to be run for
           d = d,
           # the dataframe with species, site, and each date it was seen at that site - must have a columns called Species, Site and DateTime
           all_cams = all_cams,
           # a matrix with all the dates cameras, 1s for dates when the camera was working/deployed and NAs for when it wasn't
           startDate = startDate,
           #start date in date format
           endDate = endDate) {
    # Make a vector of breaks
    brks <-
      seq(startDate, endDate, by = "day")   #makes a sequence of all the days from start to end
    
    # Create an empty matrix of dim sites x time periods
    occ <-
      matrix(0, ncol = length(unique(d$Site)), nrow = length(brks))
    colnames(occ) <- sort(unique(d$Site))
    rownames(occ) <- strftime(brks, format = "%Y-%m-%d")
    
    for (s in unique(d$Site)) {
      #this loops through each site and inserts 1's on days which there were captures
      seen <- NA
      captures <-
        na.omit(d$DateTime[d$Species == species & d$Site == s])
      # Were animals seen at the site
      seen <- which(brks %in% captures)
      # If the species was seen, occ = 1
      
      col_i <- which(colnames(occ) == s)
      occ[seen, col_i] <- 1
    }
    
    occ <- occ * all_cams[, 2:ncol(all_cams)]
    print(paste0(species, " done!"))
    species_name <- gsub(" ", "", species)
    row.names(occ) <- brks
    write.csv(occ, here::here("matrices_out", paste0(species_name, "_tt_effort.csv")))
    
    return(occ)
    
    
  }


##timestepper - creates matrices of a given timestep, can choose to include or exclude NAs
timestepper <- function(occ_in, timestep, na_mode = "include") {
  if (na_mode == "include") {
    occ_in[is.na(occ_in)] <-
      0   #replacing NAs with 0s if we want to include them in analysis.
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
      timestep_out <- colSums(occ_in[start[i]:end[i], ])
      timesteps[i, ] <- timestep_out
      timesteps[timesteps > 0] <- 1
    }
    
    timesteps <- t(timesteps)
    return(timesteps)
    
  }
  
}
###unmarker
### runs functions from the unmarked package to produce models - outputs a list of models
unmarker <- function(sp_matrix, species_name) {
  umf.sp_matrix <- unmarkedFrameOccu(y = sp_matrix, siteCovs = covs)
  m0_5d <- occu( ~ 1 ~ 1, data = umf.sp_matrix) # null model
  m1_5d <- occu( ~ 1 ~ manag, data = umf.sp_matrix) # managemt type
  m2_5d <-
    occu( ~ 1 ~ DistPark_Norm, data = umf.sp_matrix) #distance to park
  
  m_list <- list(m0_5d, m1_5d, m2_5d)
  names(m_list) <- c(
    paste0(species_name, " m0"),
    paste0(species_name, " m1"),
    paste0(species_name, " m2")
  )
  return(m_list)
}
#####
###function which gets predictions from a chosen model
predict_unmarked <- function(model, species_name) {
  # Predicting m2
  predPark <-
    predict(model,
            type = "state",
            newdata = dist_seq,
            appendData = T)
  parkReal <-
    (predPark$DistPark * sdPark) + meanPark #formula=(z-value*SD of original value)+Mean of original values
  predParkReal <-
    cbind(species_name, ParkReal, predPark) #adding the real values as a collumn to the predicted dataframe
  return(predParkReal)
}
#####


#### Testing effect of management and distance from park on species occupancy ####
#loading variables
all_cams <-
  read.csv("Surevey_effort_30minthreshold.csv", stringsAsFactors = FALSE)
species_df <-
  read.csv("Latest_species_meta.csv", stringsAsFactors = FALSE)
covs <- read.csv("site_covs.csv") # loading covariates
covs$DistPark_Norm <- scale(covs$DistPark)
sdPark <- sd(covs$DistPark)
meanPark <- mean(covs$DistPark)
dist_seq = data.frame(DistPark_Norm = seq(min(covs$DistPark_Norm), max(covs$DistPark_Norm), by =
                                            0.1)) # must be values of the scaled variable

#######
#replacing 0s with NAs - NAs indicate when the cameras were not on
#We want to build a matrix for each species where 1 = species detected, 0 = camera on but species not detected, NA =camera not on
all_cams[all_cams == 0] <- NA

d <- species_df %>%
  select(CommonName, site_cam.x , date_fixed) %>%
  arrange(CommonName) %>%
  distinct()

d$date_fixed <-
  as.Date(d$date_fixed, format = "%d/%m/%Y")
colnames(d) <- c("Species", "Site", "DateTime")


all_cams <- all_cams %>%
  select(-c("OBZ03", "OBZ12")) #getting rid of columns with cameras with no sp detections


if (!dir.exists(here::here("matrices_out"))) {
  dir.create(here::here("matrices_out"))
}

# This lapply function will create a effort matrix for each species using the calcOcc function - in occupancyfunctions.R

lapply(
  X = unique(species_df$CommonName),
  FUN = calcOcc,
  d = d,
  all_cams = all_cams,
  startDate = as.Date("2019-03-15"),
  endDate = as.Date("2019-04-15")
)


### Building matrices with 5-day survey occasion for Nepal 2019 data ####
# it uses csv files with detection/no detection for each species per day (31 days x 148 sites)

# 5-day Wildlife matrices

#BarkingDeer
BarkingDeer <-
  read.csv("BarkingDeer_tt_effort.csv") # read matrix with det/no det for each survey day
BarkingDeer <- BarkingDeer[, -1]
#5-day survey occasion matrix
BarkingDeer_5d <- timestepper(occ_in = BarkingDeer,
                              timestep = 5,
                              na_mode = "include")
sum(BarkingDeer_5d)

bd_out <- unmarker(BarkingDeer_5d, species_name = "BarkingDeer")
#ways of looking at the list of outputted models
bd_out ##full list of models
bd_out[[1]]  #first model in list
bd_out[["BarkingDeer m0"]] #also full model in list but extracted by name
bd_out[["BarkingDeer m1"]]
bd_out[["BarkingDeer m2"]]

predict_unmarked(bd_out[["BarkingDeer m2"]], species_name = "BarkingDeer")


#Chital
Chital <-
  read.csv("Chital_tt_effort.csv") # read matrix with det/no det for each survey day
Chital <- Chital[, -1]
names(Chital)
#5-day survey occasion matrix
Chital_5d <- timestepper(occ_in = Chital,
                         timestep = 5,
                         na_mode = "include")

ch_out <- unmarker(Chital_5d, species_name = "Chital")
predict_unmarked(ch_out[["Chital m2"]], species_name = "Chital")


#GreyLangur
GreyLangur <-
  read.csv("GreyLangur_tt_effort.csv") # read matrix with det/no det for each survey day
GreyLangur <- GreyLangur[, -1]
#5-day survey occasion matrix
GreyLangur_5d <- timestepper(occ_in = GreyLangur,
                             timestep = 5,
                             na_mode = "include")

gl_out <-
  unmarker(GreyLangur_5d, species_name = "GreyLangur")   #species name should be the ssame throughout
predict_unmarked(gl_out[["GreyLangur m2"]], species_name = "GreyLangur")
