library(dplyr)
source("occupancyfunctions.R")

all_cams <-
  read.csv("Surevey_effort_30minthreshold.csv", stringsAsFactors = FALSE)

species_df <-
  read.csv("Latest_species_meta.csv", stringsAsFactors = FALSE)

#replacing 0s with NAs - NAs indicate when the cameras were not on
#We want to build a matrix for each species where 1 = species detected, 0 = camera on but species not detected, NA =camera not on
all_cams[all_cams == 0] <- NA

sp_dates <- species_df %>%
  select(CommonName, site_cam.x , date_fixed) %>%
  arrange(CommonName) %>%
  distinct()

sp_dates$date_fixed <-
  as.Date(sp_dates$date_fixed, format = "%d/%m/%Y")

# no_sp <-
#   which(!colnames(all_cams)[2:ncol(all_cams)] %in% unique(sp_dates$site_cam.x)) #which cams are these?


all_cams<-all_cams %>% 
  select(-c("OBZ03", "OBZ12")) #getting rid of columns with cameras with no sp detections


d <- sp_dates
colnames(d) <- c("Species", "Site", "DateTime")

if(!dir.exists(here::here("matrices_out"))){
  dir.create(here::here("matrices_out"))
}


# This lapply function will create a effort matrix for each species using the calcOcc function - in occupancyfunctions.R

lapply(
  X = unique(species_df$CommonName),
  FUN = calcOcc,
  d = d,
  all_cams,
  startDate = as.Date("2019-03-15"),
  endDate = as.Date("2019-04-15")
)

#If you just want it for one species use this:

calcOcc(
  species = "Barking Deer",
  d = d,
  all_cams = all_cams,
  startDate = as.Date("2019-03-15"),
  endDate = as.Date("2019-04-15")
)


#####This section for compressing the matrices into difference time chunks####


chital <-
  read.csv(here::here("matrices_out", "Chital_tt_effort.csv"))

row.names(chital) <- chital$X

chital <- chital[,-1]

#na_mode = "include" means that NAs will effectively be treated as zeros.
#if na_mode = anything apart from "include" an NA in a time step will count the whole timestep as NA

#have just done an example with the Chital data but could set it up as above to create for all sp.

source("timestepper.R")


timestepper(occ_in = chital,
            timestep = 10,
            na_mode = "include")


timestepper(occ_in = chital,
            timestep = 4,
            na_mode = "exclude")
