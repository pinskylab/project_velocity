library("rgdal") # for `ogrInfo()` and `readOGR()`
library("tools") # for `file_path_sans_ext()`
library("dplyr") # for `inner_join()`, `filter()`, `summarise()`, and the pipe operator (%>%)
library("ggplot2") # for `fortify()` and for plotting
library("sp") # for `point.in.polygon()` and `spDists()`
library("tidyr") # for `gather()`
library("readr") # for `write_tsv()`
 
# Below is a function I got off the web to convert a shapefile into a data frame
fortify.shape <- function(x){
  x@data$id <- rownames(x@data)
  x.f <- fortify(x, region = "id")
  x.join <- inner_join(x.f, x@data, by = "id")
}

# Import the shapefile of the U.S. EEZ
path.eez.usa <- ("/Users/jamesmorley/Desktop/USMaritimeLimitsAndBoundariesSHP")
fnam.eez.usa <- "USMaritimeLimitsNBoundaries.shp"
eez.usa <- readOGR(dsn = path.eez.usa, layer = file_path_sans_ext(fnam.eez.usa))
dat.eez.usa1 <- fortify.shape(eez.usa) # make into a dataframe

dat.eez.usa1MOD <- dat.eez.usa1[dat.eez.usa1$REGION=='Atlantic Coast and Gulf of Mexico',]
dat.eez.usa1MOD <- dat.eez.usa1MOD[dat.eez.usa1MOD$FEAT_TYPE == 'Maritime Boundary',]
dat.eez.usa1MOD <- dat.eez.usa1MOD[dat.eez.usa1MOD$lat > 35,]

# FIRST GET THE REFERENCE COLUMNS AND CARVE OUT THE MAJORITY OF IT WITH SIMPLE ifelse STATEMENTS
# IMPORT pred.agg FOR AN EAST COAST SPECIES FROM THE 'EPA_data.R' script
pred.aggEEZ <- pred.agg[,c(1:3)]
pred.aggEEZ$EEZ <- TRUE
# Carve out Mexico_small area 
pred.aggEEZ$EEZ[pred.aggEEZ$latitude < 25.99 & pred.aggEEZ$longitude < -90] <- FALSE
# Carve out most of Canada
pred.aggEEZ$EEZ[pred.aggEEZ$longitude > -65.70] <- FALSE
pred.aggEEZ$EEZ[pred.aggEEZ$latitude > 45.4] <- FALSE
# Identify which rows need to be sorted out in the for loop
pred.aggEEZ$sort_region <- FALSE
pred.aggEEZ$sort_region[pred.aggEEZ$EEZ==TRUE & pred.aggEEZ$longitude > -67.74] <- TRUE

 
EEZ <- matrix(data = NA, nrow = nrow(pred.aggEEZ), ncol = 1)
options(warn=2)
for(i in 1:nrow(pred.aggEEZ)){
  EEZ[i,1] <- pred.aggEEZ$EEZ[i] # default to what is already assigned
  if(pred.aggEEZ$sort_region[i]==TRUE){
    abc <- dat.eez.usa1MOD  
    abc$lat_test <- pred.aggEEZ$latitude[i]
    abc$diff <- abs(abc$lat - abc$lat_test)
    lon <- abc$long[abc$diff == min(abc$diff)[1]][1]
    
    if(pred.aggEEZ$longitude[i] > lon){
      EEZ[i,1] <- FALSE
    }
  }  
}


pred.aggEEZ$EEZ_update <- EEZ
plot(latitude~longitude, cex=.01, data=pred.aggEEZ[pred.aggEEZ$year_range=='2007-2020',])
points(latitude~longitude, cex=.01, col='green', data=pred.aggEEZ[pred.aggEEZ$year_range=='2007-2020' & pred.aggEEZ$EEZ_update==TRUE,])

points(latitude~longitude, cex=.01, col='green', data=pred.aggEEZ[pred.aggEEZ$year_range=='2007-2020' & pred.aggEEZ$EEZ ==TRUE,])
points(lat~long, col='red', cex=.05, dat.eez.usa1MOD)
points(latitude~longitude, col='blue', cex=.01, data=pred.aggEEZ[pred.aggEEZ$sort_region==TRUE,])
points(lat~long, col='red', cex=.05, dat.eez.usa1MOD)

#identical(pred.agg[,c(2:3)], pred.aggEEZ[,c(2:3)])
EEZ_east <- EEZ
save(EEZ_east, file='data/EEZ_grid_east.RData')


# NOW THE WEST COAST ==============================================================================
# FIRST GET THE REFERENCE COLUMNS AND CARVE OUT THE MAJORITY OF IT WITH SIMPLE ifelse STATEMENTS
# IMPORT pred.agg FOR A WEST COAST SPECIES FROM THE 'EPA_data.R' script
pred.aggEEZ <- pred.agg[,c(1:3)]
pred.aggEEZ$EEZ <- TRUE
# Carve out Mexico_small area 
pred.aggEEZ$EEZ[pred.aggEEZ$latitude < 32.59 & pred.aggEEZ$longitude > -118] <- FALSE

# Need to specify three boundary lines, which I'll have to numerically code in the forloop to identify somehow
# First do Canada-Washington border
dat.eez.usa1MOD <- dat.eez.usa1[dat.eez.usa1$REGION=='Pacific Coast',]
dat.eez.usa1MOD <- dat.eez.usa1MOD[dat.eez.usa1MOD$FEAT_TYPE == 'Maritime Boundary',]
dat.eez.usa1MOD <- dat.eez.usa1MOD[dat.eez.usa1MOD$lat > 42 & dat.eez.usa1MOD$long < -124.5,]
# Second make Alaska-Canada border and Bering Sea border
dat.eez.usa1MODb <- dat.eez.usa1[dat.eez.usa1$REGION=='Alaska',]
dat.eez.usa1MODb <- dat.eez.usa1MODb[dat.eez.usa1MODb$FEAT_TYPE == 'Maritime Boundary',]
dat.eez.usa1MODb <- dat.eez.usa1MODb[!dat.eez.usa1MODb$lat > 62.125,]
dat.eez.usa1MODb <- dat.eez.usa1MODb[!(dat.eez.usa1MODb$lat < 60 & dat.eez.usa1MODb$lat > 57.5),]
# Separate the Bering sea border
dat.eez.usa1MODc <- dat.eez.usa1MODb[!dat.eez.usa1MODb$lat < 58,]
dat.eez.usa1MODc <- dat.eez.usa1MODc[!dat.eez.usa1MODc$lat < 60.2,]
dat.eez.usa1MODb <- dat.eez.usa1MODb[!dat.eez.usa1MODb$lat > 58,]
dat.eez.usa1MODb <- dat.eez.usa1MODb[!dat.eez.usa1MODb$long > 0,]

#NOW BASED ON THESE BOUNDARY LINES, CARVE OUT MUCH OF CANADA
pred.aggEEZ$EEZ[pred.aggEEZ$latitude > 48.515 & pred.aggEEZ$latitude < 53.465 & pred.aggEEZ$longitude > -150] <- FALSE

# Identify which rows need to be sorted out in the for loop with a numeric value to indicate which border region
pred.aggEEZ$sort_region <- 0
pred.aggEEZ$sort_region[pred.aggEEZ$EEZ==TRUE & pred.aggEEZ$latitude > 46.53 & pred.aggEEZ$latitude < 48.51] <- 1
pred.aggEEZ$sort_region[pred.aggEEZ$EEZ==TRUE & pred.aggEEZ$latitude > 53.47 & pred.aggEEZ$latitude < 54.77 & pred.aggEEZ$longitude > -150] <- 2
pred.aggEEZ$sort_region[pred.aggEEZ$EEZ==TRUE & pred.aggEEZ$latitude > 60 & pred.aggEEZ$longitude < -176.2] <- 3


EEZ <- matrix(data = NA, nrow = nrow(pred.aggEEZ), ncol = 1)
options(warn=2) # stops the loop if an error occurs
for(i in 1:nrow(pred.aggEEZ)){
  EEZ[i,1] <- pred.aggEEZ$EEZ[i] # default to what is already assigned
  # Now three logic blocks, one for each boundary region
  if(pred.aggEEZ$sort_region[i]==1){
    abc <- dat.eez.usa1MOD  
    abc$lon_test <- pred.aggEEZ$longitude[i]
    abc$diff <- abs(abc$long - abc$lon_test)
    lat <- abc$lat[abc$diff == min(abc$diff)[1]][1]
    if(pred.aggEEZ$latitude[i] > lat){
      EEZ[i,1] <- FALSE
    }
  }  
  if(pred.aggEEZ$sort_region[i]==2){
    abc <- dat.eez.usa1MODb  
    abc$lon_test <- pred.aggEEZ$longitude[i]
    abc$diff <- abs(abc$long - abc$lon_test)
    lat <- abc$lat[abc$diff == min(abc$diff)[1]][1]
    if(pred.aggEEZ$latitude[i] < lat){
      EEZ[i,1] <- FALSE
    }
  }  
  if(pred.aggEEZ$sort_region[i]==3){
    abc <- dat.eez.usa1MODc  
    abc$lon_test <- pred.aggEEZ$longitude[i]
    abc$diff <- abs(abc$long - abc$lon_test)
    lat <- abc$lat[abc$diff == min(abc$diff)[1]][1]
    if(pred.aggEEZ$latitude[i] > lat){
      EEZ[i,1] <- FALSE
    }
  }  
}

pred.aggEEZ$EEZ_update <- EEZ
plot(latitude~longitude, cex=.01, data=pred.aggEEZ[pred.aggEEZ$year_range=='2007-2020',])
points(latitude~longitude, cex=.01, col='green', data=pred.aggEEZ[pred.aggEEZ$year_range=='2007-2020' & pred.aggEEZ$EEZ_update==TRUE,])

#points(latitude~longitude, cex=.01, col='blue', data=pred.aggEEZ[pred.aggEEZ$year_range=='2007-2020' & pred.aggEEZ$EEZ==T & pred.aggEEZ$sort_region > 0,])
points(lat~long, col='red', cex=.01, dat.eez.usa1MOD)
points(lat~long, col='red', cex=.01, dat.eez.usa1MODb)
points(lat~long, col='red', cex=.01, dat.eez.usa1MODc)

EEZ_west <- EEZ
save(EEZ_west, file='data/EEZ_grid_west.RData')
 


