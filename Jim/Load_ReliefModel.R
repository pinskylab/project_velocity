### Coastal Relief
setwd("~/Documents/Collaborations/Rutgers/External Layers/")

library(rgdal)
library(raster)
ne.raster <- raster("ne_atl_crm_v1.asc")

se.raster <- raster("/Users/abigailporay/Desktop/se_atl_crm_v1.asc")




dev.new()
plot(ne.raster)

dev.new()
plot(se.raster)
