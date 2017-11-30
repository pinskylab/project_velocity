args <- commandArgs(trailingOnly=T)
print(args)
myalcnt1 <- as.numeric(args[1])
myalcnt2 <- as.numeric(args[2])
 
load('data/speciesProjectionList.RData')
# below is a list of species that were already done and need to be removed from projspp
doneSpp <- sort(projspp[c(419,600,105,310,523,517,270,438,587,465,525,526,392,17,79,135,196,26,28,498,307,164,275,37,66,126,136,151,163,170,188,264)])
projspp <- projspp[!(projspp %in% doneSpp)]
regionFreq <- regionFreq[-c(419,600,105,310,523,517,270,438,587,465,525,526,392,17,79,135,196,26,28,498,307,164,275,37,66,126,136,151,163,170,188,264)]
  
require(mgcv)

RCP <- c(26,85)
season <- 'jas'
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MPI-ESM-LR','NorESM1-ME')
pred.metric <- c('max', 'min', 'mean')

predsE <- matrix(data=NA, nrow=7763648, ncol=16)
predsW <- matrix(data=NA, nrow=6187644, ncol=16)

sppSum <- data.frame(coefnt = numeric())
for(k in myalcnt1:myalcnt2){
  species = projspp[k]
  load(paste('output/CEmodels/CEmods_fitallreg_2017_', species, '.RData', sep=''))
  mygam1 <- mods[[1]]
  sppSum[k,] <- data.frame(coefnt = coef(mygam1)[2])
  print(paste("Finished ", species, sep=""))  
}
write.csv(sppSum, file="test_shell.csv")
print('FINISHED')   
  
  
  
  
  
  