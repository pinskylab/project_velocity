# These first 3 lines are used to make arguments for using a shell script on the cluster system_if not using a shell script they are not necessary
args <- commandArgs(trailingOnly=T)
print(args)
myalcnt1 <- as.numeric(args[1])
 
require(mgcv)
 
projfolder <- 'PAmodels_proj_May2018/'
RCP <- c(26,45,85)
season <- 'jas'
# Below there are two new climate models that could be added here
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')

load('speciesProjectionList.RData') # A list of species with corresponding regionfact values for each species 

# Begin loop to make a prediction for each species and save each species individually
predsE <- matrix(data=NA, nrow=7763648, ncol=18)
predsW <- matrix(data=NA, nrow=6187644, ncol=18)

  species = projspp[myalcnt1]
  load(paste('GAM_models/CEmods_Nov2017_fitallreg_2017_', species, '.RData', sep=''))
  mygam1 <- mods[[1]]
  rm(mods)
  for(i in 1:3){
    rcp <- RCP[i]
    print(paste('Beginning predictions for rcp', rcp, '_', species, '_', Sys.time(), sep=''))
    for(j in 1:length(modelrun)){
      if(grepl('_Atl', species)){
        load(paste('Prediction_files/predictionEASTnoBias_rcp',rcp, '_', season, '_', modelrun[j], '.RData', sep=''))
        predictSpp <- pred.bathE
      } else{
        load(paste('Prediction_files/predictionWESTnoBias_rcp',rcp, '_', season, '_', modelrun[j], '.RData', sep=''))
        predictSpp <- pred.bathW
      }
      predictSpp$regionfact = regionFreq[myalcnt1]
      print(paste('Predicting for: ', species, '    model number:', j, ' ', modelrun[j], '    ', Sys.time(), sep=''))
      preds1 <- predict(mygam1, newdata = predictSpp, type='response')
      if(grepl('_Atl', species)){
        predsE[,j] <- preds1
      } else{
        predsW[,j] <- preds1
      }
    }    
    predictSpp <- predictSpp[,c(1:4, 10:12)]
    if(grepl('_Atl', species)){
      pred.grid <- data.frame(cbind(predictSpp, predsE))
    } else{
      pred.grid <- data.frame(cbind(predictSpp, predsW))
    }
       
    # Aggregate predictions into the 20year bins for plotting_saving
    pred.agg <- aggregate(list(mean1=pred.grid$X1, mean2=pred.grid$X2, mean3=pred.grid$X3, mean4=pred.grid$X4, mean5=pred.grid$X5, mean6=pred.grid$X6, mean7=pred.grid$X7, 
                mean8=pred.grid$X8, mean9=pred.grid$X9, mean10=pred.grid$X10, mean11=pred.grid$X11, mean12=pred.grid$X12, mean13=pred.grid$X13, mean14=pred.grid$X14, 
                mean15=pred.grid$X15, mean16=pred.grid$X16, mean17=pred.grid$X17, mean18=pred.grid$X18),
                by=list(year_range=pred.grid$bin, latitude=pred.grid$latBathgrid, longitude=pred.grid$lonBathgrid), FUN=mean)
      
    filenameAgg <- paste(projfolder, species, '_rcp', rcp, '_', season, '_PresAbs_projection_AGG.RData', sep='')
    save(pred.agg, file=filenameAgg)
  }
   