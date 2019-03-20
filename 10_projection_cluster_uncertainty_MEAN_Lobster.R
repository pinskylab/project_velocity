# MUTED LINES ARE FOR MY PERSONAL CPU
library(mgcv)
library(brglm)
library(dismo)
library(gbm)
library(MASS) 
       
projfolder <- 'Mean_niche_proj_Dec2018/'
   
RCP <- c(26,45,85)
season <- 'jas'
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')
mods <- c('GLM', 'GAM', 'BRT')
years <- c("2007-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")
species <- 'homarus americanus_Atl'
    
if(grepl('_Atl', species)){
  predsE_PA <- matrix(data=NA, nrow=7763648, ncol=18)
  predsE_biom <- matrix(data=NA, nrow=7763648, ncol=18)
} else{
  predsW_PA <- matrix(data=NA, nrow=6187644, ncol=18)
  predsW_biom <- matrix(data=NA, nrow=6187644, ncol=18)
}

for(i in 1:3){
  nicheMod <- mods[i]
  load(paste('Niche_mods/nicheMods_', nicheMod, 's_', species, '_uncertainty2018.RData', sep=''))# load the niche model
  # Load the niche model
  if(nicheMod=='GLM'){
    PAmod <- myglm1
    BIOMmod <- myglm2
    rm(myglm1tt, myglm2tt, myglm1, myglm2)
  }
  if(nicheMod=='GAM'){
    PAmod <- mygam1
    BIOMmod <- mygam2
    rm(mygam1tt, mygam2tt, mygam1, mygam2)
  }
  if(nicheMod=='BRT'){
    PAmod <- myBRT1
    BIOMmod <- myBRT2
    rm(myBRT1tt, myBRT2tt, myBRT1, myBRT2)
  }
  
  for(j in 1:3){
    rcp <- RCP[j]
    
    for(k in 1:18){
      gcm <- modelrun[k]
      print(paste('Beginning predictions for rcp', rcp, '_', gcm, '_', species, '_', nicheMod, '_', Sys.time(), sep=''))
      # load climate prediction file
      if(grepl('_Atl', species)){
        load(paste('Prediction_files/predictionEASTnoBias_rcp',rcp, '_', season, '_', gcm, '.RData', sep=''))
        #load(paste('data/prediction_files_Nov2017/predictionEASTnoBias_rcp',rcp, '_', season, '_', gcm, '.RData', sep=''))
        predictSpp <- pred.bathE
        rm(pred.bathE)
      } else{
        load(paste('Prediction_files/predictionWESTnoBias_rcp',rcp, '_', season, '_', gcm, '.RData', sep=''))
        #load(paste('data/prediction_files_Nov2017/predictionWESTnoBias_rcp',rcp, '_', season, '_', gcm, '.RData', sep=''))
        predictSpp <- pred.bathW
        rm(pred.bathW)
      }
      predictSpp$regionfact = as.factor(names(myregions)[myregions==max(myregions)][1])
      
      if(nicheMod=='GLM' | nicheMod=='GAM'){
        preds1 <- predict(PAmod, newdata = predictSpp, type='response')
        preds2 <- exp(predict(BIOMmod, newdata = predictSpp, type='response')) 
      }  
      if(nicheMod=='BRT'){
        preds1 <- predict.gbm(PAmod, newdata = predictSpp, type='response', n.trees=PAmod$n.trees)
        preds2 <- exp(predict.gbm(BIOMmod, newdata = predictSpp, type='response', n.trees=BIOMmod$n.trees))
      }
      
      if(grepl('_Atl', species)){
        predsE_PA[,k] <- preds1
        predsE_biom[,k] <- preds1*preds2
      } else{
        predsW_PA[,k] <- preds1
        predsW_biom[,k] <- preds1*preds2
      }
    }
  
    predictSpp <- predictSpp[,c(10:12)]
    if(grepl('_Atl', species)){
      pred.grid_PA <- data.frame(cbind(predictSpp, predsE_PA))
      pred.grid_BIOM <- data.frame(cbind(predictSpp, predsE_biom))
    } else{   
      pred.grid_PA <- data.frame(cbind(predictSpp, predsW_PA))
      pred.grid_BIOM <- data.frame(cbind(predictSpp, predsW_biom))
    }
    
    # Aggregate predictions into 20year bins for plotting_saving
    pred.agg_PA <- aggregate(list(mean1=pred.grid_PA$X1, mean2=pred.grid_PA$X2, mean3=pred.grid_PA$X3, mean4=pred.grid_PA$X4, mean5=pred.grid_PA$X5, mean6=pred.grid_PA$X6, mean7=pred.grid_PA$X7, 
            mean8=pred.grid_PA$X8, mean9=pred.grid_PA$X9, mean10=pred.grid_PA$X10, mean11=pred.grid_PA$X11, mean12=pred.grid_PA$X12, mean13=pred.grid_PA$X13, mean14=pred.grid_PA$X14, 
            mean15=pred.grid_PA$X15, mean16=pred.grid_PA$X16, mean17=pred.grid_PA$X17, mean18=pred.grid_PA$X18), by=list(year_range=pred.grid_PA$bin, latitude=pred.grid_PA$latBathgrid, longitude=pred.grid_PA$lonBathgrid), FUN=mean)
    pred.agg_biom <- aggregate(list(mean1=pred.grid_BIOM$X1, mean2=pred.grid_BIOM$X2, mean3=pred.grid_BIOM$X3, mean4=pred.grid_BIOM$X4, mean5=pred.grid_BIOM$X5, mean6=pred.grid_BIOM$X6, mean7=pred.grid_BIOM$X7, 
            mean8=pred.grid_BIOM$X8, mean9=pred.grid_BIOM$X9, mean10=pred.grid_BIOM$X10, mean11=pred.grid_BIOM$X11, mean12=pred.grid_BIOM$X12, mean13=pred.grid_BIOM$X13, mean14=pred.grid_BIOM$X14, 
            mean15=pred.grid_BIOM$X15, mean16=pred.grid_BIOM$X16, mean17=pred.grid_BIOM$X17, mean18=pred.grid_BIOM$X18), by=list(year_range=pred.grid_BIOM$bin, latitude=pred.grid_BIOM$latBathgrid, longitude=pred.grid_BIOM$lonBathgrid), FUN=mean)
      
    filenameAgg <- paste(projfolder, species, '_rcp', rcp, '_', season, '_', nicheMod, '_proj_AGG_2018.RData', sep='')
    save(pred.agg_PA, pred.agg_biom, file=filenameAgg)
  }
}
  
  