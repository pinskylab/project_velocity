args <- commandArgs(trailingOnly=T)
print(args)
myalcnt1 <- as.numeric(args[1]) # for RCP
myalcnt2 <- as.numeric(args[2]) # for GCM
     
# MUTED LINES ARE FOR MY PERSONAL CPU
library(mgcv)
projfolder <- 'CEmodels_projUncer_Mar2018/'
#setwd('/Users/jamesmorley/Documents/project_velocity')
#projfolder <- 'output/CEmodels_Uncertainty_2018/'
  
RCP <- c(26,45,85)
season <- 'jas'
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')
nicheMod <- 'GAM'
years <- c("2007-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")
species <- 'sebastes alutus_Pac'
 
rcp = RCP[myalcnt1]
modelGCM = modelrun[myalcnt2]
 
# Load niche model
load(paste('niche_models/Cluster_nicheMods_', nicheMod, '_', species, '_uncer2018.RData', sep=''))
#load(paste('output/CEmodels_Uncertainty_2018/Cluster_nicheMods_', nicheMod, '_', species, '_uncer2018.RData', sep=''))
print(paste('Beginning predictions for rcp', rcp, '_', modelGCM, '_', species, '_', nicheMod, '_', Sys.time(), sep=''))
    
      if(grepl('_Atl', species)){
        load(paste('Prediction_files/predictionEASTnoBias_rcp',rcp, '_', season, '_', modelGCM, '.RData', sep=''))
        #load(paste('data/prediction_files_Nov2017/predictionEASTnoBias_rcp',rcp, '_', season, '_', modelGCM, '.RData', sep=''))
        predictSpp <- pred.bathE
        rm(pred.bathE)
      } else{
        load(paste('Prediction_files/predictionWESTnoBias_rcp',rcp, '_', season, '_', modelGCM, '.RData', sep=''))
        #load(paste('data/prediction_files_Nov2017/predictionWESTnoBias_rcp',rcp, '_', season, '_', modelGCM, '.RData', sep=''))
        predictSpp <- pred.bathW
        rm(pred.bathW)
      }
      predictSpp$regionfact = names(myregions)[myregions==max(myregions)][1]

        # a prediction matrix, needed for ultimately getting variance estimates for quantities derived from the model
        # each row is corresponding to a row of prediction data set, these values get multiplied by individual curve coefficient
        # values and then those values get summed to get a prediction for that row of data
        
        # Need to do each year bin separate b/c files are too big otherwise_at least for my CPU
        PA_final <- NULL
        biom_final <- NULL
        for(i in 1:5){ # five time periods
          print(paste('Starting time period ', i))
          year_bin <- years[i]
          Xp.1 <- predict(mod1, predictSpp[predictSpp$bin==year_bin,], type='lpmatrix')
          Xp.2 <- predict(mod2, predictSpp[predictSpp$bin==year_bin,], type='lpmatrix')
          ilink <- family(mod1)$linkinv

          if(grepl('_Atl', species)){
            PA_matrix <- matrix(nrow=82592, ncol=40)
            biom_matrix <- matrix(nrow=82592, ncol=40)
          }else{
            PA_matrix <- matrix(nrow=65826, ncol=40)
            biom_matrix <- matrix(nrow=65826, ncol=40)
          }
          
          for(k in 1:40){# 40 iterations
            print(paste('iteration ', k))
            pred.a1 <- Xp.1%*%br.1[k,] # combine mean predictions from Xp with resamples of parameter estimates_ the %*% is a matrix algebra operator I believe
            pred.a1 <- ilink(pred.a1) # mod1 uses the logit link function_here the 'ilink' function transforms it back to the scale of the response variable (i.e. a 0-1 value)
            pred.a2 <- Xp.2%*%br.2[k,]
            pred.a2 <- pred.a1 * exp(pred.a2) # the ilink function is not needed for biomass model (mod2) as it does not use the link function (gaussian distribution)
            pred.a1a2 <- data.frame(cbind(latBathgrid=predictSpp$latBathgrid[predictSpp$bin==year_bin], lonBathgrid=predictSpp$lonBathgrid[predictSpp$bin==year_bin], PA=as.vector(pred.a1), biom=as.vector(pred.a2)))
            # aggregate to one value for this year-bin
            a1a2_agg <- aggregate(list(PA=pred.a1a2$PA, biom=pred.a1a2$biom), by=list(latitude=pred.a1a2$latBathgrid, longitude=pred.a1a2$lonBathgrid), FUN=mean)
            PA_matrix[,k] <- a1a2_agg$PA
            biom_matrix[,k] <- a1a2_agg$biom
          } 
          # Need to stack the 5 time chunks here so it doesn't get replaced
          PA_matrix <- data.frame(cbind(PA_matrix, latitude=a1a2_agg$latitude, longitude=a1a2_agg$longitude))
          PA_matrix$year_bin <- year_bin
          PA_final <- rbind(PA_final, PA_matrix)
          biom_matrix <- data.frame(cbind(biom_matrix, latitude=a1a2_agg$latitude, longitude=a1a2_agg$longitude))
          biom_matrix$year_bin <- year_bin
          biom_final <- rbind(biom_final, biom_matrix)
        }

        filename <- paste(projfolder, species, '_rcp', rcp, '_', modelGCM, '_', nicheMod, '_iter', k,'.RData', sep='')
        # filename <- paste('output/CEmodels_Proj_Uncertainty_2018/', species, '_rcp', rcp, '_', modelGCM, '_', nicheMod, '_iter', k,'.RData', sep='')
        save(PA_final, biom_final, file=filename)
           