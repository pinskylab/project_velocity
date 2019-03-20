# Below runs all projections for GLM models. GAMs take much longer and were run on the cluster

setwd('/Users/jamesmorley/Documents/project_velocity')
projfolder <- 'output/CEmodels_Uncertainty_2018/'
library(brglm)
   
RCP <- c(26,45,85)
season <- 'jas'
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')
nicheMod <- 'GLM'
years <- c("2007-2020", "2021-2040", "2041-2060", "2061-2080", "2081-2100")
species <- 'sebastes alutus_Pac'

# Load niche model
load(paste('output/CEmodels_Uncertainty_2018/Cluster_nicheMods_', nicheMod, '_', species, '_uncer2018.RData', sep=''))
      
for(j in 1:3){
  rcp = RCP[j]
  
  for(z in 1:18){
    modelGCM = modelrun[z]
    print(paste('Beginning predictions for rcp', rcp, '_', modelGCM, '_', species, '_', nicheMod, '_', Sys.time(), sep=''))
    
      if(grepl('_Atl', species)){
        load(paste('data/prediction_files_Nov2017/predictionEASTnoBias_rcp',rcp, '_', season, '_', modelGCM, '.RData', sep=''))
        predictSpp <- pred.bathE
        rm(pred.bathE)
      } else{
        load(paste('data/prediction_files_Nov2017/predictionWESTnoBias_rcp',rcp, '_', season, '_', modelGCM, '.RData', sep=''))
        predictSpp <- pred.bathW
        rm(pred.bathW)
      }
      predictSpp$regionfact = names(myregions)[myregions==max(myregions)][1]
        
        PA_final <- NULL
        biom_final <- NULL          
        orig_coefPA <- mod1$coefficients # may not be necessary to store
        orig_coefBIOM <- mod2$coefficients # may not be necessary to store

        for(i in 1:5){ # five time periods
          print(paste('Starting time period ', i))
          year_bin <- years[i]
          pred_dat <- predictSpp[predictSpp$bin==year_bin,]
          
          if(grepl('_Atl', species)){
            PA_matrix <- matrix(nrow=82592, ncol=40)
            biom_matrix <- matrix(nrow=82592, ncol=40)
          }else{
            PA_matrix <- matrix(nrow=65826, ncol=40)
            biom_matrix <- matrix(nrow=65826, ncol=40)
          }
          
          for(k in 1:40){# 40 iterations
            print(paste('iteration ', k))
            mod1$coefficients <- br.1[k,]
            mod2$coefficients <- br.2[k,]
            
            pred.a1 <- predict(mod1, newdata = pred_dat, type='response')
            pred.a2 <- predict(mod2, newdata = pred_dat, type='response')
            pred.a2 <- pred.a1 * exp(pred.a2) 
            pred.a1a2 <- data.frame(cbind(latBathgrid=pred_dat$latBathgrid, lonBathgrid=pred_dat$lonBathgrid, PA=as.vector(pred.a1), biom=as.vector(pred.a2)))
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
        filename <- paste('output/CEmodels_Proj_Uncertainty_2018_', species, '/', species, '_rcp', rcp, '_', modelGCM, '_', nicheMod, '_iter', k,'.RData', sep='')
        save(PA_final, biom_final, file=filename)
  }
}
      
      # CHECK OUT SOME GRAPHS TO MAKE SURE LOOKS RIGHT
      #xlimit = c(min(PA_final$longitude),max(PA_final$longitude))
      #ylimit = c(min(PA_final$latitude), max(PA_final$latitude))
      #scale85 = seq(0, max(PA_final[,1:5]), length.out=20)
      #scale85_biom = seq(0, max(biom_final[,1:5]), length.out=20)
      #cols = colorRampPalette(colors = c('gray90', 'blue', 'dark blue', 'black'))
      #library(latticeExtra)
      #library(maps)
      
      #Projmap <- map('world', xlim=xlimit, ylim=ylimit, plot=F) 
      #Projmap <- data.frame(lon=Projmap$x, lat=Projmap$y)
      #filename <- paste('figures/Uncertainty_mods/PAiter_', species, '_rcp', rcp, '_', modelGCM, '_', nicheMod, '_iter', k,'.pdf', sep='')
      #pdf(height=5, width=10, file=filename)
      #print(levelplot(V10 ~ longitude*latitude|year_bin, data=PA_final, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
      #     at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
      #dev.off()
       
      #filename <- paste('figures/Uncertainty_mods/BIOMiter_', species, '_rcp', rcp, '_', modelGCM, '_', nicheMod, '_iter', k,'.pdf', sep='')
      #pdf(height=5, width=10, file=filename)
      #print(levelplot(V10 ~ longitude*latitude|year_bin, data=biom_final, xlim=xlimit, ylim=ylimit, xlab=NULL, ylab=NULL, 
      #                at = scale85, col.regions=cols, layout = c(5, 1)) + xyplot(lat ~ lon, Projmap, type='l', lty=1, lwd=1, col='dark gray')) 
      #dev.off()
      
      

