# Code for the 'annotate' computer system at Rutgers  
   
#.libPaths("/local/home/jamesm/R/x86_64-pc-linux-gnu-library/3.3")
#options(repos='http://cran.rstudio.com/')
  
library(caret) # For optimizing BRT model features
library(pdp) # for partial dependency plots
library(e1071)
library(gbm)
# NOTE that muted columns are often for testing on my cpu
  
#setwd('Documents/project_velocity')
sp <- 'sebastes alutus_Pac'
#filename <- paste('output/CEmodels_Uncertainty_2018/RawData_', sp, '_uncertainty2018.RData', sep='')  
filename <- paste('species_data_uncert/RawData_', sp, '_uncertainty2018.RData', sep='')  
load(filename)

# optimize with caret package_just for full mods as it takes a long time
# Do 10-fold CV, repeated 3x
fitControl <- trainControl(method="repeatedcv", number=10, repeats=3)
paste(Sys.time(), '_Starting cpue mod')
start_time <- Sys.time() 

myBRT2_caret <- train(logwtcpue.pad ~ SBT.seasonal + SST.seasonal.mean + SBT.min + SBT.max + SST.max + rugosity + GRAINSIZE + regionfact,
                      data=spdataB,
                      method="gbm",
                      distribution = "gaussian",
                      bag.fraction = 0.6,
                      tuneGrid=expand.grid(interaction.depth=c(6,8,10), n.trees=c(50, seq(500,10000,by=500)), shrinkage=c(0.02, 0.01, 0.005), n.minobsinnode=10), 
                      trControl=fitControl, verbose=FALSE)

filename <- paste('output_caret/', sp, '_optimizeBRT_Jul2018.RData', sep='')
save(myBRT2_caret, file=filename)
end_time <- Sys.time() 

paste(Sys.time(), '_Starting P/A mod')

# Need to convert presfit to a factor for train() to work
spdata$presfit_fact <- factor(spdata$presfit)
start_time_PA <- Sys.time() 

myBRT1_caret <- train(presfit_fact ~ SBT.seasonal + SST.seasonal.mean + SBT.min + SBT.max + SST.max + rugosity + GRAINSIZE + regionfact,
                      data=spdata,
                      method="gbm", 
                      distribution = "bernoulli",
                      bag.fraction = 0.6,
                      tuneGrid=expand.grid(interaction.depth=c(6,8,10), n.trees=c(50, seq(500,10000,by=500)), shrinkage=c(0.02, 0.01, 0.005), n.minobsinnode=10), 
                      trControl=fitControl, verbose=FALSE)

end_time_PA <- Sys.time() 
save(myBRT2_caret, myBRT1_caret, file=filename)

 
