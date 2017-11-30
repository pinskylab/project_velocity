# seasonal mean sbt
  filename = paste('data/', pred.folder, modelrun[13], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[3], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMean))
  rownames(filein) <- NULL
  sbtIPSL <- filein
  tempsmeansbt <- as.vector(filein)

# min annual sbt
  filename = paste('data/', pred.folder, modelrun[13], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[2], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  sbtIPSLmin <- filein
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  tempsminsbt <- fileinVec

# max annual sbt
  filename = paste('data/', pred.folder, modelrun[13], '_rcp', rcp, '_fill.nc_2006_2100_', season,'_', pred.metric[1], '.txt', sep="")
  filein <- readLines(filename)
  filein <- t(sapply(filein, processlinesMinMax))
  rownames(filein) <- NULL
  sbtIPSLmax <- filein
  fileinVec <- as.vector(filein)
  fileinVec <- fileinVec[1:1281878]
  tempsmaxsbt <- fileinVec

  
  
  pred <- data.frame(cbind(clim.grid, year=rep(2007:2100, rep=94, each=13637), SBT.seasonal = tempsmeansbt, SBT.min = tempsminsbt, SBT.max = tempsmaxsbt))
  pred <- merge(pred, year.bin, by='year', all.x=T)  # Bin into 20 year periods to make some plots
  pred.bathE <- merge(pred, proj.grid[proj.grid$lonClimgrid > -105,], by=c('lonClimgrid', 'latClimgrid'), all.y = T, sort=F)
  pred.bathW <- merge(pred, proj.grid[proj.grid$lonClimgrid < -105,], by=c('lonClimgrid', 'latClimgrid'), all.y = T, sort=F)
  pred.bathE <- pred.bathE[order(pred.bathE$year, pred.bathE$index),] # reorder by year and index
  pred.bathW <- pred.bathW[order(pred.bathW$year, pred.bathW$index),] # reorder by year and index
  
  aggE <- aggregate(list(sbt.seasonal = pred.bathE$SBT.seasonal, SBT.min = pred.bathE$SBT.min,
                         SBT.max = pred.bathE$SBT.max), by=list(year_range=pred.bathE$bin, latitude=pred.bathE$latBathgrid, longitude=pred.bathE$lonBathgrid), FUN=mean)
  aggW <- aggregate(list(sbt.seasonal = pred.bathW$SBT.seasonal, SBT.min = pred.bathW$SBT.min,
                         SBT.max = pred.bathW$SBT.max), by=list(year_range=pred.bathW$bin, latitude=pred.bathW$latBathgrid, longitude=pred.bathW$lonBathgrid), FUN=mean)
  
  # 20yr. summary figures for each model
  pdf(width=14, height=5, file=paste('figures/Temp_projections/TEST', modelrun[13], '_rcp', rcp, '_20yr_tempProj.pdf', sep=''))
  print(levelplot(sbt.seasonal~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$sbt.seasonal-.1),max(aggE$sbt.seasonal+.1), length.out=40), 
                  col.regions=cols, layout = c(5, 1), panel = function(...) {
                    panel.fill(col = "light gray")
                    panel.levelplot(...)
                  }))  
  print(levelplot(SBT.min~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SBT.min_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.min-.1),max(aggE$SBT.min+.1), length.out=40), 
                  col.regions=cols, layout = c(5, 1), panel = function(...) {
                    panel.fill(col = "light gray")
                    panel.levelplot(...)
                  }))  
  print(levelplot(SBT.max~longitude*latitude|year_range, data=aggE, xlab=NULL, ylab=NULL, main=paste('SBT.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggE$SBT.max-.1),max(aggE$SBT.max+.1), length.out=40), 
                  col.regions=cols, layout = c(5, 1), panel = function(...) {
                    panel.fill(col = "light gray")
                    panel.levelplot(...)
                  }))  
  print(levelplot(sbt.seasonal~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('sbt.seasonal_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$sbt.seasonal-.1),max(aggW$sbt.seasonal+.1), length.out=40), 
                  col.regions=cols, layout = c(5, 1), panel = function(...) {
                    panel.fill(col = "light gray")
                    panel.levelplot(...)
                  }))  
  print(levelplot(SBT.min~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SBT.min_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.min-.1),max(aggW$SBT.min+.1), length.out=40), 
                  col.regions=cols, layout = c(5, 1), panel = function(...) {
                    panel.fill(col = "light gray")
                    panel.levelplot(...)
                  }))  
  print(levelplot(SBT.max~longitude*latitude|year_range, data=aggW, xlab=NULL, ylab=NULL, main=paste('SBT.max_rcp', rcp, '  season:',season, '  Model:', modelrun[i]), at = seq(min(aggW$SBT.max-.1),max(aggW$SBT.max+.1), length.out=40), 
                  col.regions=cols, layout = c(5, 1), panel = function(...) {
                    panel.fill(col = "light gray")
                    panel.levelplot(...)
                  }))  
dev.off()
