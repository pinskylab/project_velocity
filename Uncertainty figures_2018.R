# This script produces one of the main figure types for the uncertainty paper. It mainly uses the 'uncert_final' dataframe 
  # That is produced from another script.
   
setwd('/Users/jamesmorley/Documents/project_velocity')
library(Hmisc)
library(data.table)
library(reshape2)
library(ggplot2)
 
species <- 'sebastes alutus_Pac'
projfolder <- paste('output/CEmodels_Proj_Uncertainty_2018', species, sep='_')

filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
load(filename)
 

# Centroid plots =============================================================================
 
uncert_cent <- uncert_final[,c(1:6,10)]
uncert_cent <- melt(uncert_cent, id.vars=c('nicheMod', 'rcp', 'modelrun', 'years', 'iter'), measure.vars=c('centroid_lat', 'centroidPA_lat'), variable.name='method', value.name='centroid')
# Make factor columns of everything...just to be sure everything works right
uncert_cent$nicheModel <- as.factor(paste(uncert_cent$nicheMod, uncert_cent$method, sep='_'))
uncert_cent$rcp_fact <- as.factor(uncert_cent$rcp)
uncert_cent$modelrun_fact <- as.factor(uncert_cent$modelrun)
uncert_cent$years_fact <- as.factor(uncert_cent$years) # for regression tree only

# Partition the variance with a custom function, from 'Dominance analysis.R'
mod_cent2020 <- dom_centroid('2007-2020', uncert_cent)
mod_cent2040 <- dom_centroid('2021-2040', uncert_cent)
mod_cent2060 <- dom_centroid('2041-2060', uncert_cent)
mod_cent2080 <- dom_centroid('2061-2080', uncert_cent)
mod_cent2100 <- dom_centroid('2081-2100', uncert_cent)

SStable <- data.table(rbind(mod_cent2020, mod_cent2040, mod_cent2060, mod_cent2080, mod_cent2100))
SStable_trans <- transpose(SStable)
colnames(SStable_trans) <- SStable$years
SStable_trans <- SStable_trans[2:5,]
rownames(SStable_trans) <- names(SStable)[2:5]
SStable_trans_centroid <- data.matrix(SStable_trans)
barplot(SStable_trans_centroid, main="", ylab='Variation (Sum of Squares)', col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x='topleft')) 
box()
# For BSB, the GLM_biomass model gave wonky results with a really low centroid estimate (30.06 lat) (all other factors combined)
    # also the GAM_biomass model was probably low, at 34.03 while the others were 36.6 (all other factors combined)

cent_mean <- aggregate(list(mean_cent = uncert_cent$centroid), by=list(years = uncert_cent$years), FUN=mean) 
cent_mean$xAxis <- c(1:5)

ggplot(uncert_cent, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=3, fill='white')


# Habitat North =====================================================================================
 
uncert_habN <- uncert_final[,c(1:5,9,13)]
uncert_habN <- melt(uncert_habN, id.vars=c('nicheMod', 'rcp', 'modelrun', 'years', 'iter'), measure.vars=c('habitat_north', 'habitat_northPA'), variable.name='method', value.name='habitat_north')
uncert_habN$nicheModel <- as.factor(paste(uncert_habN$nicheMod, uncert_habN$method, sep='_'))
uncert_habN$rcp_fact <- as.factor(uncert_habN$rcp)
uncert_habN$modelrun_fact <- as.factor(uncert_habN$modelrun)
uncert_habN$years_fact <- as.factor(uncert_habN$years) # for regression tree only

uncert_habN <- uncert_habN[,c(2,3,4,5,7,8)]
colnames(uncert_habN)[5] <- 'habitat'
niche <- c('GAM_habitat_north','GAM_habitat_northPA','GLM_habitat_north','GLM_habitat_northPA')
RCP <- c(26,45,85)
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')

hab_changeN <- data.frame(nicheMod=character(), rcp=character(), modelrun=character(), iter=integer(), time1=numeric(), time2=numeric(), time3=numeric(), time4=numeric(), stringsAsFactors = F)
hab_changN_fin <- NULL

for(i in 1:3){ # RCP
  rcp <- RCP[i]
  for(j in 1:18){ # GCMs
    gcm <- modelrun[j]
    for(k in 1:4){ #nicheMod
      nicModel <- niche[k]
      uncert_habN_sub <- uncert_habN[uncert_habN$rcp==rcp & uncert_habN$modelrun==gcm & uncert_habN$nicheModel==nicModel,]
      uncert_habN_cast <- dcast(uncert_habN_sub, rcp+modelrun+iter+nicheModel~years, value.var='habitat')  
      time1 <- ((uncert_habN_cast$'2021-2040' - uncert_habN_cast$'2007-2020')/uncert_habN_cast$'2007-2020')*100
      time2 <- ((uncert_habN_cast$'2041-2060' - uncert_habN_cast$'2007-2020')/uncert_habN_cast$'2007-2020')*100
      time3 <- ((uncert_habN_cast$'2061-2080' - uncert_habN_cast$'2007-2020')/uncert_habN_cast$'2007-2020')*100
      time4 <- ((uncert_habN_cast$'2081-2100' - uncert_habN_cast$'2007-2020')/uncert_habN_cast$'2007-2020')*100
      
      hab_changeN <- data.frame(nicheMod=nicModel, rcp=rcp, modelrun=gcm, iter=1:40, time1=time1, time2=time2, time3=time3, time4=time4, stringsAsFactors = F)
      hab_changN_fin <- rbind(hab_changN_fin, hab_changeN)
    }
  }
}

hab_changN_fin <- melt(hab_changN_fin, id.vars=c('nicheMod', 'rcp', 'modelrun', 'iter'), measure.vars=c('time1','time2','time3','time4'), variable.name='years', value.name='habitat')

hab_changN_fin$nicheModel <- as.factor(hab_changN_fin$nicheMod)
hab_changN_fin$rcp_fact <- as.factor(hab_changN_fin$rcp)
hab_changN_fin$modelrun_fact <- as.factor(hab_changN_fin$modelrun)
hab_changN_fin$years_fact <- as.factor(hab_changN_fin$years) # for regression tree only

# fit linear models
mod_cent2020 <- dom_percChan('time1', hab_changN_fin)
mod_cent2040 <- dom_percChan('time2', hab_changN_fin)
mod_cent2060 <- dom_percChan('time3', hab_changN_fin)
mod_cent2080 <- dom_percChan('time4', hab_changN_fin)

SStable <- data.table(rbind(mod_cent2020, mod_cent2040, mod_cent2060, mod_cent2080))
SStable_trans <- transpose(SStable)
colnames(SStable_trans) <- SStable$years
SStable_trans <- SStable_trans[2:5,]
rownames(SStable_trans) <- names(SStable)[2:5]
SStable_trans_habN <- data.matrix(SStable_trans)
barplot(SStable_trans_habN, main="", ylab='Variation (Sum of Squares)', col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x='topleft')) 
box()

# Violin plots
habN_mean <- aggregate(list(habitat = hab_changN_fin$habitat), by=list(years = hab_changN_fin$years), FUN=mean) 
habN_mean$xAxis <- c(1,2,3,4)

ggplot(hab_changN_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=habN_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changN_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habN_mean, aes(x=xAxis, y=habitat), shape=21, size=3, fill='white')



# Habitat South =====================================================================================

uncert_habS <- uncert_final[,c(1:5,8,12)]
uncert_habS <- melt(uncert_habS, id.vars=c('nicheMod', 'rcp', 'modelrun', 'years', 'iter'), measure.vars=c('habitat_south', 'habitat_southPA'), variable.name='method', value.name='habitat_south')
uncert_habS$nicheModel <- as.factor(paste(uncert_habS$nicheMod, uncert_habS$method, sep='_'))
uncert_habS$rcp_fact <- as.factor(uncert_habS$rcp)
uncert_habS$modelrun_fact <- as.factor(uncert_habS$modelrun)
uncert_habS$years_fact <- as.factor(uncert_habS$years) # for regression tree only

uncert_habS <- uncert_habS[,c(2,3,4,5,7,8)]
colnames(uncert_habS)[5] <- 'habitat'
niche <- c('GAM_habitat_south','GAM_habitat_southPA','GLM_habitat_south','GLM_habitat_southPA')
RCP <- c(26,45,85)
modelrun <- c('bcc-csm1-1-m','bcc-csm1-1','CanESM2','CCSM4','CESM1-CAM5','CNRM-CM5','GFDL-CM3','GFDL-ESM2M','GFDL-ESM2G','GISS-E2-R','GISS-E2-H','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM','MIROC5','MPI-ESM-LR','NorESM1-ME')

hab_changeS <- data.frame(nicheMod=character(), rcp=character(), modelrun=character(), iter=integer(), time1=numeric(), time2=numeric(), time3=numeric(), time4=numeric(), stringsAsFactors = F)
hab_changS_fin <- NULL

for(i in 1:3){ # RCP
  rcp <- RCP[i]
  for(j in 1:18){ # GCMs
    gcm <- modelrun[j]
    for(k in 1:4){ #nicheMod
      nicModel <- niche[k]
      uncert_habS_sub <- uncert_habS[uncert_habS$rcp==rcp & uncert_habS$modelrun==gcm & uncert_habS$nicheModel==nicModel,]
      uncert_habS_cast <- dcast(uncert_habS_sub, rcp+modelrun+iter+nicheModel~years, value.var='habitat')  
      time1 <- ((uncert_habS_cast$'2021-2040' - uncert_habS_cast$'2007-2020')/uncert_habS_cast$'2007-2020')*100
      time2 <- ((uncert_habS_cast$'2041-2060' - uncert_habS_cast$'2007-2020')/uncert_habS_cast$'2007-2020')*100
      time3 <- ((uncert_habS_cast$'2061-2080' - uncert_habS_cast$'2007-2020')/uncert_habS_cast$'2007-2020')*100
      time4 <- ((uncert_habS_cast$'2081-2100' - uncert_habS_cast$'2007-2020')/uncert_habS_cast$'2007-2020')*100
      
      hab_changeS <- data.frame(nicheMod=nicModel, rcp=rcp, modelrun=gcm, iter=1:40, time1=time1, time2=time2, time3=time3, time4=time4, stringsAsFactors = F)
      hab_changS_fin <- rbind(hab_changS_fin, hab_changeS)
    }
  }
}

hab_changS_fin <- melt(hab_changS_fin, id.vars=c('nicheMod', 'rcp', 'modelrun', 'iter'), measure.vars=c('time1','time2','time3','time4'), variable.name='years', value.name='habitat')

hab_changS_fin$nicheModel <- as.factor(hab_changS_fin$nicheMod)
hab_changS_fin$rcp_fact <- as.factor(hab_changS_fin$rcp)
hab_changS_fin$modelrun_fact <- as.factor(hab_changS_fin$modelrun)
hab_changS_fin$years_fact <- as.factor(hab_changS_fin$years) # for regression tree only

# fit linear models
mod_cent2020 <- dom_percChan('time1', hab_changS_fin)
mod_cent2040 <- dom_percChan('time2', hab_changS_fin)
mod_cent2060 <- dom_percChan('time3', hab_changS_fin)
mod_cent2080 <- dom_percChan('time4', hab_changS_fin)

SStable <- data.table(rbind(mod_cent2020, mod_cent2040, mod_cent2060, mod_cent2080))
SStable_trans <- transpose(SStable)
colnames(SStable_trans) <- SStable$years
SStable_trans <- SStable_trans[2:5,]
rownames(SStable_trans) <- names(SStable)[2:5]
SStable_trans_habS <- data.matrix(SStable_trans)
barplot(SStable_trans_habS, main="", ylab='Variation (Sum of Squares)', col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x='topleft')) 
box()

# Violin plots
habS_mean <- aggregate(list(habitat = hab_changS_fin$habitat), by=list(years = hab_changS_fin$years), FUN=mean) 
habS_mean$xAxis <- c(1,2,3,4)

ggplot(hab_changS_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=habS_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changS_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habS_mean, aes(x=xAxis, y=habitat), shape=21, size=3, fill='white')


# Plot up the manuscript layout
library(gridBase)
library(grid)
library(gridExtra)
 
  
# =================================================================================================================
# HALIBUT =========================================================================================================
# =================================================================================================================
filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
save(uncert_final, uncert_cent, SStable_trans_centroid, cent_mean, uncert_habN, hab_changN_fin, SStable_trans_habN, habN_mean, uncert_habS, hab_changS_fin, SStable_trans_habS, habS_mean, file=filename) 
  
centN <- ggplot(uncert_cent, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Centroid (latitude)') +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10), panel.border=element_blank(), axis.line=element_line(colour='black'))

habN <- ggplot(hab_changN_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habN_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changN_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habN_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

habS <- ggplot(hab_changS_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habS_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changS_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habS_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))
 
filename <- paste('figures/Uncertainty_MS_2018/', species, '_main_fig.pdf', sep='')
 
pdf(width=6.5, height=5.5, file=filename)
par(mfrow=c(2,3),mar=c(4,4.5,.2,.1), oma=c(.1,.1,.1,.1), tcl=-.25, cex.axis=1, cex.lab=1.3, mgp=c(2.4,.5,0))
barplot(SStable_trans_centroid, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"100)"), col=c("dodgerblue2","orange","green3","red"))
axis(2,at=c(0,200,400,600,800,1000,1200),labels=c(0,2,4,6,8,10,12), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habN, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^5~')'), col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x=3.6, y=1510000, cex=1.2, bty='n')) 
axis(2,at=c(0,200000,400000,600000,800000,1000000,1200000,1400000),labels=c(0,2,4,6,8,10,12,14), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habS, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^6~')'), col=c("dodgerblue2","orange","green3","red")) 
axis(2,at=c(0,500000,1000000,1500000,2000000,2500000,3000000),labels=c(0,0.5,1,1.5,2,2.5,3), lwd=1.3, cex.axis=1.1, las=1) #for halibut

vp <- viewport(height = unit(.59,"npc"), width=unit(0.33, "npc"), just = c("left","top"), y = .59, x = 0.015)
print(centN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.333, "npc"), just = c("left","top"), y = .59, x = 0.347)
print(habN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.34, "npc"), just = c("left","top"), y = .59, x = 0.668)
print(habS, vp = vp)

grid.text("(a)", x = unit(0.033, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(b)", x = unit(0.364, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(c)", x = unit(0.696, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(d)", x = unit(0.033, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(e)", x = unit(0.364, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(f)", x = unit(0.696, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
dev.off()

  
# =================================================================================================================
# SABLEFISH =========================================================================================================
# =================================================================================================================
filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
save(uncert_final, uncert_cent, SStable_trans_centroid, cent_mean, uncert_habN, hab_changN_fin, SStable_trans_habN, habN_mean, uncert_habS, hab_changS_fin, SStable_trans_habS, habS_mean, file=filename) 

centN <- ggplot(uncert_cent, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Centroid (latitude)') +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10), panel.border=element_blank(), axis.line=element_line(colour='black'))

habN <- ggplot(hab_changN_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habN_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changN_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habN_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

habS <- ggplot(hab_changS_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habS_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changS_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habS_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

filename <- paste('figures/Uncertainty_MS_2018/', species, '_main_fig.pdf', sep='')

pdf(width=6.5, height=5.5, file=filename)
par(mfrow=c(2,3),mar=c(4,4.5,.2,.1), oma=c(.1,.1,.1,.1), tcl=-.25, cex.axis=1, cex.lab=1.3, mgp=c(2.4,.5,0))
barplot(SStable_trans_centroid, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"100)"), col=c("dodgerblue2","orange","green3","red"))
axis(2,at=c(0,500,1000,1500,2000,2500,3000,3500),labels=c(0,5,10,15,20,25,30,35), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habN, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^5~')'), col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x=3.6, y=2710000, cex=1.2, bty='n')) 
axis(2,at=c(0,500000,1000000,1500000,2000000,2500000),labels=c(0,5,10,15,20,25), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habS, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^6~')'), col=c("dodgerblue2","orange","green3","red")) 
axis(2,at=c(0,1000000,2000000,3000000,4000000,5000000),labels=c(0,1,2,3,4,5), lwd=1.3, cex.axis=1.1, las=1) #for halibut

vp <- viewport(height = unit(.59,"npc"), width=unit(0.33, "npc"), just = c("left","top"), y = .59, x = 0.015)
print(centN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.333, "npc"), just = c("left","top"), y = .59, x = 0.347)
print(habN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.34, "npc"), just = c("left","top"), y = .59, x = 0.668)
print(habS, vp = vp)

grid.text("(a)", x = unit(0.033, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(b)", x = unit(0.364, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(c)", x = unit(0.696, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(d)", x = unit(0.033, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(e)", x = unit(0.364, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(f)", x = unit(0.696, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
dev.off()


# =================================================================================================================
# MARKET SQUID =========================================================================================================
# =================================================================================================================
species <- 'doryteuthis opalescens_Pac'
filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
save(uncert_final, uncert_cent, SStable_trans_centroid, cent_mean, uncert_habN, hab_changN_fin, SStable_trans_habN, habN_mean, uncert_habS, hab_changS_fin, SStable_trans_habS, habS_mean, file=filename) 
#load(filename)

centN <- ggplot(uncert_cent, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Centroid (latitude)') +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10), panel.border=element_blank(), axis.line=element_line(colour='black'))

habN <- ggplot(hab_changN_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habN_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changN_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habN_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

habS <- ggplot(hab_changS_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habS_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changS_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habS_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

filename <- paste('figures/Uncertainty_MS_2018/', species, '_main_fig.pdf', sep='')

pdf(width=6.5, height=5.5, file=filename)
par(mfrow=c(2,3),mar=c(4,4.5,.2,.1), oma=c(.1,.1,.1,.1), tcl=-.25, cex.axis=1, cex.lab=1.3, mgp=c(2.4,.5,0))
barplot(SStable_trans_centroid, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^4~')'), col=c("dodgerblue2","orange","green3","red"))
axis(2,at=c(0,50000,100000,150000,200000),labels=c(0,5,10,15,20), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habN, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^9~')'), col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x=3.6, y=12000000000, cex=1.2, bty='n')) 
axis(2,at=c(0,4000000000,8000000000,12000000000),labels=c(0,4,8,12), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habS, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^7~')'), col=c("dodgerblue2","orange","green3","red")) 
axis(2,at=c(0,5000000,10000000,15000000,20000000,25000000),labels=c(0,0.5,1,1.5,2,2.5), lwd=1.3, cex.axis=1.1, las=1) #for halibut
# 
vp <- viewport(height = unit(.59,"npc"), width=unit(0.33, "npc"), just = c("left","top"), y = .59, x = 0.015)
print(centN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.35, "npc"), just = c("left","top"), y = .59, x = 0.328)
print(habN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.34, "npc"), just = c("left","top"), y = .59, x = 0.668)
print(habS, vp = vp)

grid.text("(a)", x = unit(0.033, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(b)", x = unit(0.345, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(c)", x = unit(0.696, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(d)", x = unit(0.033, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(e)", x = unit(0.345, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(f)", x = unit(0.696, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
dev.off()


# =================================================================================================================
# PACIFIC OCEAN PERCH =========================================================================================================
# =================================================================================================================
species <- 'sebastes alutus_Pac'
filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
save(uncert_final, uncert_cent, SStable_trans_centroid, cent_mean, uncert_habN, hab_changN_fin, SStable_trans_habN, habN_mean, uncert_habS, hab_changS_fin, SStable_trans_habS, habS_mean, file=filename) 
#load(filename)

centN <- ggplot(uncert_cent, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Centroid (latitude)') +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10), panel.border=element_blank(), axis.line=element_line(colour='black'))

habN <- ggplot(hab_changN_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habN_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changN_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habN_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

habS <- ggplot(hab_changS_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habS_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changS_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habS_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

filename <- paste('figures/Uncertainty_MS_2018/', species, '_main_fig.pdf', sep='')

pdf(width=6.5, height=5.5, file=filename)
par(mfrow=c(2,3),mar=c(4,4.5,.2,.1), oma=c(.1,.1,.1,.1), tcl=-.25, cex.axis=1, cex.lab=1.3, mgp=c(2.4,.5,0))
barplot(SStable_trans_centroid, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^3~')'), col=c("dodgerblue2","orange","green3","red"))
axis(2,at=c(0,2000,4000,6000,8000),labels=c(0,2,4,6,8), lwd=1.3, cex.axis=1.1, las=1) 
barplot(SStable_trans_habN, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^5~')'), col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x=3.6, y=2850000, cex=1.2, bty='n')) 
axis(2,at=c(0,500000,1000000,1500000,2000000,2500000),labels=c(0,5,10,15,20,25), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habS, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^6~')'), col=c("dodgerblue2","orange","green3","red")) 
axis(2,at=c(0,1000000,2000000,3000000,4000000),labels=c(0,1,2,3,4), lwd=1.3, cex.axis=1.1, las=1) #for halibut


vp <- viewport(height = unit(.59,"npc"), width=unit(0.33, "npc"), just = c("left","top"), y = .59, x = 0.015)
print(centN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.333, "npc"), just = c("left","top"), y = .59, x = 0.347)
print(habN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.34, "npc"), just = c("left","top"), y = .59, x = 0.668)
print(habS, vp = vp)

grid.text("(a)", x = unit(0.033, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(b)", x = unit(0.364, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(c)", x = unit(0.696, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(d)", x = unit(0.033, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(e)", x = unit(0.364, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(f)", x = unit(0.696, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
dev.off()

 

# =================================================================================================================
# LOBSTER =========================================================================================================
# =================================================================================================================
filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
save(uncert_final, uncert_cent, SStable_trans_centroid, cent_mean, uncert_habN, hab_changN_fin, SStable_trans_habN, habN_mean, uncert_habS, hab_changS_fin, SStable_trans_habS, habS_mean, file=filename) 
 
centN <- ggplot(uncert_cent, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Centroid (latitude)') +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10), panel.border=element_blank(), axis.line=element_line(colour='black'))

habN <- ggplot(hab_changN_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habN_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changN_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habN_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

habS <- ggplot(hab_changS_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habS_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changS_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habS_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

filename <- paste('figures/Uncertainty_MS_2018/', species, '_main_fig.pdf', sep='')

pdf(width=6.5, height=5.5, file=filename)
par(mfrow=c(2,3),mar=c(4,4.5,.2,.1), oma=c(.1,.1,.1,.1), tcl=-.25, cex.axis=1, cex.lab=1.3, mgp=c(2.4,.5,0))
barplot(SStable_trans_centroid, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^3~')'), col=c("dodgerblue2","orange","green3","red"))
axis(2,at=c(0,2000,4000,6000,8000,10000,12000),labels=c(0,2,4,6,8,10,12), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habN, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^7~')'), col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x=3.65, y=27000000, cex=1.2, bty='n')) 
axis(2,at=c(0,5000000,10000000,15000000,20000000,25000000),labels=c(0,0.5,1,1.5,2,2.5), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habS, xaxt='n', yaxt="n",  main="", ylab=expression("Sum of squares ("%*%"1e"^6~')'), col=c("dodgerblue2","orange","green3","red")) 
axis(2,at=c(0,2000000,4000000,6000000,8000000),labels=c(0,2,4,6,8), lwd=1.3, cex.axis=1.1, las=1) #for halibut
# 

vp <- viewport(height = unit(.59,"npc"), width=unit(0.33, "npc"), just = c("left","top"), y = .59, x = 0.015)
print(centN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.34, "npc"), just = c("left","top"), y = .59, x = 0.336)
print(habN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.34, "npc"), just = c("left","top"), y = .59, x = 0.668)
print(habS, vp = vp)

grid.text("(a)", x = unit(0.033, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(b)", x = unit(0.364, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(c)", x = unit(0.696, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(d)", x = unit(0.033, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(e)", x = unit(0.364, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(f)", x = unit(0.696, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
dev.off()


# =================================================================================================================
# SUMMER FLOUNDER =========================================================================================================
# =================================================================================================================
species <- 'paralichthys dentatus_Atl'
filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
save(uncert_final, uncert_cent, SStable_trans_centroid, cent_mean, uncert_habN, hab_changN_fin, SStable_trans_habN, habN_mean, uncert_habS, hab_changS_fin, SStable_trans_habS, habS_mean, file=filename) 
load(filename)

centN <- ggplot(uncert_cent, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Centroid (latitude)') +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10), panel.border=element_blank(), axis.line=element_line(colour='black'))

habN <- ggplot(hab_changN_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habN_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changN_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habN_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

habS <- ggplot(hab_changS_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habS_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changS_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habS_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

filename <- paste('figures/Uncertainty_MS_2018/', species, '_main_fig.pdf', sep='')
 
pdf(width=6.5, height=5.5, file=filename)
par(mfrow=c(2,3),mar=c(4,4.5,.2,.1), oma=c(.1,.1,.1,.1), tcl=-.25, cex.axis=1, cex.lab=1.3, mgp=c(2.4,.5,0))
barplot(SStable_trans_centroid, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^3~')'), col=c("dodgerblue2","orange","green3","red"))
axis(2,at=c(0,2000,4000,6000,8000,10000,12000),labels=c(0,2,4,6,8,10,12), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habN, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^8~')'), col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x=3.65, y=175000000, cex=1.2, bty='n')) 
axis(2,at=c(0,50000000,100000000,150000000),labels=c(0,0.5,1,1.5), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habS, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^6~')'), col=c("dodgerblue2","orange","green3","red")) 
axis(2,at=c(0,1000000,2000000,3000000,4000000,5000000,6000000),labels=c(0,1,2,3,4,5,6), lwd=1.3, cex.axis=1.1, las=1) #for halibut
# yaxt="n", 

vp <- viewport(height = unit(.59,"npc"), width=unit(0.33, "npc"), just = c("left","top"), y = .59, x = 0.015)
print(centN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.333, "npc"), just = c("left","top"), y = .59, x = 0.34)
print(habN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.33, "npc"), just = c("left","top"), y = .59, x = 0.680)
print(habS, vp = vp)

grid.text("(a)", x = unit(0.033, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(b)", x = unit(0.357, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(c)", x = unit(0.696, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(d)", x = unit(0.033, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(e)", x = unit(0.357, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(f)", x = unit(0.696, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
dev.off()


# =================================================================================================================
# BLACK SEA BASS =========================================================================================================
# =================================================================================================================
species <- 'centropristis striata_Atl'
filename <- paste('output/Uncert_summary_', species, '_Nov2018.RData', sep='')
save(uncert_final, uncert_cent, SStable_trans_centroid, cent_mean, uncert_habN, hab_changN_fin, SStable_trans_habN, habN_mean, uncert_habS, hab_changS_fin, SStable_trans_habS, habS_mean, file=filename) 
load(filename)

centN <- ggplot(uncert_cent, aes(x=years, y=centroid)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_line(data=cent_mean, aes(x=xAxis, y=mean_cent)) + 
  geom_violin(data=uncert_cent, aes(x=years, y=centroid), fill='gray80') + geom_point(data=cent_mean, aes(x=xAxis, y=mean_cent), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Centroid (latitude)') +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10), panel.border=element_blank(), axis.line=element_line(colour='black'))

habN <- ggplot(hab_changN_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habN_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changN_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habN_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

habS <- ggplot(hab_changS_fin, aes(x=years_fact, y=habitat)) + geom_violin() + 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept=0, lwd=.5, color='gray', linetype='longdash') +
  geom_line(data=habS_mean, aes(x=xAxis, y=habitat)) + 
  geom_violin(data=hab_changS_fin, aes(x=years_fact, y=habitat), fill='gray80') + geom_point(data=habS_mean, aes(x=xAxis, y=habitat), shape=21, size=2, fill='white') +
  xlab('Years') + ylab('Habitat change (%)') +
  scale_x_discrete(labels=c('2021-2040','2041-2060','2061-2080','2081-2100')) +
  theme(axis.text.x = element_text(angle=290, hjust=0, vjust=1, colour='black'), axis.text.y = element_text(colour='black'), axis.title.y=element_text(size=10, margin=margin(l=0)), panel.border=element_blank(), axis.line=element_line(colour='black'))

filename <- paste('figures/Uncertainty_MS_2018/', species, '_main_fig.pdf', sep='')

pdf(width=6.5, height=5.5, file=filename)
par(mfrow=c(2,3),mar=c(4,4.5,.2,.1), oma=c(.1,.1,.1,.1), tcl=-.25, cex.axis=1, cex.lab=1.3, mgp=c(2.4,.5,0))
barplot(SStable_trans_centroid, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^4~')'), col=c("dodgerblue2","orange","green3","red"))
axis(2,at=c(0,20000,40000,60000,80000,100000,120000),labels=c(0,2,4,6,8,10,12), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habN, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^12~')'), col=c("dodgerblue2","orange","green3","red"), legend = c('Niche model','GCM','RCP','Parameter'), args.legend=list(x=3.6, y=2000000000000, cex=1.2, bty='n')) 
axis(2,at=c(0,500000000000,1000000000000,1500000000000,2000000000000),labels=c(0,0.5,1,1.5,2), lwd=1.3, cex.axis=1.1, las=1) #for halibut
barplot(SStable_trans_habS, xaxt='n', yaxt="n", main="", ylab=expression("Sum of squares ("%*%"1e"^9~')'), col=c("dodgerblue2","orange","green3","red")) 
axis(2,at=c(0,1000000000,2000000000,3000000000),labels=c(0,1,2,3), lwd=1.3, cex.axis=1.1, las=1) #for halibut
# 
vp <- viewport(height = unit(.59,"npc"), width=unit(0.33, "npc"), just = c("left","top"), y = .59, x = 0.015)
print(centN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.348, "npc"), just = c("left","top"), y = .59, x = 0.327)
print(habN, vp = vp)
vp <- viewport(height = unit(.59,"npc"), width=unit(0.34, "npc"), just = c("left","top"), y = .59, x = 0.668)
print(habS, vp = vp)
 
grid.text("(a)", x = unit(0.033, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(b)", x = unit(0.343, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(c)", x = unit(0.696, "npc"), y = unit(0.98, "npc"), gp=gpar(fontface='bold'))
grid.text("(d)", x = unit(0.033, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(e)", x = unit(0.343, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
grid.text("(f)", x = unit(0.696, "npc"), y = unit(0.573, "npc"), gp=gpar(fontface='bold'))
dev.off()

