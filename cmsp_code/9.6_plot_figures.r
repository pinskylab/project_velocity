# Plot figures for paper

## Set working directory
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages (chron and ncdf4)
	}
# could add code for Lauren's working directory here

#####################
## Useful functions
#####################
require(maps)
require(mapdata) # for hires world map
#require(rworldmap)
#require(rworldxtra)
require(lattice)
require(data.table)
require(beanplot)
require(RColorBrewer)

convcol <- function(x, colfun){
	x <- x-min(x) # set min to 0
	x <- x/max(x) # set max to 1
	cs <- colfun(x)
	return(rgb(cs[,1], cs[,2], cs[,3], maxColorValue=255))
}

ramp2hex <- function(x){
	cs <- colfun(x)
	return(rgb(cs[,1], cs[,2], cs[,3], maxColorValue=255))
}

addtemps <- function(x, pos){
	usr <- par('usr')
	if(pos=='right'){
		x1 <- usr[1] + (usr[2]-usr[1])*.7
		x2 <- usr[1] + (usr[2]-usr[1])*.9
	}
	if(pos=='left'){
		x1 <- usr[1] + (usr[2]-usr[1])*0.1
		x2 <- usr[1] + (usr[2]-usr[1])*0.3
	}
	y1 <- usr[3] + (usr[4]-usr[3])*0.1
	text(x1,y1,round(min(x),1), col=ramp2hex(0), cex=0.5)
	text(x2,y1,round(max(x),1), col=ramp2hex(1), cex=0.5)
}

se <- function(x,na.rm=FALSE){ # standard error
	if(!na.rm){
		return(sd(x, na.rm=FALSE)/sqrt(length(x)))
	}
	if(na.rm){
		return(sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
	}
}

##################
## Study regions
##################
clim <- read.csv('data/climGrid.csv')

# set regions and expansion for dots
regs <- c('AFSC_EBS', 'AFSC_Aleutians', 'AFSC_GOA', 'NWFSC_WCAnn', 'SEFSC_GOMex', 'NEFSC_NEUSSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'DFO_NewfoundlandFall')
cexs <- c(0.3, 0.12, 0.15, 0.6, 0.5, 0.3, 0.35, 0.6, 0.2)
regsnice = c('Eastern Bering Sea', 'Aleutian Islands', 'Gulf of Alaska', 'West Coast U.S.', 'Gulf of Mexico', 'Northeast U.S.', 'Scotian Shelf', 'So. Gulf of St. Lawrence', 'Newfoundland')
ylabs = c('Latitude (°N)', '', '', '', 'Latitude (°N)', '', '', 'Latitude (°N)', '')
xlabs = c('', '', '', 'Longitude (°E)', '', '', 'Longitude (°E)', 'Longitude (°E)', 'Longitude (°E)')
ylims = list(ebs = c(54,62.5), al = c(50, 57), goa = c(50, 65), wc = c(32.2, 48.5), gom = c(25.5,30.5), ne = c(33, 45), ss = c(41, 48), sl = c(45.5, 49.5), nf = c(42, 62))
xlims = list(ebs = c(-179.5,-154), al = c(169, 198), goa = c(-171, -132), wc = c(-126.5, -117), gom = c(-97.5,-86.5), ne = c(-77, -64.5), ss = c(-69, -56), sl = c(-66, -60), nf = c(-65, -43))
pos <- c('right', 'right', 'right', 'left', 'right', 'right', 'right', 'left', 'right')
bcol <- 'dark grey' # background color

# trim to these regions
clim <- droplevels(clim[clim$region %in% regs,])

# convert a version to lon to -180 to 180
clim2 <- clim
clim2$lon[clim2$lon>180] <- clim2$lon[clim2$lon>180] - 360

# get a map from rworldmap
#newmap <- getMap(resolution = "high")

# plot map
	colfun <- colorRamp(colors = c('blue', 'white', 'red'))
	
	# quartz(width=8.7/2.54,height=5/2.54)
	pdf(width=8.7/2.54,height=5/2.54, file=paste('cmsp_figures/study_regions_w_BT.pdf', sep=''))
	par(mai=c(0.15, 0.08, 0.1, 0.1), omi=c(0.15, 0.15, 0, 0), tck=-0.06, mgp=c(1.2,0.4,0), las=1, cex.main=0.5, cex.axis=0.5)
	layout(mat=matrix(c(1,2,3,5,4,6,7,5,8,9,10,11), byrow=TRUE, nrow=3))

	# Add continent-scale map
	with(clim[!is.na(clim$bottemp.clim.int),], plot(lon, lat, col=convcol(bottemp.clim.int, colfun), pch=15, cex=0.1, xlab='', ylab='', main='', xaxt='n', xlim=c(170,320)))
	axis(1, mgp=c(1.2, 0.02, 0), at=c(200,250,300), labels=c(-160, -110, -60))
	map('world2Hires', add=TRUE, xlim=c(170,320), col=bcol, lwd=0.2)
	addtemps(clim$bottemp.clim.int[!is.na(clim$bottemp.clim.int)], 'left')

	for(i in 1:length(regs)){
		inds <- clim2$region==regs[i] & !is.na(clim2$bottemp.clim.int)
		if(regs[i] =='AFSC_Aleutians'){
			with(clim[inds,], plot(lon, lat, col=convcol(bottemp.clim.int, colfun), pch=15, cex=cexs[i], xlab='', ylab='', xlim=xlims[[i]], ylim=ylims[[i]], main=regsnice[i], xaxt='n', asp=1))
			axis(1, mgp=c(1.2, 0.02, 0), at=seq(170,195,by=5), labels=c('170', '', '180', '', '-170', ''))
			mymap <- map('world2Hires', plot=FALSE, xlim=xlims[[i]])
			mymap$x[mymap$x > xlims[[i]][2] | mymap$x < xlims[[i]][1]] <- NA
			polygon(mymap, col=bcol, border=NA)
			addtemps(clim$bottemp.clim.int[inds], pos[i])
		}else{
			with(clim2[inds,], plot(lon, lat, col=convcol(bottemp.clim.int, colfun), pch=15, cex=cexs[i], xlab='', ylab='', xlim=xlims[[i]], ylim=ylims[[i]], main=regsnice[i], xaxt='n', asp=1))
			axis(1, mgp=c(1.2, 0.02, 0))
			map('worldHires',add=TRUE, col=bcol, fill=TRUE, border=FALSE)
			addtemps(clim$bottemp.clim.int[inds], pos[i])
		}
	}
	mtext('Longitude (°E)', side=1, outer=TRUE, cex=0.5)
	mtext('Latitude (°N)', side=2, outer=TRUE, cex=0.5, las=0, line=0.3)

	
	dev.off()
	
##########################
## MPA gains and losses
##########################
wdpaturnbynetbymod <- read.csv(paste('data/wdpaturnbynetbymod_fitallreg_xreg_rcp85&45.csv', sep=''), row.names=1) # network results
wdpaturnbyMPAbymod <- read.csv(paste('data/wdpaturnbyMPAbymod_fitallreg_xreg_rcp85&45.csv', sep=''), row.names=1) # individual MPA results

regs <- c('AFSC_EBS', 'AFSC_Aleutians', 'AFSC_GOA', 'NWFSC_WCAnn', 'SEFSC_GOMex', 'NEFSC_NEUSSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'DFO_NewfoundlandFall')
#regsnice = c('Eastern Bering Sea', 'Aleutian Islands', 'Gulf of Alaska', 'West Coast U.S.', 'Gulf of Mexico', 'Northeast U.S.', 'Scotian Shelf', 'So. Gulf of St. Lawrence', 'Newfoundland')
regsnice = c('EBS', 'AI', 'GoA', 'WC', 'GoM', 'Neast', 'SS', 'SGoSL', 'Newf')
cols <- list(c('#67a9cf','white','black','#2166ac'), c('#ef8a62', 'white','black','#b2182b')) # colors from Colorbrewer2 7-class RdBu
# trim data to focal regions
wdpaturnbyMPAbymod <- droplevels(wdpaturnbyMPAbymod[wdpaturnbyMPAbymod$region %in% regs,])

# re-order the regions
wdpaturnbyMPAbymod$region <- factor(wdpaturnbyMPAbymod$region, levels=regs)

# calculate means within models/rcps (across regions)
means <- with(wdpaturnbyMPAbymod, aggregate(list(flost=flost, fgainedalt=fgainedalt, beta_sor=beta_sor), by=list(rcp=rcp, model=model), FUN=mean, na.rm=TRUE)) # mean across MPAs within climate models

# calculate means within regions/models/rcps
means.reg <- with(wdpaturnbyMPAbymod, aggregate(list(flost=flost, fgainedalt=fgainedalt, beta_sor=beta_sor), by=list(rcp=rcp, model=model, region=region), FUN=mean, na.rm=TRUE)) # mean across MPAs within climate models

# calculate network means within models/rcps
means.net <- with(wdpaturnbynetbymod, aggregate(list(flost=flost, fgainedalt=fgainedalt, beta_sor=beta_sor), by=list(rcp=rcp, model=model), FUN=mean, na.rm=TRUE))

# calculate means within models/rcps for individual MPAs in the networks
means.net.indiv <- with(wdpaturnbyMPAbymod[!is.na(wdpaturnbyMPAbymod$network),], aggregate(list(flost=flost, fgainedalt=fgainedalt, beta_sor=beta_sor), by=list(rcp=rcp, model=model), FUN=mean, na.rm=TRUE)) # mean 

# combine network and individual MPA results
means.net$type <- 'Net'
means.net.indiv$type <- 'Ind'
means.net <- rbind(means.net, means.net.indiv)

	# means
	aggregate(list(flost=means.net$flost, fgained=means.net$fgainedalt, beta_sor=means.net$beta_sor), by=list(rcp=means.net$rcp, type=means.net$type), FUN=mean)
	aggregate(list(flost=means.net$flost, fgained=means.net$fgainedalt, beta_sor=means.net$beta_sor), by=list(rcp=means.net$rcp, type=means.net$type), FUN=se)

# plot of mean MPA change and network change
	# quartz(width=17.8/2.54,height=4/2.54)
	pdf(width=17.8/2.54,height=4/2.54, file=paste('cmsp_figures/MPA_turnover_by_region.pdf', sep=''))

	par(mfrow=c(1,4), las=2, mai=c(0.5,0.4,0.1, 0.05), omi=c(0,0,0,0), mgp=c(1.6,0.4,0), tcl=-0.2, cex.axis=0.8)

	beanplot(flost ~ rcp+region, data=means.reg, what=c(0,1,1,0), side='both', col=cols, border=NA, wd=0.18, names=regsnice, ylim=c(0,1), cut=0.01, ylab='Fraction lost') # flost
	mtext(side=3, text='a)', adj=0.05, line=-1, las=1, cex=0.7)
	abline(h=mean(means$flost[means$rcp==45]), col=cols[[1]][4], lty=3)
	abline(h=mean(means$flost[means$rcp==85]), col=cols[[2]][4], lty=3)

	beanplot(fgainedalt ~ rcp+region, data=means.reg, what=c(0,1,1,0), side='both', col=cols, border=NA, wd=0.18, names=regsnice, ylim=c(0,1), log='', cut=0.01, ylab='Fraction gained') # fgained
	mtext(side=3, text='b)', adj=0.05, line=-1, las=1, cex=0.7)
	abline(h=mean(means$fgainedalt[means$rcp==45]), col=cols[[1]][4], lty=3)
	abline(h=mean(means$fgainedalt[means$rcp==85]), col=cols[[2]][4], lty=3)

	beanplot(beta_sor ~ rcp+region, data=means.reg, what=c(0,1,1,0), side='both', col=cols, border=NA, wd=0.15, names=regsnice, ylim=c(0,1), cut=0.01, ylab='Similarity') #similarity
	mtext(side=3, text='c)', adj=0.05, line=-1, las=1, cex=0.7)
	abline(h=mean(means$beta_sor[means$rcp==45]), col=cols[[1]][4], lty=3)
	abline(h=mean(means$beta_sor[means$rcp==85]), col=cols[[2]][4], lty=3)


	beanplot(flost ~ rcp+type, data=means.net, what=c(0,1,1,0), at=c(1,2), side='both', col=cols, border=NA, wd=0.1, cut=0.01, ylab='Fraction', xlim=c(0.5,6.5), ylim=c(0,1)) #similarity
	beanplot(fgainedalt ~ rcp+type, data=means.net, what=c(0,1,1,0), at=c(3,4), side='both', col=cols, border=NA, wd=0.08, cut=0.01, add=TRUE) #similarity
	beanplot(beta_sor ~ rcp+type, data=means.net, what=c(0,1,1,0), at=c(5,6), side='both', col=cols, border=NA, wd=0.1, cut=0.01, add=TRUE) #similarity
	mtext(side=3, text='d)', adj=0.05, line=-1, las=1, cex=0.7)

	mtext(side=1, text='Lost', adj=0.12, line=2, las=1, cex=0.7)
	mtext(side=1, text='Gained', adj=0.5, line=2, las=1, cex=0.7)
	mtext(side=1, text='Sim.', adj=0.9, line=2, las=1, cex=0.7)
	
	dev.off()
	

####################################################################
# Compare the planning approaches against each climate projection
####################################################################
goalsmetbymod1out <- read.csv('output/goalsmetbymod_fitallreg_xreg_cmsphistonly_all.csv')
goalsmetbymod2out <- read.csv('output/goalsmetbymod_fitallreg_xreg_cmsp2per_all.csv')


myregs <- c('AFSC_EBS', 'AFSC_Aleutians', 'AFSC_GOA', 'NWFSC_WCAnn', 'SEFSC_GOMex', 'NEFSC_NEUSSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'DFO_NewfoundlandFall')
regnames = c('Eastern Bering Sea', 'Aleutian Islands', 'Gulf of Alaska', 'West Coast U.S.', 'Gulf of Mexico', 'Northeast U.S.', 'Scotian Shelf', 'So. Gulf of St. Lawrence', 'Newfoundland')
mods <- sort(unique(goalsmetbymod1out$model))
rcps <- sort(unique(goalsmetbymod1out$rcp))


	
# plot %goals met (solution #1 and #2)
	colmat <- t(col2rgb(brewer.pal(4, 'Paired')))
	cols <- rgb(red=colmat[,1], green=colmat[,2], blue=colmat[,3], alpha=c(255,255,255,255), maxColorValue=255)
	yaxts <- c('s', 'n', 'n', 's', 'n', 'n', 's', 'n', 'n')
	xaxts <- c('n', 'n', 'n', 'n', 'n', 'n', 's', 's', 's')
	outfile <- paste('cmsp_figures/MarZone_allregs_goalsmetbymod_forfigure.pdf')
		outfile
	ylims <- c(0.3,1)

	# quartz(width=8.7/2.54, height=8.7/2.54)
	pdf(width=8.7/2.54, height=8.7/2.54, file=outfile)
	par(mfrow=c(3,3), mai=c(0.05, 0.05, 0.2, 0.05), omi=c(0.4,0.4,0,0), cex.main=0.8, tcl=-0.3, mgp=c(2,0.7,0))

	for(i in 1:length(myregs)){ # for each region
		inds <- goalsmetbymod1out$model == mods[1] & goalsmetbymod1out$rcp == rcps[1] & goalsmetbymod1out$region==myregs[i]
		plot(goalsmetbymod1out$mid[inds], goalsmetbymod1out$pmet[inds], xlab='', ylab='', ylim=ylims, type='l', pch=16, las=1, col=cols[1], main=regnames[i], yaxt='n', xaxt='n')
		
		if(yaxts[i]=='s') axis(2, mgp=c(2,0.6,0))
		else axis(2, labels=FALSE)

		if(xaxts[i]=='s') axis(1, mgp=c(2,0.5,0))
		else axis(1, labels=FALSE)

		# plot histonly
		for(k in 1:length(mods)){
			for(j in 1:length(rcps)){
				if(!(k==1 & j==1)){
					inds <- goalsmetbymod1out$model == mods[k] & goalsmetbymod1out$rcp == rcps[j] & goalsmetbymod1out$region==myregs[i]
					points(goalsmetbymod1out$mid[inds], goalsmetbymod1out$pmet[inds], type='l', pch=16, col=cols[1])	
				}
			}
		}
		inds <- goalsmetbymod1out$region==myregs[i]
		ensmean <- aggregate(list(nmet=goalsmetbymod1out$nmet[inds], pmet=goalsmetbymod1out$pmet[inds]), by=list(mid=goalsmetbymod1out$mid[inds]), FUN=mean)
		lines(ensmean$mid, ensmean$pmet, col=cols[2], lwd=2)

		# plot 2per
		for(k in 1:length(mods)){
			for(j in 1:length(rcps)){
				inds <- goalsmetbymod2out$model == mods[k] & goalsmetbymod2out$rcp == rcps[j] & goalsmetbymod2out$region==myregs[i]
				points(goalsmetbymod2out$mid[inds], goalsmetbymod2out$pmet[inds], type='l', pch=16, col=cols[3])
			}	
		}
		inds <- goalsmetbymod2out$region==myregs[i]
		ensmean2 <- aggregate(list(nmet=goalsmetbymod2out$nmet[inds], pmet=goalsmetbymod2out$pmet[inds]), by=list(mid=goalsmetbymod2out$mid[inds]), FUN=mean)
		lines(ensmean2$mid, ensmean2$pmet, col=cols[4], lwd=2)
	}
	
	mtext(side=1,text='Year',line=1.6, outer=TRUE)
	mtext(side=2,text='Fraction goals met',line=1.8, outer=TRUE)
	
	dev.off()
	
###################################
## Compare costs for each plan
###################################

# Read in plans
marxfolder <- '../MarZone_runs/'
runname1s <- c('cmsphistonly_neus', 'cmsphistonly_ebs', 'cmsphistonly_newf', 'cmsphistonly_gmex', 'cmsphistonly_scot', 'cmsphistonly_sgulf', 'cmsphistonly_goa', 'cmsphistonly_ai', 'cmsphistonly_wc')
runname2s <- c('cmsp2per_neus', 'cmsp2per_ebs', 'cmsp2per_newf', 'cmsp2per_gmex', 'cmsp2per_scot', 'cmsp2per_sgulf', 'cmsp2per_goa', 'cmsp2per_ai', 'cmsp2per_wc')

outputfolder1s <- paste(marxfolder, runname1s, '_output', sep='')
outputfolder2s <- paste(marxfolder, runname2s, '_output', sep='')
consplan1s <- vector('list', length(runname1s))
for(i in 1:length(consplan1s)) consplan1s[[i]] <- read.csv(paste(outputfolder1s[i], '/', runname1s[i], '_best.csv', sep=''))
consplan2s <- vector('list', length(runname2s))
for(i in 1:length(consplan2s)) consplan2s[[i]] <- read.csv(paste(outputfolder2s[i], '/', runname2s[i], '_best.csv', sep=''))

# Make matrix of proportion in each zone in each region
mathist <- matrix(rep(NA,4*length(runname1s)), ncol=length(runname1s))
colnames(mathist) <- sapply(strsplit(runname1s, split='_'), '[[', 2) # note: [[ is to index into the list
rownames(mathist) <- 1:4
mat2per <- mathist
mathistraw <- mathist
mat2perraw <- mathist
for(i in 1:length(runname1s)){
	for(j in 1:4){
		mathist[j,i] <- sum(consplan1s[[i]]$zone==j)/nrow(consplan1s[[i]])
		mat2per[j,i] <- sum(consplan2s[[i]]$zone==j)/nrow(consplan2s[[i]])
		mathistraw[j,i] <- sum(consplan1s[[i]]$zone==j)
		mat2perraw[j,i] <- sum(consplan2s[[i]]$zone==j)
	}
}
	colSums(mathist)
	colSums(mat2per)
	colSums(mathistraw)
	colSums(mat2perraw)
	
l <- list(mathist=mathist, mat2per=mat2per)
mat <- do.call(cbind, l)[,order(sequence(sapply(l, ncol)))]

# change in zones
mat2perraw - mathistraw
round(mat2per - mathist,3)

# modify to trick barplot into letting me use 8 colors instead of 4
# see http://r.789695.n4.nabble.com/barplot-colors-td4662538.html
matmod <- cbind(c(mat[,1], rep(0,nrow(mat))), c(rep(0,nrow(mat)), mat[,2])) 
for(i in seq(3,ncol(mat),by=2)){
	matmod <- cbind(matmod, c(mat[,i], rep(0,nrow(mat))), c(rep(0,nrow(mat)), mat[,i+1])) 
}

# plot
#regnames = c('Eastern Bering Sea', 'Aleutian Is.', 'Gulf of Alaska', 'West Coast U.S.', 'Gulf of Mexico', 'Northeast U.S.', 'Scotian Shelf', 'So. Gulf of St. Lawrence', 'Newfoundland')
regnames = c('EBS', 'AI', 'GoA', 'WC', 'GoM', 'Neast', 'SS', 'SGoSL', 'Newf')
cols <- brewer.pal(8, 'RdYlBu')
cols <- cols[c(1:4, 8:5)]

# quartz(width=8.7/2.54, height=5/2.54)
pdf(width=8.7/2.54, height=5/2.54, file='cmsp_figures/MarZone_plans.pdf')
par(mai=c(0.6, 0.6, 0.2, 0.1), mgp=c(1.8, 0.4, 0), tcl=-0.3, las=1, cex.axis=0.8)
barplot(height=matmod, space=rep(c(1,0), ncol(mathist)), xaxt='n', col=cols, ylab='Proportion', xlim=c(1.5,37))
axis(1, at=seq(2,26,by=3), labels=regnames, cex.axis=0.8, las=2, mgp=c(1.8, 0.5, 0))
legend(x=28.5, y=1, fill=cols[4:1], legend=c('Energy', 'Fishing', 'Conservation','Free'), cex=0.5, bty='n')
legend(x=27, y=1, fill=cols[8:5], legend=rep("",4), cex=0.5, bty='n')

dev.off()