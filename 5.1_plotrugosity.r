## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	}
if(Sys.info()["user"] == "lauren"){
	setwd('~/backup/NatCap/proj_ranges/')
}

####################
## helper functions
####################
require(RColorBrewer)
require(maps)

# normalize to 0-1
norm01 <- function(x){
	mn <- min(x)
	mx <- max(x)
	return((x-mn)/(mx-mn))
}

###############
# Plot rugosity
###############
rugfile<-read.csv("data/trawl_latlons_rugosity_forMalin_2015_02_10.csv")
rugfile$lr <- log(rugfile$rugosity+0.01)

colfun <- colorRamp(rev(brewer.pal(11, 'Spectral')))
xlims <- c(1.1, 0.9) * range(rugfile$lon)
ylims <- c(0.9, 1.1) * range(rugfile$lat)
# quartz(width=7, height=5)
pdf(width=7, height=5, file='figures/rugosity_map.pdf')
par(mai=c(0.7,0.7,0.3, 0.1), las=1, mgp=c(2,1,0))
map(database='world', fill=TRUE, col='dark gray', border=FALSE, xlim=xlims, ylim=ylims, xlab='Longitude', ylab='Latitude')
map.axes()
title('Rugosity')
points(rugfile$lon, rugfile$lat, col=rgb(colfun(norm01(rugfile$lr)), maxColorValue=255), pch=16, cex=0.2)
legend('bottomleft', legend=round(seq(min(rugfile$lr), max(rugfile$lr), length.out=10),2), col=rgb(colfun(norm01(seq(min(rugfile$lr), max(rugfile$lr), length.out=10))), maxColorValue=255), pch=16, cex=0.8, title='Log rugosity', bty='n')

dev.off()
