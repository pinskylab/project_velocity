# Choose which species count as 'commercial' o

## Set working directories
if(Sys.info()["nodename"] == "pinsky-macbookair"){
	setwd('~/Documents/Rutgers/Range projections/proj_ranges/')
	marxfolder <- '../MarZone_runs/'
	}
if(Sys.info()["nodename"] == "amphiprion.deenr.rutgers.edu"){
	setwd('~/Documents/range_projections/')
	.libPaths(new='~/R/x86_64-redhat-linux-gnu-library/3.1/') # so that it can find my old packages
	marxfolder <- 'MarZone_runs/'
	}
# could add code for Lauren's working directory here

############
## Flags
############
runtype <- 'fitallreg'; projtype='_xreg'

# choose the rcp (for runs using just one)
rcp <- 85

# sort regions by ocean
pacregs <- c('AFSC_Aleutians', 'AFSC_EBS', 'AFSC_GOA', 'AFSC_WCTri', 'NWFSC_WCAnn')
atlregs <- c('DFO_NewfoundlandFall', 'DFO_NewfoundlandSpring', 'DFO_ScotianShelf', 'DFO_SoGulf', 'NEFSC_NEUSFall', 'NEFSC_NEUSSpring', 'SEFSC_GOMex')

#################
# Load functions
#################
## Install rfishbase if needed
#install.packages("rfishbase", repos = c("http://packages.ropensci.org", "http://cran.rstudio.com"), type="source")
#require('rfishbase')


###########################################################
# Determine which species are commercial or recreational
# from FishBase
###########################################################
#
## get list of projected species
#load(paste('data/taxthresh_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads biomassavemap data.frame
#	head(taxthresh)
#spps <- gsub('_Atl|_Pac', '', taxthresh$sppocean)
#
## validate names (a few minutes)
#spps_fish <- validate_names(spps) # will return warnings for any inverts
#spps_invert <- validate_names(spps, server='http://fishbase.ropensci.org/sealifebase')
#
## did we get them all?
#length(spps_fish)
#length(spps_invert)
#length(spps) == length(spps_fish) + length(spps_invert) # NO. lists are missing some species.
#
#species(spps_fish[1], fields=c('Importance')) # Importance seems very broad and not useful to me


###########################################################
# Determine which species are commercial or recreational
# from Sea Around Us
###########################################################
## get list of projected species
load(paste('data/presmap_', runtype, projtype, '_rcp', rcp, '.RData', sep='')) # loads biomassavemap data.frame
head(presmap)
spps <- sort(unique(as.character(presmap$sppocean)))
sppstrim <- gsub('_Atl|_Pac', '', spps)

# get tables of landings by LME from Sea Around Us
ebs <- read.csv('cmsp_data/SAU/SAU LME 1 v1-40.csv')
al <- read.csv('cmsp_data/SAU/SAU LME 65 v1-40.csv')
goa <- read.csv('cmsp_data/SAU/SAU LME 2 v1-40.csv')
wcann <- read.csv('cmsp_data/SAU/SAU EEZ 848 v1-40.csv')
wctri <- wcann # west coast annual and trienniel: same LME
gmex <- read.csv('cmsp_data/SAU/SAU LME 5 v1-40.csv')
neuss <- read.csv('cmsp_data/SAU/SAU LME 7 v1-40.csv')
neusf <- neuss # NEUS spring and fall: same LME
scot <- read.csv('cmsp_data/SAU/SAU LME 8 v1-40.csv')
sgulf <- scot # southern Gulf of St. Lawrence, same LME as Scotian Shelf 
newff <- read.csv('cmsp_data/SAU/SAU LME 9 v1-40.csv')
newfs <- newff # newfoundland spring and fall: same LME

# add region
ebs$region <- 'AFSC_EBS'
al$region <- 'AFSC_Aleutians'
goa$region <- 'AFSC_GOA'
wcann$region <- 'NWFSC_WCAnn'
wctri$region <- 'AFSC_WCTri'
gmex$region <- 'SEFSC_GOMex'
neuss$region <- 'NEFSC_NEUSSpring'
neusf$region <- 'NEFSC_NEUSFall'
scot$region <- 'DFO_ScotianShelf'
sgulf$region <- 'DFO_SoGulf'
newff$region <- 'DFO_NewfoundlandFall'
newfs$region <- 'DFO_NewfoundlandSpring'

# combine into a list
sau <- list(ebs=ebs, al=al, goa=goa, wcann=wcann, wctri=wctri, gmex=gmex, neuss=neuss, neusf=neusf, scot=scot, sgulf=sgulf, newff=newff, newfs=newfs)

# aggregate recent reported landings by species
sppston <- sau
for(i in 1:length(sau)){
	sppston[[i]] <- with(sau[[i]][sau[[i]]$reporting_status=='Reported' & sau[[i]]$catch_type=='Landings' & sau[[i]]$year>=1990,], aggregate(list(tonnes = tonnes), by=list(region=region, scientific_name = scientific_name, common_name=common_name), FUN=sum))
}

# turn SAU names to lower case to match projections
for(i in 1:length(sppston)){
	sppston[[i]]$scientific_name <- tolower(sppston[[i]]$scientific_name)
}

# set up vector hold names that match projections
for(i in 1:length(sppston)){
	sppston[[i]]$projname <- NA
}


# match names from SAU to names from projections
options(warn=1)
for(i in 1:length(sppston)){
	for(j in 1:length(sppston[[i]]$scientific_name)){
		ind <- agrep(sppston[[i]]$scientific_name[j], sppstrim) # find index into spps and sppstrim
		if(length(ind)==1) sppston[[i]]$projname[j] <- spps[ind] # only enter if exactly one match
	}
}


# order by decreasing landings
for(i in 1:length(sppston)){
	sppston[[i]] <- sppston[[i]][order(sppston[[i]]$tonnes, decreasing=TRUE),]
}

# manually fix a few that agrep missed or filled in mistakenly
for(i in 1:length(sppston)){
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'clupea pallasii pallasii' & sppston[[i]]$region %in% atlregs] <- 'clupea pallasii_Atl'
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'clupea pallasii pallasii' & sppston[[i]]$region %in% pacregs] <- 'clupea pallasii_Pac'
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'sebastes'] <- NA
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'doryteuthis pealeii' & sppston[[i]]$region %in% atlregs] <- 'loligo pealeii_Atl'
	sppston[[i]]$projname[sppston[[i]]$scientific_name == 'doryteuthis pealeii' & sppston[[i]]$region %in% pacregs] <- 'loligo pealeii_Pac'
}

# top 10 that are in projections for each region
sppston[[1]][1:(which(cumsum(!is.na(sppston[[1]]$projname))==10)[1]),]
sppston[[2]][1:(which(cumsum(!is.na(sppston[[2]]$projname))==10)[1]),]
sppston[[3]][1:(which(cumsum(!is.na(sppston[[3]]$projname))==10)[1]),]
sppston[[4]][1:(which(cumsum(!is.na(sppston[[4]]$projname))==10)[1]),]
sppston[[5]][1:(which(cumsum(!is.na(sppston[[5]]$projname))==10)[1]),]
sppston[[6]][1:(which(cumsum(!is.na(sppston[[6]]$projname))==10)[1]),]
sppston[[7]][1:(which(cumsum(!is.na(sppston[[7]]$projname))==10)[1]),]
sppston[[8]][1:(which(cumsum(!is.na(sppston[[8]]$projname))==10)[1]),]
sppston[[9]][1:(which(cumsum(!is.na(sppston[[9]]$projname))==10)[1]),]
sppston[[10]][1:(which(cumsum(!is.na(sppston[[10]]$projname))==10)[1]),]
sppston[[11]][1:(which(cumsum(!is.na(sppston[[11]]$projname))==10)[1]),]
sppston[[12]][1:(which(cumsum(!is.na(sppston[[12]]$projname))==10)[1]),]

# combine into table to write out
temp <- sppston[[1]][1:(which(cumsum(!is.na(sppston[[1]]$projname))==10)[1]),]
temp <- temp[!is.na(temp$projname),]
out <- temp
for(i in 2:length(sppston)){
	temp <- sppston[[i]][1:(which(cumsum(!is.na(sppston[[i]]$projname))==10)[1]),]
	temp <- temp[!is.na(temp$projname),]
	out <- rbind(out,temp)
}

# write out
write.csv(out, file='cmsp_data/fishery_spps.csv')