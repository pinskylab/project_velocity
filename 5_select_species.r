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

library(data.table)
# Select species for climate-envelope models
 
#################
### Load data ###	
#################
#load('data/trawl_allregionsforprojections_2015-02-02.RData') # load dat data.frame. Has all trawl observations from all regions. wtcpue has the standardized biomass estimates. They are standardized within regions, but not across regions.
load('data/trawl_allregionsforprojections_wSST_2015-06-02.RData') # load dat data.frame. Has all trawl observations from all regions. wtcpue has the standardized biomass estimates. They are standardized within regions, but not across regions.
load('trawl_allregionsforprojections_wSEUS_wSST_2016-09-26.RData') # Currently on laptop

source("CSquareCode.r") #taken from VMStools in googlecode

############################################
# Standardize species names across regions #
############################################
# see spptaxonomy_2015-02-09_plusManual.csv for a useful conversion table
Spptax<-read.csv("data/spptaxonomy_plusManual.csv", stringsAsFactors=F) #note: new column in CSV file 
#spptax<-apply(Spptax,2,tolower) # makes everything characters strings, and also conforms to lower case text
#dat$sppl<-tolower(dat$spp)
#datspp<-unique(dat$sppl)
#spptax<-as.data.frame(spptax)
#sum(datspp %in% spptax$taxon)# 814 of 5077 spp matched

# Make a list of haul info, in case any entire hauls get dropped later on when dropping rare species
hauls <- dat
hauls$spp <- hauls$wtcpue <- hauls$sppl <- NULL
hauls <- unique(hauls)
 
taxa <- data.table(read.csv("Jim/spp.key.csv", stringsAsFactors=F)) # Brings in taxonomy file from Ryan's trawlData package
dattax <- unique(dat$spp)
sum(dattax %in% taxa$ref) # 4531 of 5306 spp matched with Ryan Batt's taxonomy file
sum(dattax %in% Spptax$taxon) #828 of 5306 matched with Malin's taxonomy file

# Create vectors that show where taxons match up with rows in both Ryan's and Malin's taxonomy files
matchesRy <- as.numeric(dattax %in% taxa$ref)
matchesMa <- as.numeric(dattax %in% Spptax$taxon)

# Loop that standardizes species names and identifies what ones were adjusted
newnames<-character(length(dattax))
how <-numeric(length(dattax))
for (i in 1:length(dattax)){
  if(matchesMa[i]==1){
    newnames[i] <- Spptax$newname[Spptax$taxon==dattax[i]]
    how[i] <- 1
  } else if(matchesRy[i]==1){
    newnames[i] <- taxa$spp[taxa$ref==dattax[i]]
    how[i] <- 2
  } else { 
    newnames[i] <- dattax[i] # if there is no match for a taxon, then use original name
    how[i] <- 0
  }
}
 
names <- data.frame(spp=dattax, sppnew=newnames, how=how, stringsAsFactors=F)  
write.csv(names, file="taxons.csv")  # Fix the remaining errors, at the species level anyway, manually in Excel
# I fixed/checked around half of the unmatched entries; mostly GMex entries with no 'N'. Most of the unchanged ones will
# probably get dropped anyway
names2 <- read.csv("Jim/taxons_Sep2016.csv", stringsAsFactors=F)
names2$sppnew.lower <- tolower(names2$sppnew) # make everything lower case, to further ensure taxa match up
abc <- data.frame(spp=names2$spp, sppnew=names2$sppnew.lower, stringsAsFactors=F)
length(unique(abc$sppnew)) # down to 4168 unique spp, from 5306

#Now merge new names with old names in dat file
dat<-merge(dat, abc, all.x=T, by="spp", sort=F) #dat now contains "sppnew"
rm(abc, Spptax, taxa)

# I stumbled on this error with pink shrimp after the above
dat$sppnew[dat$sppnew == "penaeus duorarum"] <- "farfantepenaeus duorarum"
dat$sppocean[dat$sppocean == "penaeus duorarum_Atl"] <- "farfantepenaeus duorarum_Atl"


######################
# Add useful columns #
######################

# Have a biomass that never goes to zero (useful for fitting log-links with a stratum effect) 
#dat$wtcpuena = dat$wtcpue
#dat$wtcpuena[dat$wtcpuena == 0] = 1e-4
#dat$wtcpuenal = log(dat$wtcpuena)

# other useful columns
#dat$presfit = dat$wtcpue > 0 # indicator for where present
# Jim: this isn't right, as all rows should be 'present' here, and many rows have wtcpue as zero, probably rounding or too little weight to measure at sea

#dat$stratumfact = as.factor(dat$stratum)
#dat$yrfact = as.factor(dat$year)
#dat$regionfact<-as.factor(dat$region)
#dat$sppregion = paste(dat$sppnew, dat$region, sep='_')
dat$ocean[dat$region %in% c("AFSC_Aleutians", "AFSC_EBS", "AFSC_GOA", "AFSC_WCTri", "NWFSC_WCAnn")] <- "Pac"
dat$ocean[dat$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring", "DFO_ScotianShelf","DFO_SoGulf","NEFSC_NEUSFall", "NEFSC_NEUSSpring",
                            "SCDNR_SEUSFall", "SCDNR_SEUSSummer", "SCDNR_SEUSSpring", "SCDNR_SEUSReef", "SEFSC_GOMex")] <- "Atl"
#dat$ocean[dat$region == "SEFSC_GOMex"] <- "Gulf" #Or should Gulf of Mex group with Altantic?
dat$sppocean = paste(dat$sppnew, dat$ocean, sep='_') 

##############################################################################
## Pick the spp to model/drop for thermal envelopes
##############################################################################

# First tag and drop rare taxa in each of Atl and Pac
freq <- data.frame(table(dat$sppocean))
ggplot(freq, aes(x=Freq)) + geom_histogram(binwidth=50)
freq$sppocean <- as.character(freq$Var1)
freq$Var1 <- NULL
dat <- merge(dat, freq, all.x=T, by="sppocean", sort=F)
length(unique(dat$sppocean)) # 4489 species-ocean combinations
dat <- dat[dat$Freq > 250,] # 250 represents ~10% of total hauls performed in surveys with least total hauls --> SEUS seasonal trawl surveys
length(unique(dat$sppocean)) # Down to 922 species-ocean combinations

# Drop taxa that are not to the species level, i.e., consist of only one word, or have an 'spp.'
myspp <- unique(dat$sppnew) # 892 total species
drop <- myspp[grep("spp.", myspp)] # 87 taxa
myspp <- myspp[!(myspp %in% drop)] # Down to 805 species
drop <- myspp[!grepl(" ", myspp)] # 82 taxa
myspp <- myspp[!(myspp %in% drop)] # Down to 723 species
# Below from a manual check of spp
drop <- c("gastropod eggs", "rajiformes egg case", "trash species in catch", "purple striated anemone","crustacea shrimp",
          "artediellus ","etropus ","paguridae f.","penaeus ","sebastes melanostictus or sebastes aleutianus","urophycis ","symphurus ","zero catch")
myspp <- myspp[!(myspp %in% drop)] # Down to 709 species
myspp <- sort(myspp)
# Malin had narrowed it down to 551 species; 158 fewer, but SEUS has been added since then, which would account for ~half of the difference
# I'll leave the other, potentially uncommon, species in...the frequency of capture in total, and in a region, can be preditors in 
# The uncertainty analysis, and I can come back and remove spp if the analysis suggests a threshold of "common-ness" needed.

dat <- dat[dat$sppnew %in% myspp,] 
save(dat, hauls, file="Jim/dat_selectedspp_Sep2016.Rdata")

#############################################################################################
# This code not used when Jim used this script to trim spp; including the loop

#droplist<-c("sp.", "egg","anemone","unident", "unknown", "teuthida","liparidinae","bathylagus sp.","lampanyctus sp.","caridea", "carinariidae","antipatharia","annelida" ,"crustacea shrimp", "artediellus _Atl", "asteroidea s.c._Atl", "bivalvia c._Atl", "cephalopoda c._Atl", "decapoda o._Atl", "gastropoda o._Atl", "mytilidae f._Atl", "ophiuroidea s.c._Atl", "paguridae f._Atl", "paguroidea s.f._Atl", "penaeus _Atl", "polychaeta c._Atl", "porifera p._Atl", "pycnogonida s.p._Atl", "scyphozoa c._Atl", "sepiolodae f._Atl", "anthozoa c._Atl", "beryciformes (order)_Atl", "cirripedia s.c._Atl", "clypeasteroida o._Atl", "etropus _Atl", "euphausiacea o._Atl", "fistularia _Atl", "halichondria cf. sitiens_Pac", "holothuroidea c._Atl", "macrouriformes (order)_Atl", "melanostomiidae (stomiatidae)_Atl", "octopoda o._Atl", "pandalidae f._Atl", "paralepididae _Atl", "seapen (order)_Atl", "symphurus _Atl", "tunicata s.p._Atl", "urophycis _Atl", "melanostomiidae (stomiatidae)_Atl", "gorgonocephalidae,asteronychidae f._Atl", "trash species in catch_Atl")
#for (i in 1:length(droplist)) {drop<-c(drop,myspp[grep(droplist[i],myspp, fixed=TRUE)])}

# Identify taxa caught in at least 10 survey years and at least 300 total observations > 0 in a given survey (hopefully avoiding changes in species classification through time).
# Species are identified by species+ocean for later trimming of dat.
# Some surveys have zero hauls, but not so systematically: these are very small catches (<1kg)

myregions<-unique(dat$region)
myspp<-NULL
for(r in myregions){
#	surveyyrs<-names(table(dat$year[dat$region==r]))
	regdat<-dat[dat$region==r & dat$wtcpue > 0 & !is.na(dat$wtcpue),] # trim to focal region and rows with catch > 0
	yrocc<-table(regdat$year,regdat$sppocean) #table of number of occurrances each year
	sumyrs<-apply(yrocc,2,function(x) sum(x>0)) #identify years with more than zero catch of each taxon
	sumobs<-colSums(yrocc)
	min1<-colnames(yrocc)[sumyrs>=10 & sumobs >= 300]  #identify taxa with at least 10 years of catch and at least 300 total catch records
	myspp<-c(myspp,min1) #now with myspp plus ocean
}
#table(myspp) #shows which species are selected in multiple surveys
myspp<-unique(myspp) 
length(myspp) # 663 unique taxa_ocean (up from 621 if we require presence in all years of a survey)

# However, we'll need to go back to the full dat file to determine where zero hauls occur, unless we assume each haul had at least one of these species. Check this here:
#dat2<-dat[dat$sppocean %in% myspp,]

########################################################################################

 
# add rugosity
rugfile<-read.csv("data/trawl_latlons_rugosity_forMalin_2015_02_10.csv")
# Need to add SEUS to this; need to do this as separate file to match up SEUS trawls with rugosity values

seus <- dat[dat$region=="SCDNR_SEUSFall" | dat$region=="SCDNR_SEUSSpring" | dat$region=="SCDNR_SEUSSummer",]
seus <- data.frame(lat=seus$lat, lon=seus$lon)
seus <- unique(seus) # a file with all unique lat/lon combos for SEUS
# 35.229 & -75.592
# 28.756 & -81.44
load("Jim/Depth.RData") # from trawlData package, I think it is depth data from Global Relief Model...that rugosity can be constructed with

 




#===============
# Some potentially old-eraseable code below
#===============

# Make subset of rugosity file that overlaps with SEUS box; then round to same number of sigfigs.
rug.seus <- rugfile[rugfile$lat < 35.5 & rugfile$lat > 27 & rugfile$lon > -81.7,]
# Had to round to .01 decimals to get overlap; this seems justified as resolution of rugosity is at .017 decimal degrees
rug.seus$lat.round <- signif(rug.seus$lat, digits=3)
rug.seus$lon.round <- signif(rug.seus$lon, digits=3)
seus$lat.round <- signif(seus$lat, digits=3)
seus$lon.round <- signif(seus$lon, digits=3)
# Plot below shows that the rugosity data file has values overlapping with entire seus trawl survey, but lat/lons still don't match up after rounding
plot(lat.round~lon.round, data=rug.seus)
points(lat.round~lon.round, col="red", cex=.2, data=seus)

cde <- rug.seus; cde$lat <- cde$lon <- cde$depth <- NULL
cde <- unique(cde)
abc <- seus; abc$lat <- abc$lon <- NULL
abc <- unique(abc)
blah <- merge(abc, cde, all.x=T)
# Need to bin coordinates, at a similar scale to what rugosity data are, which is 1-minute, or 0.017 decimal degrees
# I'll try binning to .010 decimal degrees to see if everything matches up, this can be done just by rounding

#====================




dat<-merge(dat,rugfile) #lose 69 instances of lumpenus lampretaeformis b/c missing lat/lon

rm(rugfile)
dat$logrugosity<-log(dat$rugosity+0.01) #log-transformed rugosity gave better model fits in initial tests of rugosity covariate

dat$cs1<-as.factor(CSquare(dat$lon,dat$lat,1)) #classify into 1 degree squares
#dat$cs6m<-as.factor(CSquare(dat$lon,dat$lat,0.1)) #classify into 6 arcminute squares


save(dat,file="data/dat_selectedspp.Rdata")


