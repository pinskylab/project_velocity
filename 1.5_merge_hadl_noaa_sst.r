
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
 
#######################################################
#  Surftemp observations are missing almost entirely from three surveys (Newfoundland Fall and Spring, WC_Ann), as well as missing from individuals hauls in the other surveys. We used the global Merged Hadley-NOAA/OI Sea Surface Temperature and Sea-Ice Concentration dataset to fill in these missing observations when possible. We validated this method by comparing the merged Hadley-NOAA values to observed survey surface temperatures when possible, and used regressions to correct for any observed bias.
# Merged Hadley-NOAA/OI Sea Surface Temperature & Sea-Ice Concentration (Hurrell et al, 2008) - See more at: https://climatedataguide.ucar.edu/climate-data/merged-hadley-noaaoi-sea-surface-temperature-sea-ice-concentration-hurrell-et-al-2008#sthash.T8eqg21N.dpuf
# Citation: James W. Hurrell, James J. Hack, Dennis Shea, Julie M. Caron, and James Rosinski, 2008: A New Sea Surface Temperature and Sea Ice Boundary Dataset for the Community Atmosphere Model. J. Climate, 21, 5145–5153. doi: http://dx.doi.org/10.1175/2008JCLI2292.1
# Data obtained from here:

####################################################### 
# Validation/Calibration
# Compare existing surftemp records to Hadley-NOAA SST to determine if there is bias. Correct for this bias using an appropriate regression:
# Use WCAnn to validate WCTri survey (same region). 
# Since only 9 surftemp measurements exist from Newfoundland (fall only), this is not enough to calibrate a relationship. Thus, use Hadley-NOAA SST for all measurements in this region.
# For all other regions, use regression with sst*region to translate Hadley-NOAA SST into "surftemp".
########################################################
 

library(lattice)
library(ggplot2)

# Read in file with HadlSST values for all hauls.
# setwd("/Users/abigailporay/Documents/Work/proj_ranges")
allSST<-read.csv("data/allhauls_withHadlSST_150527.csv",row.names=1)

# Load dat to later merge with new SST.
load('data/trawl_allregionsforprojections_2015-02-02.RData') # load dat data.frame. 
# dat = rbind(dat, seusspr2, seussum2, seusfal2, seus.shelf2) # attach seus data when allSST has been updated?

# Correct some surftemp measurements before regression 
plot(allSST$surftemp,allSST$sst) #obvious issues with surftemp = 0. Many of these also miss bottomtemp, but still an issue. Maybe some of them are true zeros, but not many.
ggplot(allSST, aes(x=surftemp, y=sst)) + geom_point() + facet_wrap(~region)# Going by region allows more errors to be seen
# Code outliers in observed temps_based on previous plot_with NAs so they can be filled in with OISST
allSST$surftemp[allSST$region=="AFSC_GOA" & allSST$surftemp<2] <- NA # 1 random outlier
allSST$surftemp[allSST$region=="DFO_SoGulf" & allSST$surftemp<2.5] <- NA # 28 outliers, also seems way to cold for Septemper, even up there
allSST$surftemp[allSST$surftemp>38 &is.na(allSST$surftemp)==F]<-NA # Corrects a single outlier on Scotian Shelf
allSST$surftemp[allSST$region=="DFO_ScotianShelf" & allSST$surftemp<5 & allSST$sst > 11 & is.na(allSST$surftemp)==F] <- NA # 23 outliers, also seems way to cold for July, even up there
# Zero records of temp on NEUS all seem suspect for both seasons
allSST$surftemp[allSST$surftemp==0 & is.na(allSST$surftemp)==F & allSST$region=="NEFSC_NEUSFall"]<-NA
allSST$surftemp[allSST$surftemp==0 & is.na(allSST$surftemp)==F & allSST$region=="NEFSC_NEUSSpring"]<-NA
# allSST$surftemp[allSST$surftemp==0]<-NA # Jim muted this, as many of the surveys do realistically extend into zero territory

# See if any outliers in dat file
abc <- unique(data.frame(haulid=dat$haulid,surftemp=dat$surftemp,bottemp=dat$bottemp, region=dat$region, stringsAsFactors=F))
ggplot(abc, aes(x=surftemp, y=bottemp)) + geom_point() + facet_wrap(~region)

# Use complete cases for regression
allSSTc<-allSST[complete.cases(allSST[,c("sst","surftemp")]),]

########################################################
# Regress trawl surftemp on Hadl-NOAA SST
########################################################

#All regions
lmAll<-lm(surftemp~sst*region,data=allSSTc)
# Region specific:
lmAleut<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="AFSC_Aleutians",]) #r2=.42; slope=.662
lmEBS<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="AFSC_EBS",]) # .55; .900
lmGOA<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="AFSC_GOA",]) # .86; 1.024
lmWC<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="AFSC_WCTri",]) # .32; .885
lmScShel<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="DFO_ScotianShelf",]) # .89; 1.006
lmSoGulf<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="DFO_SoGulf",]) # .52; 1.209
lmNeusF<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="NEFSC_NEUSFall",]) # .89; .919
lmNeusS<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="NEFSC_NEUSSpring",]) # .64; .643
lmGmex<-lm(surftemp~sst,data=allSSTc[allSSTc$region=="SEFSC_GOMex",]) # .38; .902

# For some regions-seasons there is a tight relationship, but not all....
# Further, based on NEUS, there may be a seasonal effect on the relationship, better in fall

########################################################
# Use regressions to bias-correct Hadl-Noaa SST
########################################################

# Use regional/survey regressions to "correct" regional data. 
# Use lmWC to bias-correct NWFSC_WCAnn (and WCTri) SST values
allSST$adjSST[allSST$region %in% c("AFSC_WCTri", "NWFSC_WCAnn")] <- predict(lmWC,newdata=allSST[allSST$region %in% c("AFSC_WCTri", "NWFSC_WCAnn"),])
allSST$adjSST[allSST$region %in% c("AFSC_Aleutians")] <- predict(lmAleut,newdata=allSST[allSST$region %in% c("AFSC_Aleutians"),])
allSST$adjSST[allSST$region %in% c("AFSC_EBS")] <- predict(lmEBS,newdata=allSST[allSST$region %in% c("AFSC_EBS"),])
allSST$adjSST[allSST$region %in% c("AFSC_GOA")] <- predict(lmGOA,newdata=allSST[allSST$region %in% c("AFSC_GOA"),])
allSST$adjSST[allSST$region %in% c("DFO_ScotianShelf")] <- predict(lmScShel,newdata=allSST[allSST$region %in% c("DFO_ScotianShelf"),])
allSST$adjSST[allSST$region %in% c("DFO_SoGulf")] <- predict(lmSoGulf,newdata=allSST[allSST$region %in% c("DFO_SoGulf"),])
allSST$adjSST[allSST$region %in% c("NEFSC_NEUSFall")] <- predict(lmNeusF,newdata=allSST[allSST$region %in% c("NEFSC_NEUSFall"),])
allSST$adjSST[allSST$region %in% c("NEFSC_NEUSSpring")] <- predict(lmNeusS,newdata=allSST[allSST$region %in% c("NEFSC_NEUSSpring"),])
allSST$adjSST[allSST$region %in% c("SEFSC_GOMex")] <- predict(lmGmex,newdata=allSST[allSST$region %in% c("SEFSC_GOMex"),])

# allSST$adjSST[!allSST$region %in% c("AFSC_WCTri", "NWFSC_WCAnn", "DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")] <- predict(lmAll, newdata=allSST[!allSST$region %in% c("AFSC_WCTri", "NWFSC_WCAnn", "DFO_NewfoundlandFall", "DFO_NewfoundlandSpring"),])

# Use Hadley-NOAA SST directly for all Newfoundland hauls, and don't use any measured surftemp (n=7)
allSST$adjSST[allSST$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")] <- allSST$sst[allSST$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")]

########################################################
# Some simple diagnostics
########################################################
allSST$diff <- allSST$adjSST - allSST$surftemp  
ggplot(allSST, aes(x=diff)) + geom_histogram(binwidth=1) + facet_wrap(~region)
cor(allSST$surftemp,allSST$adjSST,use="c")
# [1] 0.9769

plot(allSST$surftemp,allSST$adjSST,col=rainbow(length(unique(allSST$region)))[allSST$region],xlab="Survey surftemp",ylab="Bias-corrected SST",pch=1,cex=0.5)
plot(allSST$surftemp,allSST$sst,col=rainbow(length(unique(allSST$region)))[allSST$region],xlab="Survey surftemp",ylab="Hadley model SST",pch=1,cex=0.5)

########################################################
# Create a new surftemp2 column to merge with dat file:
########################################################

# Default is:
allSST$surftemp2 <- allSST$surftemp
# If Newfoundland, use adjSST (HADL-NOAA)
allSST$surftemp2[allSST$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")] <- allSST$adjSST[allSST$region %in% c("DFO_NewfoundlandFall", "DFO_NewfoundlandSpring")]
# If WCAnn, use adjSST (HADL-NOAA)
allSST$surftemp2[allSST$region == "NWFSC_WCAnn"] <- allSST$adjSST[allSST$region == "NWFSC_WCAnn"]
# If otherwise missing, using adjSST (HADL-NOAA)
allSST$surftemp2[is.na(allSST$surftemp)==T] <- allSST$adjSST[is.na(allSST$surftemp)==T]

########################################################
# Merge allSST with dat file to create a new column
########################################################

#Haul IDs got tweaked in allSST - WCAnn were rounded. 

allhauls<-unique(dat$haulid) #114879
# Since order is preserved,
# replace haulid with allhauls - test below with full dataframes to be sure data are correctly matched
allSST$haulid<-allhauls

#change colnames in allSST (lat/lon were altered to merge with Hadl-NOAA SST)
colnames(allSST)[3:4]<-c("lon_new","lat_new")
colnames(allSST)[10:11]<-c("lat","lon")

#test<-merge(dat,allSST,by="haulid",sort=F,all=T)
#plot(test$surftemp.x,test$surftemp.y) # perfect. order was preserved.
#test$diffs <- test$surftemp.y - test$surftemp.x
#rm(test)

allSST2<-allSST[,c("haulid","surftemp2")]
dat2<-merge(dat,allSST2,by="haulid",sort=F,all=T)

########################################################
# Rename columns so that we don't have to re-write code.
########################################################

dat2$surftemp_orig<-dat2$surftemp
dat2$surftemp<-dat2$surftemp2
dat2<-dat2[,colnames(dat2)!="surftemp2"]
dat<-dat2
rm(dat2, abc, allSST2, allSSTc, allhauls)

#add seus data to dat, although there is still no correction of missing surftemp data in seus. Not too many missing values though.
load("Jim/SEUS_output_Sep23_2016.Rdata")
rm(seus.shelf2, seusfal2, seusspr2, seussum2)
dat.seus$surftemp_orig <- dat.seus$surftemp # for now, these two temperature columns for seus are identical
nm = c('haulid', 'region', 'year', 'yearsurv', 'month', 'stratum', 'lat', 'lon', 'depth', 'surftemp', 'bottemp', 'spp', 'wtcpue', 'surftemp_orig')
dat.seus = dat.seus[,nm]
dat <- rbind(dat, dat.seus)

########################################################
# Save new dat file with surftemp
########################################################

#save(dat,file=paste("data/trawl_allregionsforprojections_wSST_",Sys.Date(),".RData",sep=""))
save(dat,file=paste("Jim/trawl_allregionsforprojections_wSEUS_wSST_",Sys.Date(),".RData",sep=""))
