#Test how much variation/deviance is explained by temperature in abundance GAMs

library(mgcv)
source("CSquareCode.r") #taken from VMStools in googlecode
#start with dat file with selected species from fit_climate_envelopes.R

moddev<-data.frame()
sppocean = sort(unique(dat$sppocean))
for(i in 1:length(sppocean)){
sp<-sppocean[i]
mydat<-dat[dat$sppocean==sp,] 
myocean<-mydat[1,"ocean"]

#Need to add records for when there was a haul but no fish
zs<-dat[!dat$haulid %in% mydat$haulid & dat$ocean==myocean,] #extract haulids where missing
# Extract one record of each haulid
matchpos<-match(unique(zs$haulid),zs$haulid)
zeros<-zs[matchpos,] #create new dataframe with an entry for each haulid in that ocean where the species was not caught (need to update if we want regional models instead). Dummy data gets replaced next.

#Add/change relevant columns --> zero catch for target spp.
zeros$spp<-mydat$spp[1]
zeros$sppl<-mydat$sppl[1]
zeros$sppnew<-mydat$sppnew[1]
zeros$sppregion<-paste(zeros$sppnew,zeros$region,sep="_")
zeros$sppocean<-mydat$sppocean[1] #may need to add "ocean"
zeros$wtcpue<-0
zeros$presfit<-F
zeros$wtcpuena<-1e-4 #for consistency, but should be zero
zeros$wtcpuenal<-log(zeros$wtcpuena)

mydatw0<-rbind(mydat,zeros) #combine positive hauls with zero hauls

#mydatw0$cs30m<-as.factor(CSquare(mydatw0$lon,mydatw0$lat,1)) #1 degree squares

#Calculate mean catch per haul (by region??). For true biomass estimate, would need to stratify
ave.catch.wt<-tapply(mydatw0$wtcpue,list(mydatw0$year,mydatw0$region),mean,na.rm=T)
if(dim(ave.catch.wt)[2]<2) {
	avecatchyrreg<-cbind(as.data.frame(ave.catch.wt), rep(colnames(ave.catch.wt),dim(ave.catch.wt)[2]),rownames(ave.catch.wt))
	colnames(avecatchyrreg)<-c("biomassmean","region","year")
	}else{
#if(sum(ave.catch.wt>0,na.rm=T)<30) next
#avecatchyr<-cbind("year"=as.numeric(names(ave.catch.wt)),"biomassmean"=ave.catch.wt)
#avecatchyr<-cbind("year"=as.numeric(rownames(ave.catch.wt)),ave.catch.wt)
	avecatchyrreg<-cbind(stack(as.data.frame(ave.catch.wt)),rep(rownames(ave.catch.wt),dim(ave.catch.wt)[2]))
	colnames(avecatchyrreg)<-c("biomassmean","region","year")
}
spdata<-merge(mydatw0,avecatchyrreg)

spdata<-spdata[complete.cases(spdata[,c("surftemp","bottemp","rugosity")]),] #Creates problems for regions missing surf temp...

# mygam <- gam(presfit~s(bottemp, by=regionfact,k=4)+s(surftemp,by=regionfact,k=4) + s(logrugosity,by=regionfact,k=4) ,family="binomial",data=spdata) #+ s(biomassmean,k=4)
# preds<- predict(mygam,spdata,type="response")
# plot(spdata$bottemp,preds,col=spdata$regionfact)

### Abundance GAMs

#For Gulf of Mexico don't include factor for Region.
if(myocean=="Gulf"){
	
	gamha2 = try(gam(wtcpuenal ~ s(surftemp,k=4) + s(bottemp,k=4) +s(rugosity,k=4) + s(biomassmean,k=4), family=gaussian, data=spdata[spdata$presfit,]) ) # 
	gamha3 = try(gam(wtcpuenal ~ s(rugosity,k=4) + s(biomassmean,k=4), family=gaussian, data=spdata[spdata$presfit,]) ) # 
	gamha4 = try(gam(wtcpuenal ~ s(surftemp,k=4) + s(bottemp,k=4), family=gaussian, data=spdata[spdata$presfit,])) # 	

} else {

	gamha2 = try(gam(wtcpuenal ~ s(surftemp,k=4) + s(bottemp,k=4) +s(rugosity,k=4) + regionfact + s(biomassmean,k=4), family=gaussian, data=spdata[spdata$presfit,]) ) # 
	gamha3 = try(gam(wtcpuenal ~ s(rugosity,k=4) + regionfact + s(biomassmean,k=4), family=gaussian, data=spdata[spdata$presfit,]) ) # 
	gamha4 = try(gam(wtcpuenal ~ s(surftemp,k=4) + s(bottemp,k=4), family=gaussian, data=spdata[spdata$presfit,]))
}

moddev[i,1]<-sp
moddev[i,2]<-sum(spdata$presfit)
if(class(gamha2)[1] != "try-error"){
	moddev[i,3]<-round(summary(gamha2)$dev.expl,3) #Full model
	moddev[i,4]<-round(summary(gamha3)$dev.expl,3) #Full model excluding temp terms
	moddev[i,5]<-round(summary(gamha4)$dev.expl,3) #Temp terms only
	}
}
colnames(moddev)<-c("spp_ocean","n_pres","full_model","full_minus_temp","temp_only")
write.table(moddev,file="AbundanceGAM_Deviance_summary.txt",row.names=F,quote=F,sep=",")