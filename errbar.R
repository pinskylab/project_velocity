#######################
# Errorbar function
errbar <- function (x, y, xl=0, xu=xl, yl=0, yu = yl, length = 0.08, plot='', xlim=c(NA,NA), ylim=c(NA,NA),...){
	# xl, yl, xu, and yu are distances
	xlim2 = c(min(x-xl, na.rm=TRUE), 1.1*max(x+xl, na.rm=TRUE))
	ylim2 = c(min(y-yl, na.rm=TRUE), 1.1*max(y+yl, na.rm=TRUE))
	if(xlim2[1]<xlim[1] | is.na(xlim[1])) xlim[1] <- xlim2[1]
	if(xlim2[2]>xlim[2] | is.na(xlim[2])) xlim[2] <- xlim2[2]
	#if(ylim2[1]<ylim[1] | is.na(ylim[1])) ylim[1] <- ylim2[1]
	#if(ylim2[2]>ylim[2] | is.na(ylim[2])) ylim[2] <- ylim2[2]
	if(is.na(ylim[1])) ylim[1] <- ylim2[1]
	if(is.na(ylim[2])) ylim[2] <- ylim2[2]

	if(plot=='plot'){
		plot(x,y, xlim=xlim, ylim=ylim, ...)	
	}
	if(plot=='lines'){
		lines(x,y, xlim=xlim, ylim=ylim, ...)	
	}
	if(!all(yl==0)){
    	arrows(x, y + yu, x, y - yl, angle = 90, code = 3, length = length, ...)
	}
	if(!all(xl==0)){
    	arrows(x+xu, y, x-xl, y, angle = 90, code = 3, length = length, ...)
	}
}