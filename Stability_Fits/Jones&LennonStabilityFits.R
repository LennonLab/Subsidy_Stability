### Maximum likelihood fits for stability response variables
### Jones & Lennon 
rm(list=ls())

# Terrestrial DOC supply for ponds
DOCsupply=c(25.7,115.2,151.6,84.1,55.7,0,184.3,201.5,133.4,212.7) #gC/m2

###########################################################
#### FITTING SIGMOIDAL FUNCTION TO SRP RESILIENCE DATA ####
###########################################################
SRPrt=c(0,3.5,7.8,0,1.1,0,11.8,9.5,10,7.5)

### Zero doesn't work for a return time (forces divide by zero), replace with 1 day as this is where we sampled next following nutrient perturbation
SRPrt4res=SRPrt
SRPrt4res[SRPrt4res==0]=1
SRPres=1/SRPrt4res	# SRP resilience

dev.new()
plot(DOCsupply,SRPres,xlab="Load (gC m-2)",ylab="SRP Resilience",cex=1.5,pch=19,ylim=c(0,1.1))

###### Fitting sigmoidal to SRP resilience
# flexible logistic NLL function (offset and scaling included)
flexLogistNLL<-function(p,x,y){
	a=p[1]
	b=p[2]
	c=p[3]
	d=p[4]
	varx=exp(p[5])
	
	yhat=d*(exp(a+b*x)/(1+exp(a+b*x)))+c
	res=y-yhat
	nres=length(res)
	nll=(0.5*nres*log(2*pi*varx))+((res%*%res)/(2*varx))
	
	return(nll)
}

# initial parameter values for optimization
guess=c(12,-0.12,0.1,0.9,-2)
# find most likely parameters
logFIT=optim(guess,flexLogistNLL,x=DOCsupply,y=SRPres)

# plot model fit
xpred=seq(0,215)
lines(xpred,logFIT$par[4]*(exp(logFIT$par[1]+logFIT$par[2]*xpred)/(logFIT$par[4]+exp(logFIT$par[1]+logFIT$par[2]*xpred)))+logFIT$par[3],lwd=3,col='blue',lty=2)



##################################################################
#### FITTING GAUSSIAN FUNCTION TO FOOD WEB DISCPLACEMENT DATA ####
##################################################################
disp=read.table("Jones&LennonFoodwebDisplacement.txt",header=TRUE,row.names=1,sep="\t")

#disp[,3:4]=disp[,3:4]*32/1000

# Gaussian NLL function
gaussNLL<-function(p,x,y){
	a=p[1]
	b=p[2]
	c=p[3]
	varx=exp(p[4])
	
	yhat=a*exp(-(x-b)^2/(2*c^2))
	res=y-yhat
	nres=length(res)
	nll=(0.5*nres*log(2*pi*varx))+((res%*%res)/(2*varx))
	
	return(nll)
}

### Estimate parameters for the four food web response variables
CHLguess=c(100,75,75,2)
CHLfit=optim(CHLguess,gaussNLL,x=DOCsupply,y=disp[,1])

ZOOPguess=c(10,50,25,4)
ZOOPfit=optim(ZOOPguess,gaussNLL,x=DOCsupply,y=disp[,2])

GPPguess=c(0.6,75,25,3)
GPPfit=optim(GPPguess,gaussNLL,x=DOCsupply,y=disp[,3])

Rguess=c(0.5,60,25,3)
Rfit=optim(Rguess,gaussNLL,x=DOCsupply,y=disp[,4])

#### plot observations and model fits
pl=seq(0,225)

dev.new()
par(mfrow=c(2,2))

plot(DOCsupply,disp[,1],type='p',pch=16,cex=1.5,xlab='DOC supply (gC m-2)',ylab='Chlorophyll displacement (ug L-1)',ylim=c(0,100))
lines(pl,CHLfit$par[1]*exp(-(pl-CHLfit$par[2])^2/(2*CHLfit$par[3]^2)),lwd=2,lty=2)

plot(DOCsupply,disp[,2],type='p',pch=16,cex=1.5,xlab='DOC supply (gC m-2)',ylab='Zooplankton displacement (mg L-1)',ylim=c(0,10))
lines(pl,ZOOPfit$par[1]*exp(-(pl-ZOOPfit$par[2])^2/(2*ZOOPfit$par[3]^2)),lwd=2,lty=2)

plot(DOCsupply,disp[,3],type='p',pch=16,cex=1.5,xlab='DOC supply (gC m-2)',ylab='GPP displacement (mgO2 L-1 day-1)',ylim=c(-0.1,0.7))
lines(pl,GPPfit$par[1]*exp(-(pl-GPPfit$par[2])^2/(2*GPPfit$par[3]^2)),lwd=2,lty=2)

plot(DOCsupply,disp[,4],type='p',pch=16,cex=1.5,xlab='DOC supply (gC m-2)',ylab='R displacement (mgO2 L-1 day-1)',ylim=c(-0.1,0.7))
lines(pl,Rfit$par[1]*exp(-(pl-Rfit$par[2])^2/(2*Rfit$par[3]^2)),lwd=2,lty=2)

