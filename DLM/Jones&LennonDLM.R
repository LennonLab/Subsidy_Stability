rm(list=ls())

#################################################################################
#### CUSTOM FUNCTIONS FOR ESTIMATING "BEST" DISCOUNT FACTORS AND RUNNING DLM ####
#################################################################################

###### runDLM with given discount factors 
runDLM<-function(dfw,dfv,dataIn){
	
	Y=dataIn
	
	output=matrix(0,11,length(Y))
	rownames(output)=c('mt','Ct','at','Rt','nt','dt','St','ft','Qt','et','At')
	
	### initial prior
	# level
	m0=mean(Y,na.rm=TRUE)
	C0=var(Y,na.rm=TRUE)
	at=m0
	Rt=C0/dfw

	# precision
	n0=2
	d0=2
	S0=d0/n0

	# missing observation
	if(is.na(Y[1])){
		mt=m0
		at=mt
		ft=at
		Rt=C0/dfw
		Qt=Rt+S0
		At=Rt/Qt
		Ct=C0/dfw
		#dt=d0
		#nt=n0
		#St=S0
		dt=d0*dfv
		nt=n0*dfv
		St=dt/nt

		# store model prediction
		output[,1]=c(mt,Ct,at,Rt,nt,dt,St,ft,Qt,NA,At)
	# observation exists
	}else{
		### initial forecast
		ft=at
		Qt=Rt+S0

		### Posterior or Updating
		et=Y[1]-ft
		At=Rt/Qt

		nt=dfv*n0+1
		dt=dfv*d0+S0*et*et/Qt
		St=dt/nt

		mt=at+At*et
		Ct=(St/S0)*(Rt-At*At*Qt)
		
		# store model prediction
		output[,1]=c(mt,Ct,at,Rt,nt,dt,St,ft,Qt,et,At)
	}
	
	# iterate through time series
	for(i in 2:length(Y)){
		# missing observation
		if(is.na(Y[i])){
			mt=mt
			at=mt
			ft=at
			Rt=Ct/dfw
			Qt=Rt+St
			At=Rt/Qt
			Ct=Ct/dfw
			dt=dfv*dt		
			nt=dfv*nt		
			St=dt/nt		
			# store model prediction
			output[,i]=c(mt,Ct,at,Rt,nt,dt,St,ft,Qt,NA,At)
		# observation exists
		}else{
			# prior level
			m1=mt
			C1=Ct
			at=mt
			Rt=Ct/dfw
	
			# prior precision
			nt1=nt
			dt1=dt
			St1=St
	
			# Forecast
			ft=at
			Qt=Rt+St1
	
			# Update
			et=Y[i]-ft
			At=Rt/Qt
			nt=dfv*nt1+1
			dt=dfv*dt1+St1*et*et/Qt
			St=dt/nt
	
			mt=at+At*et
			Ct=(St/St1)*(Rt-At*At*Qt)
			
			# store model prediction
			output[,i]=c(mt,Ct,at,Rt,nt,dt,St,ft,Qt,et,At)
		}
	}
	
	return(output)
}


##### identifies best discount factors for a local level DLM by calculating cumulative predictive density (CPD)
# student T cumulative predictive density function
studentTcpd<-function(x,n,f,q){
	px=((gamma((n+1)/2)*n^(n/2))/(gamma(n/2)*(pi*q)^0.5))*(n+((x-f)^2)/q)^(-1*(n+1)/2)
	return(px)
}

#dlmFITpd<-function(p,dataIn){
dlmFITpd<-function(dataIn){
	# vector of discount factors to do grid search over
	dfw_s=c(seq(0.8,0.975,0.025),0.99,1)	#
	dfv_s=c(seq(0.8,0.975,0.025),0.99,1)	#
	
	Y=dataIn
	fits=matrix(0,length(dfw_s),length(dfv_s))	#
	for(a in 1:nrow(fits)){		
		for(b in 1:ncol(fits)){		
			dfw=dfw_s[a]		
			dfv=dfv_s[b]		
	
			output=matrix(0,11,length(Y))
			rownames(output)=c('mt','Ct','at','Rt','nt','dt','St','ft','Qt','et','At')
	
			### initial prior level
			m0=mean(Y,na.rm=TRUE)
			C0=var(Y,na.rm=TRUE)
			at=m0
			Rt=C0/dfw

			# precision
			n0=2
			d0=2
			S0=d0/n0
	
			# if a missing observation
			if(is.na(Y[1])){
				mt=m0
				at=mt
				ft=at
				Rt=C0/dfw
				Qt=Rt+S0
				At=Rt/Qt
				Ct=C0/dfw
				dt=d0
				nt=n0
				St=S0
				output[,1]=c(mt,Ct,at,Rt,nt,dt,St,ft,Qt,NA,At)
			# if have an observation
			}else{
				### initial forecast
				ft=at
				Qt=Rt+S0

				### Posterior or Updating
				et=Y[1]-ft
				At=Rt/Qt

				nt=dfv*n0+1
				dt=dfv*d0+S0*et*et/Qt
				St=dt/nt

				mt=at+At*et
				Ct=(St/S0)*(Rt-At*At*Qt)

				output[,1]=c(mt,Ct,at,Rt,nt,dt,St,ft,Qt,et,At)
			}
	
			# iterate through time series
			for(i in 2:length(Y)){
				# if observation missing
				if(is.na(Y[i])){
					mt=mt
					at=mt
					ft=at
					Rt=Ct/dfw
					Qt=Rt+St
					At=Rt/Qt
					Ct=Ct/dfw
					output[,i]=c(mt,Ct,at,Rt,nt,dt,St,ft,Qt,NA,At)
				# if observation exists
				}else{
					# prior level
					m1=mt
					C1=Ct
					at=mt
					Rt=Ct/dfw
	
					# prior precision
					nt1=nt
					dt1=dt
					St1=St
	
					# Forecast
					ft=at
					Qt=Rt+St1
	
					# Update
					et=Y[i]-ft
					At=Rt/Qt
					nt=dfv*nt1+1
					dt=dfv*dt1+St1*et*et/Qt
					St=dt/nt
	
					mt=at+At*et
					Ct=(St/St1)*(Rt-At*At*Qt)
					
					# fill output matrix
					output[,i]=c(mt,Ct,at,Rt,nt,dt,St,ft,Qt,et,At)
				}
			}
	
			# number of observations
			N=sum(!is.na(output[10,]))
			# modeled
			PDin=output[,!is.na(output[10,])]
			# observed
			PDobs=Y[!is.na(output[10,])]
	
			# estimate predictive density
			PD=studentTcpd(PDobs[1],PDin[5,1],PDin[8,1],PDin[9,1])
			for(i in 2:N){
				PD=PD*studentTcpd(PDobs[i],PDin[5,i],PDin[8,i],PDin[9,i])
			}
			
			fits[a,b]=log(PD,10)
		}
	}		
	
	# use discount factors that maximize predictive density
	out=c(dfw_s[rowSums(fits==max(fits))==1],dfv_s[colSums(fits==max(fits))==1])#
	return(out)
}

### probability density function for a student t distribution
studentT<-function(x){
px=((gamma((n+1)/2)*n^(n/2))/(gamma(n/2)*(pi*q)^0.5))*(n+((x-f)^2)/q)^(-1*(n+1)/2)
return(px)
}

### calculates difference between 0.9 and integration of student T with given endpoints
find90<-function(d){
	diff=abs(0.9-integrate(studentT,f-d*sqrt(q),f+d*sqrt(q))$value)
	return(diff)
}

### dlm confidence interval
OSAF_cl<-function(dlm){
	osaf_cl=matrix(0,2,ncol(dlm))
	rownames(osaf_cl)=c('lower','upper')
	
	for(i in 1:ncol(osaf_cl)){
		n<<-dlm[5,i]	# maybe should be from the previous time (t-1)
		q<<-dlm[9,i]
		f<<-dlm[8,i]
		
		mult=optimize(find90,lower=0.5,upper=10)$minimum
		
		osaf_cl[,i]=c(f-mult*sqrt(q),f+mult*sqrt(q))
	}
	return(osaf_cl)
}




##############################
#### FORECASTING FROM DLM ####
##############################
Kforecast<-function(k,mt,Ct,dfw,St){
	out=matrix(0,4,k)
	rownames(out)=c('at','Rt','ft','Qt')
	out[1,]=mt
	out[3,]=mt
	out[2,1]=Ct+(Ct*(1-dfw)/dfw)
	out[4,1]=out[2,1]+St
	
	for(i in 2:k){
		out[2,i]=out[2,(i-1)]+(out[2,(i-1)]*(1-dfw)/dfw)
		out[4,i]=out[2,i]+St
	}
	
	return(out)
}

#### finds 90% probability limits using find90
Kfore_cl<-function(fore,nt){	
	fore_cl=matrix(0,2,ncol(fore))
	rownames(fore_cl)=c('lower','upper')
	n=nt
	for(i in 1:ncol(fore)){
		q<<-fore[4,i]
		f<<-fore[3,i]
		mult=optimize(find90,lower=0.5,upper=10)$minimum
		
		fore_cl[1,i]=f-mult*sqrt(q)
		fore_cl[2,i]=f+mult*sqrt(q)
	}		
	return(fore_cl)
}


###################################################################################
#### IMPLEMENT DLM FUNCTIONS (USES SOLUBLE REACTIVE PHOSPHORUS DATA AS EXAMPLE ####
###################################################################################
data=read.table("Jones&LennonSRP.txt",header=TRUE,sep="\t")

pulse=239	# Nutrient perturbation occurred on doy=239

DOCsupply=c(25.7,115.2,151.6,84.1,55.7,0,184.3,201.5,133.4,212.7) #gC/m2
ponds=unique(data$Pond)


#Set up pdf device for graphical output
pdf(file="SRP_DLMfits.pdf",width=11,height=8.5)

SRPparam_summ=matrix(0,length(ponds),5)
rownames(SRPparam_summ)=ponds
colnames(SRPparam_summ)=c('DFw','DFv')

for(ii in 1:length(ponds)){
	# isolate a given ponds lm data
	cur=data[data[,1]==ponds[ii],]
	doy=cur[,2]

	plot(doy,cur[,3],type='o',col='darkgrey',ylab="Soluble Reactive Phosphorus (ug L-1)",xlab='Day of Year',lwd=2,ylim=c(0,75))
	lines(c(pulse,pulse),c(0,100),col='black',lty=3,lwd=2)
	text(170,70,paste("DOC supply = ",DOCsupply[ii],sep=""))
	
	# pre spike data
	doy_pre=doy[doy<=pulse]
	
	Y=cur[1:length(doy_pre),3]
	
	df_s=dlmFITpd(Y)
	dfw=df_s[1]
	dfv=df_s[2]

	dlm=runDLM(dfw,dfv,Y)
	conf=OSAF_cl(dlm)

	lines(doy_pre,dlm[8,],lwd=2,col='blue')
	lines(doy_pre,conf[1,],lwd=2,lty=2,col='blue')
	lines(doy_pre,conf[2,],lwd=2,lty=2,col='blue')

	pred=Kforecast((length(doy)-length(doy_pre)),dlm[1,ncol(dlm)],dlm[2,ncol(dlm)],dfw,dlm[7,ncol(dlm)])
	foreCL=Kfore_cl(pred,dlm[5,ncol(dlm)])

	lines(doy[doy>pulse],pred[3,],lwd=2,col='red')
	lines(doy[doy>pulse],foreCL[1,],lwd=2,col='red',lty=2)
	lines(doy[doy>pulse],foreCL[2,],lwd=2,col='red',lty=2)
		
	post_pulse_obs=cur[doy>pulse,5]
	SRPparam_summ[ii,1]=dfw
	SRPparam_summ[ii,2]=dfv
		
}

print(SRPparam_summ)

dev.off()
