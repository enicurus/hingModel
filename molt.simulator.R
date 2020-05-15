###Simulate molt data

#Author: Ryan S. Terrill
#ornithoterrill at gmail dot com
#Last update: 17 September 2019


#n=numer of samples
#init = julian day of molt initation
#duration = molt duration
#var = SD of molt initiation
#start= julian day to start sampling
#end = julian day to end sampling


molt.simulate<-function(n,init,duration,var,start=1,end=365){
	num<-sample(start:end,n,replace=TRUE)
	out<-matrix(ncol=2,nrow=n)
	term=init+duration
	for(i in 1:n){
		initiation<-rnorm(n,mean=init,sd=var)
			if(num[i]<initiation[i]){
				out[i,]<-c(num[i],0)
			}
		else {
			out[i,]<-c(num[i],(((100/(term-init))*num[i]))-initiation[i]*(100/(term-init)))
			}
						termination<-list()
		termination[i]<-term-(init-initiation[i])	
		if(num[i]>termination[i]){
			out[i,]<-c(num[i],100)
	}
	}
	out<-data.frame(out);colnames(out)<-c("julian","pscore")
return(out)
}

