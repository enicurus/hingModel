########################################################
###Simulate molt data to test UZ vs threshold models###
########################################################

library(ggplot2)
library(dplyr)
library(moult)
library(chngpt)
library(reshape)
library(gridExtra)
library(grid)
library(cowplot)
library(ggformula)

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

source("~/Dropbox/HingeModels/molt.simulator.R")

### cols can be input init, input duration, input var,model init, model duration, model var

out<-matrix(nrow=9,ncol=6)

colnames(out)<-c("input_init","input_duration","input_var","model_init","model_duration","model_var")
rownames(out)<-c("5","10","15","20","30","40","50","75","100","150","200","250")

##start with init=150, dur = 50, var = 10

samples<-c(5,10,15,20,30,40,50,75,100,150,200,250)


simData<-list()



molt.oneSim<-function(){
for (i in samples){
	simData[[as.character(i)]]<-molt.simulate(n=i,init=150,duration=50,var=10)
	}
for (i in samples){
	simData[[as.character(i)]][,3]<-as.character(i)
	}
	return(simData)
}


sims<-list()

system.time(
for(j in 1:100){
sims[[as.character(j)]]<-molt.oneSim()
}
)

for(i in 1:100){
	sims[[i]]<-bind_rows(sims[[i]])
}

simsAll<-bind_rows(sims)

write.csv(simsAll,"~/Dropbox/HingeModels/simulatedData.csv")


sims$V3<-factor(sims$V3,levels=c(5,10,15,20,30,40,50,75,100,150,200,250))

ggplot(simsAll,aes(x=julian,y=pscore))+geom_point(alpha=.01)+facet_wrap(~V3)+theme_bw()

#ggsave("~/Dropbox/HingeModels/Figures/simulated_molt_points.pdf")


###fitting UZ models to data



for(i in 1:100){
sims[[i]]$pscore<-sims[[i]]$pscore/100
}


#this is wrapped in tryCatch to keep the loop going after "moult" spits an error - because too few molting birds
#will break the function and the loop

fits<-list()
fit<-list()
fit$'5'<-"na"
fit$'10'<-"na"

for(j in 1:100){
for(i in unique(sims[[j]]$V3)){
	tryCatch({
fit[[i]]<-moult(pscore~julian,data=sims[[j]][sims[[j]]$V3==i,])	
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	fits[[j]]<-fit
	}
}



# no need to define a list for it, see below
#se<-list()
ci<-list()
ci_end<-list()
#se_end<-list()
#ses<-list()
cis<-list()
cis_end<-list()
#ses_end<-list()
for(j in 1:length(fits)){
    for(i in 1:12){
    					tryCatch({
        # define a fit variable and reuse it later
        fit=fits[[j]][[i]]
        # se is a vector and you can use it directly. 
        se=sqrt(diag(vcov(fit)))
        ci[[i]]<-cbind(coef(fit)-1.96*se, coef(fit)+1.96*se)
       
        comb=c(1,1,0)        # the following matrix operation is simply doing coef(fit)[1] + coef(fit)[2] , but it is convenient to do it in matrix operation because the variance can be computed with it as well
        est =  comb%*%coef(fit)
        se_end=sqrt(comb%*%vcov(fit)%*%comb)      
        # do you need as.numeric?
        ci_end[[i]]<-cbind(as.numeric(est-1.96*se_end), as.numeric(est+1.96*se_end))
             			        		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    #ses[[j]]<-se

    names(ci)<-names(fits[[1]])
    cis[[j]]<-ci
    names(ci_end)<-names(fits[[1]])
    cis_end[[j]]<-ci_end
}




UZ.mat<-list()

for (i in 1:100){
UZ.mat[[i]]<-matrix(nrow=12,ncol=8)
rownames(UZ.mat[[i]])<-unique(sims[[i]]$V3)
colnames(UZ.mat[[i]])<-c("duration_min","duration_max","start_min","start_max","SD_min","SD_max","term_min","term_max")
}


for(j in 1:100){
	for(i in 1:12){
				tryCatch({
    		 	UZ.mat[[j]][i,1]<-cis[[j]][[i]][,1][1]
      			UZ.mat[[j]][i,2]<-cis[[j]][[i]][,2][1]
     			UZ.mat[[j]][i,3]<-cis[[j]][[i]][,1][2]
     			UZ.mat[[j]][i,4]<-cis[[j]][[i]][,2][2]
     			UZ.mat[[j]][i,5]<-cis[[j]][[i]][,1][3]
     			UZ.mat[[j]][i,6]<-cis[[j]][[i]][,2][3]
     			UZ.mat[[j]][i,7]<-cis_end[[j]][[i]][,1]
     			UZ.mat[[j]][i,8]<-cis_end[[j]][[i]][,2]
     			        		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
					}
		UZ.mat[[j]]<-data.frame(UZ.mat[[j]])
		UZ.mat[[j]]$n<-rownames(UZ.mat[[j]])	
}




UZ.results<-bind_rows(UZ.mat)



pdf("~/Dropbox/HingeModels/UZ.meanSims.pdf")
UZ.results$n<-factor(UZ.results$n,levels=c(5,10,15,20,30,40,50,75,100,150,200,250))
ggplot(UZ.results,aes(x=startDate))+
geom_histogram(bins=50,col="gray",fill="gray")+
facet_wrap(.~UZ.results$n)+
geom_segment(aes(x=150,y=0,xend=150,yend=40),linetype="dashed")+
theme_bw()
dev.off()


write.table(UZ.results,"~//Dropbox/HingeModels/UZ.csv",sep=",")







#Now calculate how many UZ models got the answer within the CI



UZ.results$correct<-UZ.results$start_min<151&UZ.results$start_max>149
UZ.results$termCorrect<-UZ.results$term_min<201&UZ.results$term_max>199


#count up corrects by n

Corrects<-data.frame(matrix(nrow=length(unique(UZ.results$n)),ncol=6))
colnames(Corrects)<-c("Correct","Incorrect","termCorrect","termIncorrect","model","n")
Corrects$n<-c(5,10,15,20,30,40,50,75,100,150,200,250)
Corrects$model<-"UZ"

for (i in unique(UZ.results$n)){
		tryCatch({
	num<-as.numeric(i)
	Corrects$Correct[Corrects$n==i]<-table(UZ.results$correct[UZ.results$n==i])[2]
	Corrects$Incorrect[Corrects$n==i]<-table(UZ.results$correct[UZ.results$n==i])[1]
	Corrects$termCorrect[Corrects$n==i]<-table(UZ.results$termCorrect[UZ.results$n==i])[2]
	Corrects$termIncorrect[Corrects$n==i]<-table(UZ.results$termCorrect[UZ.results$n==i])[1]
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



Corrects$noModel<-100-(as.numeric(Corrects$Correct)+as.numeric(Corrects$Incorrect))
Corrects$TermnoModel<-100-(as.numeric(Corrects$termCorrect)+as.numeric(Corrects$termIncorrect))
Corrects[,1]<-as.numeric(Corrects[,1])
Corrects[,2]<-as.numeric(Corrects[,2])
Corrects[,3]<-as.numeric(Corrects[,3])
Corrects[,4]<-as.numeric(Corrects[,4])
Corrects[,7]<-as.numeric(Corrects[,7])




Corrects.melt<-melt(Corrects,id=c("n","model"))
Corrects.melt$value<-as.numeric(as.character(Corrects.melt$value))/100
Corrects.melt$param<-"Initiation"
Corrects.melt$param[25:48]<-"Termination"
Corrects.melt$param[61:72]<-"Termination"
Corrects.melt$param<-factor(Corrects.melt$param)
Corrects.melt$variable[25:36]<-"Correct"
Corrects.melt$variable[37:48]<-"Incorrect"
Corrects.melt$variable[61:72]<-"noModel"

ggplot(Corrects.melt,aes(x=n,y=value,colour=variable))+geom_point()+
facet_grid(Corrects.melt$param~.)+
 binomial_smooth(se=FALSE)+
 xlim(c(0,100))+
theme_bw()

### also plot difference between model parameter estimate and true parameter



####fitting hinge models to data



hinge.mat<-matrix(nrow=length(unique(sims[[1]]$V3)),ncol=7)
rownames(hinge.mat)<-unique(simDataLong$V3)
colnames(hinge.mat)<-c("n","startDate","endDate","startmin","startmax","endmin","endmax")
hinge.mat<-data.frame(hinge.mat)
hinge.mat$n<-c(5,10,15,20,30,40,50,75,100,150,200,250)


hinge<-list()

hinges<-list()
for(i in 1:100){
	hinges[[i]]<-list()
}

for(j in 1:100){
	cat("\n")
	cat("Sim",j,"/100","\n")
for(i in unique(sims[[j]]$V3)){
	tryCatch({
		cat(i," ")
      hinges[[j]][i]<-summary(double.hinge(x=sims[[j]][sims[[j]]$V3==i,]$julian,y=sims[[j]][sims[[j]]$V3==i,]$pscore, lower.y=0, upper.y=1, var.type="bootstrap"))
  		})
	}
}

hinge.mats<-list()

for (i in 1:100){
	hinge.mats[[i]]<-matrix(nrow=length(unique(sims[[1]]$V3)),ncol=7)
	colnames(hinge.mats[[i]])<-c("n","startDate","endDate","startmin","startmax","endmin","endmax")
		for (j in 1: nrow(hinge.mats[[i]])){
			hinge.mats[[i]][j,2]<-unlist(hinges[[i]][j])[1]
			hinge.mats[[i]][j,3]<-unlist(hinges[[i]][j])[2]
			hinge.mats[[i]][j,4]<-unlist(hinges[[i]][j])[7]
			hinge.mats[[i]][j,5]<-unlist(hinges[[i]][j])[10]
			hinge.mats[[i]][j,6]<-unlist(hinges[[i]][j])[8]
			hinge.mats[[i]][j,7]<-unlist(hinges[[i]][j])[11]
	}
hinge.mats[[i]]<-data.frame(hinge.mats[[i]])
hinge.mats[[i]]$n<-c(5,10,15,20,30,40,50,75,100,150,200,250)
hinge.mats[[i]]$startCorrect<-hinge.mats[[i]]$startmin<151&hinge.mats[[i]]$startmax>149
hinge.mats[[i]]$endCorrect<-hinge.mats[[i]]$endmin<201&hinge.mats[[i]]$endmax>199
}	




hinges.all<-bind_rows(hinge.mats)

Corrects.hinges<-data.frame(matrix(nrow=length(unique(hinges.all$n)),ncol=6))
colnames(Corrects.hinges)<-c("Correct","Incorrect","termCorrect","termIncorrect","model","n")
Corrects.hinges$n<-c(5,10,15,20,30,40,50,75,100,150,200,250)
Corrects.hinges$model<-"hinge"

for (i in unique(hinges.all$n)){
		tryCatch({
	num<-as.numeric(i)
	Corrects.hinges$Correct[Corrects.hinges$n==i]<-sum(hinges.all$startCorrect[hinges.all$n==i]==TRUE)
	Corrects.hinges$Incorrect[Corrects.hinges$n==i]<-sum(hinges.all$startCorrect[hinges.all$n==i]==FALSE)
	Corrects.hinges$termCorrect[Corrects.hinges$n==i]<-sum(hinges.all$endCorrect[hinges.all$n==i]==TRUE)
	Corrects.hinges$termIncorrect[Corrects.hinges$n==i]<-sum(hinges.all$endCorrect[hinges.all$n==i]==FALSE)
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

Corrects.hinges$"NA"<-"0"
Corrects.hinges$noModel<-"0"
Corrects.hinges$TermnoModel<-"0"

Corrects.all<-data.frame(matrix(nrow=nrow(Corrects)+nrow(Corrects.hinges),ncol=8));colnames(Corrects.all)<-colnames(Corrects)
Corrects.all[1:12,]<-Corrects
Corrects.all[13:24,]<-Corrects.hinges


Corrects.all[,1]<-as.numeric(Corrects.all[,1])
Corrects.all[,2]<-as.numeric(Corrects.all[,2])
Corrects.all[,3]<-as.numeric(Corrects.all[,3])
Corrects.all[,4]<-as.numeric(Corrects.all[,4])
Corrects.all[,5]<-factor(Corrects.all[,5])
Corrects.all[,6]<-as.numeric(Corrects.all[,6])
Corrects.all[,7]<-as.numeric(Corrects.all[,7])
Corrects.all[,8]<-as.numeric(Corrects.all[,8])


write.csv(Corrects.all,file="~/Dropbox/HingeModels/Simuation_numbers_correct.csv",row.names=FALSE)

CA.n<-read.csv("~/Dropbox/HingeModels/Simuation_numbers_correct.csv")

n.melt<-melt(CA.n,id=c("n","model"))


n.melt$value<-as.numeric(as.character(n.melt$value))/100
n.melt$param<-"Initiation"


n.melt$variable<-"Correct"
n.melt$variable[25:48]<-"Incorrect"
n.melt$variable[73:96]<-"Incorrect"
n.melt$variable[97:144]<-"noModel"
n.melt$param<-"Initiation"
n.melt$param[49:96]<-"Termination"



p<-ggplot(n.melt,aes(x=n,y=value,colour=variable))+geom_point()+
facet_grid(n.melt$param~model)+
#binomial_smooth(se=FALSE,lwd=.5)+
#geom_spline(alpha=.5,all.knots=FALSE,nknots=6)+ 
geom_line(stat="smooth",method = "loess",alpha = 0.3,lwd=.5)+
ylim(0,1)+
ylab("% of models")+
xlab("sample size")+
 #xlim(c(0,100))+
theme_bw()

pdf("~/Dropbox/HingeModels/UZ_hinge_correct_param.pdf")
print(p)
dev.off()





UZ.start<-data.frame(UZ.results$startDate,UZ.results$n);colnames(UZ.start)<-c("value","n");UZ.start$param<-"start";UZ.start$model="UZ"
UZ.end<-data.frame(UZ.results$startDate+UZ.results$duration,UZ.results$n);colnames(UZ.end)<-c("value","n");UZ.end$param<-"end";UZ.end$model="UZ"
hinge.start<-data.frame(hinges.all$startDate,hinges.all$n);colnames(hinge.start)<-c("value","n");hinge.start$param<-"start";hinge.start$model="hinge"
hinge.end<-data.frame(hinges.all$endDate,hinges.all$n);colnames(hinge.end)<-c("value","n");hinge.end$param<-"end";hinge.end$model="hinge"

params<-rbind(UZ.start,UZ.end,hinge.start,hinge.end)


q<-ggplot(params,aes(x=n,y=value))+geom_point(alpha=.05)+
facet_grid(param~model)+
theme_bw()

pdf("~/Dropbox/HingeModels/UZ_hinge_param_estimates.pdf")
print(q)
dev.off()



##########################################################
##########################################################
######## Ok now the same thing but for different variances
##########################################################
##########################################################

v.out<-matrix(nrow=12,ncol=6)

colnames(v.out)<-c("input_init","input_duration","input_var","model_init","model_duration","model_var")
rownames(v.out)<-c("1","2","5","10","15","20","30","40","50","100","150","200")

##start with init=150, dur = 50, var = 10

v.samples<-c(1,2,3,5,8,13,21,34,55,89,144)


v.simData<-list()



v.molt.oneSim<-function(){
for (i in v.samples){
	v.simData[[as.character(i)]]<-molt.simulate(n=50,init=150,duration=100,var=i)
	}
for (i in v.samples){
	v.simData[[as.character(i)]][,3]<-as.character(i)
	}
	return(v.simData)
}


v.sims<-list()

system.time(
for(j in 1:100){
v.sims[[as.character(j)]]<-v.molt.oneSim()
}
)

v.sim1<-bind_rows(v.sims[[1]])

for(i in 1:100){
	v.sims[[i]]<-bind_rows(v.sims[[i]])
}

v.simsAll<-bind_rows(v.sims)
colnames(v.simsAll)[3]<-"variance"

write.csv(v.simsAll,"~/Dropbox/HingeModels/simulatedVarianceData.csv")

v.simsAll<-read.csv("~/Dropbox/HingeModels/simulatedVarianceData.csv")

v.sim1$V3<-factor(v.sim1$V3,levels=c(1,2,3,5,8,13,21,34,55,89,144))

r<-ggplot(v.sim1,aes(x=julian,y=pscore))+geom_point(alpha=.5)+facet_wrap(~V3)+theme_bw()

pdf("~/Dropbox/HingeModels/varianceSimulations.pdf")
print(r)
dev.off()


sims<-v.sims

for(i in 1:100){
sims[[i]]$pscore<-sims[[i]]$pscore/100
}


#this is wrapped in tryCatch to keep the loop going after "moult" spits an error - because too few molting birds
#will break the function and the loop

fits<-list()
fit<-list()
fit$'1'<-"na"
fit$'2'<-"na"
fit$'3'<-"na"
fit$'5'<-"na"
fit$'8'<-"na"
fit$'13'<-"na"
fit$'21'<-"na"
fit$'34'<-"na"
fit$'55'<-"na"
fit$'89'<-"na"
fit$'144'<-"na"

for(j in 1:100){
for(i in unique(sims[[j]]$V3)){
	tryCatch({
fit[[i]]<-moult(pscore~julian,data=sims[[j]][sims[[j]]$V3==i,])	
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	fits[[j]]<-fit
	}
}



# no need to define a list for it, see below
#se<-list()
ci<-list()
ci_end<-list()
#se_end<-list()
#ses<-list()
cis<-list()
cis_end<-list()
#ses_end<-list()
for(j in 1:length(fits)){
    for(i in 1:11){
    					tryCatch({
        # define a fit variable and reuse it later
        fit=fits[[j]][[i]]
        # se is a vector and you can use it directly. 
        se=sqrt(diag(vcov(fit)))
        ci[[i]]<-cbind(coef(fit)-1.96*se, coef(fit)+1.96*se)
       
        comb=c(1,1,0)        # the following matrix operation is simply doing coef(fit)[1] + coef(fit)[2] , but it is convenient to do it in matrix operation because the variance can be computed with it as well
        est =  comb%*%coef(fit)
        se_end=sqrt(comb%*%vcov(fit)%*%comb)      
        # do you need as.numeric?
        ci_end[[i]]<-cbind(as.numeric(est-1.96*se_end), as.numeric(est+1.96*se_end))
             			        		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    #ses[[j]]<-se

    names(ci)<-names(fits[[1]])
    cis[[j]]<-ci
    names(ci_end)<-names(fits[[1]])
    cis_end[[j]]<-ci_end
}




UZ.mat<-list()

for (i in 1:100){
UZ.mat[[i]]<-matrix(nrow=11,ncol=8)
rownames(UZ.mat[[i]])<-unique(sims[[i]]$V3)
colnames(UZ.mat[[i]])<-c("duration_min","duration_max","start_min","start_max","SD_min","SD_max","term_min","term_max")
}


for(j in 1:100){
	for(i in 1:11){
				tryCatch({
    		 	UZ.mat[[j]][i,1]<-cis[[j]][[i]][,1][1]
      			UZ.mat[[j]][i,2]<-cis[[j]][[i]][,2][1]
     			UZ.mat[[j]][i,3]<-cis[[j]][[i]][,1][2]
     			UZ.mat[[j]][i,4]<-cis[[j]][[i]][,2][2]
     			UZ.mat[[j]][i,5]<-cis[[j]][[i]][,1][3]
     			UZ.mat[[j]][i,6]<-cis[[j]][[i]][,2][3]
     			UZ.mat[[j]][i,7]<-cis_end[[j]][[i]][,1]
     			UZ.mat[[j]][i,8]<-cis_end[[j]][[i]][,2]
     			        		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
					}
		UZ.mat[[j]]<-data.frame(UZ.mat[[j]])
		UZ.mat[[j]]$var<-rownames(UZ.mat[[j]])	
}




UZ.results<-bind_rows(UZ.mat)





write.table(UZ.results,"~//Dropbox/HingeModels/UZ_var.csv",sep=",")







#Now calculate how many UZ models got the answer within the CI



UZ.results$correct<-UZ.results$start_min<151&UZ.results$start_max>149
UZ.results$termCorrect<-UZ.results$term_min<251&UZ.results$term_max>249


#count up corrects by n

Corrects<-data.frame(matrix(nrow=length(unique(UZ.results$var)),ncol=6))
colnames(Corrects)<-c("Correct","Incorrect","termCorrect","termIncorrect","model","var")
Corrects$var<-c(1,2,3,5,8,13,21,34,55,89,144)
Corrects$model<-"UZ"

for (i in unique(UZ.results$var)){
		tryCatch({
	num<-as.numeric(i)
	Corrects$Correct[Corrects$var==i]<-table(UZ.results$correct[UZ.results$var==i])[2]
	Corrects$Incorrect[Corrects$var==i]<-table(UZ.results$correct[UZ.results$var==i])[1]
	Corrects$termCorrect[Corrects$var==i]<-table(UZ.results$termCorrect[UZ.results$var==i])[2]
	Corrects$termIncorrect[Corrects$var==i]<-table(UZ.results$termCorrect[UZ.results$var==i])[1]
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



Corrects$noModel<-100-(as.numeric(Corrects$Correct)+as.numeric(Corrects$Incorrect))
Corrects$TermnoModel<-100-(as.numeric(Corrects$termCorrect)+as.numeric(Corrects$termIncorrect))
Corrects[,1]<-as.numeric(Corrects[,1])
Corrects[,2]<-as.numeric(Corrects[,2])
Corrects[,3]<-as.numeric(Corrects[,3])
Corrects[,4]<-as.numeric(Corrects[,4])
Corrects[,7]<-as.numeric(Corrects[,7])




Corrects.melt<-melt(Corrects,id=c("var","model"))
Corrects.melt$value<-as.numeric(as.character(Corrects.melt$value))/100
Corrects.melt$param<-"Initiation"
Corrects.melt$param[23:44]<-"Termination"
Corrects.melt$param[56:66]<-"Termination"
Corrects.melt$param<-factor(Corrects.melt$param)
Corrects.melt$variable<-"Correct"
Corrects.melt$variable[12:21]<-"Incorrect"
Corrects.melt$variable[34:44]<-"Incorrect"
Corrects.melt$variable[45:66]<-"noModel"

ggplot(Corrects.melt,aes(x=var,y=value,colour=variable))+geom_point()+
facet_grid(Corrects.melt$param~.)+
 binomial_smooth(se=FALSE,lwd=.5)+
 xlim(c(0,100))+
theme_bw()

### also plot difference between model parameter estimate and true parameter



####fitting hinge models to data



hinge.mat<-matrix(nrow=length(unique(sims[[1]]$V3)),ncol=7)
rownames(hinge.mat)<-unique(simDataLong$V3)
colnames(hinge.mat)<-c("var","startDate","endDate","startmin","startmax","endmin","endmax")
hinge.mat<-data.frame(hinge.mat)
hinge.mat$var<-c(1,2,3,5,8,13,21,34,55,89,144)


hinge<-list()

hinges<-list()
for(i in 1:100){
	hinges[[i]]<-list()
}

for(j in 1:100){
	cat("\n")
	cat("Sim",j,"/100","\n")
for(i in unique(sims[[j]]$V3)){
	tryCatch({
		cat(i," ")
      hinges[[j]][i]<-summary(double.hinge(x=sims[[j]][sims[[j]]$V3==i,]$julian,y=sims[[j]][sims[[j]]$V3==i,]$pscore, lower.y=0, upper.y=1, var.type="bootstrap"))
  		})
	}
}

hinge.mats<-list()

for (i in 1:100){
	hinge.mats[[i]]<-matrix(nrow=length(unique(sims[[1]]$V3)),ncol=7)
	colnames(hinge.mats[[i]])<-c("n","startDate","endDate","startmin","startmax","endmin","endmax")
		for (j in 1: nrow(hinge.mats[[i]])){
			hinge.mats[[i]][j,2]<-unlist(hinges[[i]][j])[1]
			hinge.mats[[i]][j,3]<-unlist(hinges[[i]][j])[2]
			hinge.mats[[i]][j,4]<-unlist(hinges[[i]][j])[7]
			hinge.mats[[i]][j,5]<-unlist(hinges[[i]][j])[10]
			hinge.mats[[i]][j,6]<-unlist(hinges[[i]][j])[8]
			hinge.mats[[i]][j,7]<-unlist(hinges[[i]][j])[11]
	}
hinge.mats[[i]]<-data.frame(hinge.mats[[i]])
hinge.mats[[i]]$var<-c(1,2,3,5,8,13,21,34,55,89,144)
hinge.mats[[i]]$startCorrect<-hinge.mats[[i]]$startmin<151&hinge.mats[[i]]$startmax>149
hinge.mats[[i]]$endCorrect<-hinge.mats[[i]]$endmin<251&hinge.mats[[i]]$endmax>149
}	




hinges.all<-bind_rows(hinge.mats)

Corrects.hinges<-data.frame(matrix(nrow=length(unique(hinges.all$var)),ncol=6))
colnames(Corrects.hinges)<-c("Correct","Incorrect","termCorrect","termIncorrect","model","var")
Corrects.hinges$var<-c(1,2,3,5,8,13,21,34,55,89,144)
Corrects.hinges$model<-"hinge"

for (i in unique(hinges.all$var)){
		tryCatch({
	num<-as.numeric(i)
	Corrects.hinges$Correct[Corrects.hinges$var==i]<-sum(hinges.all$startCorrect[hinges.all$var==i]==TRUE)
	Corrects.hinges$Incorrect[Corrects.hinges$var==i]<-sum(hinges.all$startCorrect[hinges.all$var==i]==FALSE)
	Corrects.hinges$termCorrect[Corrects.hinges$var==i]<-sum(hinges.all$endCorrect[hinges.all$var==i]==TRUE)
	Corrects.hinges$termIncorrect[Corrects.hinges$var==i]<-sum(hinges.all$endCorrect[hinges.all$var==i]==FALSE)
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

Corrects.hinges$noModel<-"0"
Corrects.hinges$TermnoModel<-"0"

Corrects.all<-data.frame(matrix(nrow=nrow(Corrects)+nrow(Corrects.hinges),ncol=8));colnames(Corrects.all)<-colnames(Corrects)
Corrects.all[1:11,]<-Corrects
Corrects.all[12:22,]<-Corrects.hinges


Corrects.all[,1]<-as.numeric(Corrects.all[,1])
Corrects.all[,2]<-as.numeric(Corrects.all[,2])
Corrects.all[,3]<-as.numeric(Corrects.all[,3])
Corrects.all[,4]<-as.numeric(Corrects.all[,4])
Corrects.all[,5]<-factor(Corrects.all[,5])
Corrects.all[,6]<-as.numeric(Corrects.all[,6])
Corrects.all[,7]<-as.numeric(Corrects.all[,7])
Corrects.all[,8]<-as.numeric(Corrects.all[,8])


write.csv(Corrects.all,file="~/Dropbox/HingeModels/Simuation_variance_numbers_correct.csv",row.names=FALSE)

CA.v<-read.csv("~/Dropbox/HingeModels/Simuation_variance_numbers_correct.csv")

var.melt<-melt(CA.v,id=c("var","model"))


var.melt$value<-as.numeric(as.character(var.melt$value))/100
var.melt$param<-"Initiation"
var.melt$param[45:88]<-"Termination"


var.melt$variable<-"Correct"
var.melt$variable[23:44]<-"Incorrect"
var.melt$variable[67:88]<-"Incorrect"
var.melt$variable[89:132]<-"noModel"






v<-ggplot(var.melt,aes(x=var,y=value,colour=variable))+geom_point()+
facet_grid(var.melt$param~model)+
#binomial_smooth(se=FALSE)+
#geom_spline(alpha=.5,all.knots=FALSE,nknots=6)+ 
geom_line(stat="smooth",method = "loess",alpha = 0.3,lwd=.5)+
ylim(0,1)+
ylab("% of models")+
xlab("initiation variance")+
 #xlim(c(0,100))+
theme_bw()

pdf("~/Dropbox/HingeModels/UZ_hinge_correct_param.pdf")
print(v)
dev.off()





UZ.start<-data.frame(UZ.results$startDate,UZ.results$n);colnames(UZ.start)<-c("value","n");UZ.start$param<-"start";UZ.start$model="UZ"
UZ.end<-data.frame(UZ.results$startDate+UZ.results$duration,UZ.results$n);colnames(UZ.end)<-c("value","n");UZ.end$param<-"end";UZ.end$model="UZ"
hinge.start<-data.frame(hinges.all$startDate,hinges.all$n);colnames(hinge.start)<-c("value","n");hinge.start$param<-"start";hinge.start$model="hinge"
hinge.end<-data.frame(hinges.all$endDate,hinges.all$n);colnames(hinge.end)<-c("value","n");hinge.end$param<-"end";hinge.end$model="hinge"

params<-rbind(UZ.start,UZ.end,hinge.start,hinge.end)


q.v<-ggplot(params,aes(x=n,y=value))+geom_point(alpha=.05)+
facet_grid(param~model)+
theme_bw()

pdf("~/Dropbox/HingeModels/UZ__variance_hinge_param_estimates.pdf")
print(q.v)
dev.off()



################
################
### now duration
################
###############

d.out<-matrix(nrow=12,ncol=6)

colnames(d.out)<-c("input_init","input_duration","input_var","model_init","model_duration","model_dur")
rownames(d.out)<-c("1","2","5","10","15","20","30","40","50","100","150","200","250")

##start with init=150, dur = 50, var = 10

d.samples<-c(1,2,5,10,20,30,40,50,100,150,200,250)


d.simData<-list()



d.molt.oneSim<-function(){
for (i in d.samples){
	d.simData[[as.character(i)]]<-molt.simulate(n=50,init=50,duration=i,var=10)
	}
for (i in d.samples){
	d.simData[[as.character(i)]][,3]<-as.character(i)
	}
	return(d.simData)
}


d.sims<-list()

system.time(
for(j in 1:100){
d.sims[[as.character(j)]]<-d.molt.oneSim()
}
)

d.sim1<-bind_rows(d.sims[[1]])

for(i in 1:100){
	d.sims[[i]]<-bind_rows(d.sims[[i]])
}

d.simsAll<-bind_rows(d.sims)
colnames(d.simsAll)[3]<-"durations"

write.csv(d.simsAll,"~/Dropbox/HingeModels/simulatedDurationData.csv")

d.sim1$V3<-factor(d.sim1$V3,levels=c(1,2,5,10,20,30,40,50,100,150,200,250))

s<-ggplot(d.sim1,aes(x=julian,y=pscore))+geom_point(alpha=.5)+facet_wrap(~V3)+theme_bw()

pdf("~/Dropbox/HingeModels/durationSimulations.pdf")
print(s)
dev.off()


sims<-d.sims

for(i in 1:100){
sims[[i]]$pscore<-sims[[i]]$pscore/100
}


#this is wrapped in tryCatch to keep the loop going after "moult" spits an error - because too few molting birds
#will break the function and the loop

fits<-list()
fit<-list()
fit$'1'<-"na"
fit$'2'<-"na"
fit$'5'<-"na"
fit$'10'<-"na"
fit$'20'<-"na"
fit$'30'<-"na"
fit$'40'<-"na"
fit$'50'<-"na"
fit$'100'<-"na"
fit$'150'<-"na"
fit$'200'<-"na"
fit$'250'<-"na"

for(j in 1:100){
for(i in unique(sims[[j]]$V3)){
	tryCatch({
fit[[i]]<-moult(pscore~julian,data=sims[[j]][sims[[j]]$V3==i,])	
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	fits[[j]]<-fit
	}
}



# no need to define a list for it, see below
#se<-list()
ci<-list()
ci_end<-list()
#se_end<-list()
#ses<-list()
cis<-list()
cis_end<-list()
#ses_end<-list()
for(j in 1:length(fits)){
    for(i in 1:12){
    					tryCatch({
        # define a fit variable and reuse it later
        fit=fits[[j]][[i]]
        # se is a vector and you can use it directly. 
        se=sqrt(diag(vcov(fit)))
        ci[[i]]<-cbind(coef(fit)-1.96*se, coef(fit)+1.96*se)
       
        comb=c(1,1,0)        # the following matrix operation is simply doing coef(fit)[1] + coef(fit)[2] , but it is convenient to do it in matrix operation because the variance can be computed with it as well
        est =  comb%*%coef(fit)
        se_end=sqrt(comb%*%vcov(fit)%*%comb)      
        # do you need as.numeric?
        ci_end[[i]]<-cbind(as.numeric(est-1.96*se_end), as.numeric(est+1.96*se_end))
             			        		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    #ses[[j]]<-se

    names(ci)<-names(fits[[1]])
    cis[[j]]<-ci
    names(ci_end)<-names(fits[[1]])
    cis_end[[j]]<-ci_end
}




UZ.mat<-list()

for (i in 1:100){
UZ.mat[[i]]<-matrix(nrow=12,ncol=8)
rownames(UZ.mat[[i]])<-unique(sims[[i]]$V3)
colnames(UZ.mat[[i]])<-c("duration_min","duration_max","start_min","start_max","SD_min","SD_max","term_min","term_max")
}


for(j in 1:100){
	for(i in 1:12){
				tryCatch({
    		 	UZ.mat[[j]][i,1]<-cis[[j]][[i]][,1][1]
      			UZ.mat[[j]][i,2]<-cis[[j]][[i]][,2][1]
     			UZ.mat[[j]][i,3]<-cis[[j]][[i]][,1][2]
     			UZ.mat[[j]][i,4]<-cis[[j]][[i]][,2][2]
     			UZ.mat[[j]][i,5]<-cis[[j]][[i]][,1][3]
     			UZ.mat[[j]][i,6]<-cis[[j]][[i]][,2][3]
     			UZ.mat[[j]][i,7]<-cis_end[[j]][[i]][,1]
     			UZ.mat[[j]][i,8]<-cis_end[[j]][[i]][,2]
     			        		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
					}
		UZ.mat[[j]]<-data.frame(UZ.mat[[j]])
		UZ.mat[[j]]$n<-rownames(UZ.mat[[j]])	
}




UZ.results<-bind_rows(UZ.mat)



pdf("~/Dropbox/HingeModels/UZ.meanSims.pdf")
UZ.results$n<-factor(UZ.results$n,levels=c(1,2,5,10,20,30,40,50,100,150,200,250)
ggplot(UZ.results,aes(x=startDate))+
geom_histogram(bins=50,col="gray",fill="gray")+
facet_wrap(.~UZ.results$n)+
geom_segment(aes(x=150,y=0,xend=150,yend=40),linetype="dashed")+
theme_bw()
dev.off()


write.table(UZ.results,"~//Dropbox/HingeModels/UZ.csv",sep=",")







#Now calculate how many UZ models got the answer within the CI



UZ.results$Correct<-NA
UZ.results$termCorrect<-NA

for (i in unique(UZ.results$n)){
	UZ.results[UZ.results$n==i,]$Correct<-UZ.results[UZ.results$n==i,]$start_min<50&UZ.results[UZ.results$n==i,]$start_max>50
	UZ.results[UZ.results$n==i,]$termCorrect<-UZ.results[UZ.results$n==i,]$duration_min<as.numeric(i)&UZ.results[UZ.results$n==i,]$duration_max>as.numeric(i)
		}



#count up corrects by n

Corrects<-data.frame(matrix(nrow=length(unique(UZ.results$n)),ncol=6))
colnames(Corrects)<-c("Correct","Incorrect","model","n","termCorrect","termIncorrect")
Corrects$n<-c(1,2,5,10,20,30,40,50,100,150,200,250)
Corrects$model<-"UZ"

for (i in unique(UZ.results$n)){
		tryCatch({
	num<-as.numeric(i)
	Corrects$Correct[Corrects$n==i]<-table(UZ.results$Correct[UZ.results$n==i])[2]
	Corrects$Incorrect[Corrects$n==i]<-table(UZ.results$Correct[UZ.results$n==i])[1]
	Corrects$termCorrect[Corrects$n==i]<-table(UZ.results$termCorrect[UZ.results$n==i])[2]
	Corrects$termIncorrect[Corrects$n==i]<-table(UZ.results$termCorrect[UZ.results$n==i])[1]
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



Corrects$noModel<-100-(as.numeric(Corrects$Correct)+as.numeric(Corrects$Incorrect))
Corrects[,1]<-as.numeric(Corrects[,1])
Corrects[,2]<-as.numeric(Corrects[,2])
Corrects[,3]<-as.numeric(Corrects[,3])
Corrects[,4]<-as.numeric(Corrects[,4])
Corrects$model="UZ"




Corrects.melt<-melt(Corrects,id=c("n","model"))
Corrects.melt$value<-as.numeric(as.character(Corrects.melt$value))/100
Corrects.melt$variable<-"Correct"
Corrects.melt$variable[13:24]<-"Incorrect"
Corrects.melt$variable[25:36]<-"noModel"

ggplot(Corrects.melt,aes(x=n,y=value,colour=variable))+geom_point()+
 binomial_smooth(se=FALSE)+
 xlim(c(0,100))+
theme_bw()

### also plot difference between model parameter estimate and true parameter



####fitting hinge models to data



hinge.mat<-matrix(nrow=length(unique(sims[[1]]$V3)),ncol=7)
rownames(hinge.mat)<-unique(simDataLong$V3)
colnames(hinge.mat)<-c("n","startDate","endDate","startmin","startmax","endmin","endmax")
hinge.mat<-data.frame(hinge.mat)
hinge.mat$n<-c(1,2,5,10,20,30,40,50,100,150,200,250)


hinge<-list()

hinges<-list()
for(i in 1:100){
	hinges[[i]]<-list()
}

for(j in 1:100){
	cat("\n")
	cat("Sim",j,"/100","\n")
for(i in unique(sims[[j]]$V3)){
	tryCatch({
		cat(i," ")
      hinges[[j]][i]<-summary(double.hinge(x=sims[[j]][sims[[j]]$V3==i,]$julian,y=sims[[j]][sims[[j]]$V3==i,]$pscore, lower.y=0, upper.y=1, var.type="bootstrap"))
  		})
	}
}

hinge.mats<-list()

for (i in 1:100){

	hinge.mats[[i]]<-matrix(nrow=length(unique(sims[[1]]$V3)),ncol=7)
	colnames(hinge.mats[[i]])<-c("n","startDate","endDate","startmin","startmax","endmin","endmax")
		for (j in 1: nrow(hinge.mats[[i]])){
			hinge.mats[[i]][j,2]<-unlist(hinges[[i]][j])[1]
			hinge.mats[[i]][j,3]<-unlist(hinges[[i]][j])[2]
			hinge.mats[[i]][j,4]<-unlist(hinges[[i]][j])[7]
			hinge.mats[[i]][j,5]<-unlist(hinges[[i]][j])[10]
			hinge.mats[[i]][j,6]<-unlist(hinges[[i]][j])[8]
			hinge.mats[[i]][j,7]<-unlist(hinges[[i]][j])[11]
	}
	hinge.mats[[i]]<-data.frame(hinge.mats[[i]])
	hinge.mats[[i]]$n<-c(1,2,5,10,20,30,40,50,100,150,200,250)
	hinge.mats[[i]]$correct<-NA
	hinge.mats[[i]]$termCorrect<-NA

}	







for(i in 1:100){
		for(k in unique(hinge.mats[[1]]$n)){
			hinge.mats[[i]][hinge.mats[[i]]$n==k,]$correct<-hinge.mats[[i]][hinge.mats[[i]]$n==k,]$startmin<50&hinge.mats[[i]][hinge.mats[[i]]$n==k,]$startmax>50
			hinge.mats[[i]][hinge.mats[[i]]$n==k,]$termCorrect<-hinge.mats[[i]][hinge.mats[[i]]$n==k,]$endmin<(50+as.numeric(k))&hinge.mats[[i]][hinge.mats[[i]]$n==k,]$endmax>(50+as.numeric(k))
		}
}


hinges.all<-bind_rows(hinge.mats)

Corrects.hinges<-data.frame(matrix(nrow=length(unique(hinges.all$n)),ncol=6))
colnames(Corrects.hinges)<-c("Correct","Incorrect","model","n","termCorrect","termIncorrect")
Corrects.hinges$n<-c(1,2,5,10,20,30,40,50,100,150,200,250)



Corrects.hinges$model<-"hinge"

for (i in unique(hinges.all$n)){
		tryCatch({
	num<-as.numeric(i)
	Corrects.hinges$Correct[Corrects.hinges$n==i]<-sum(hinges.all$correct[hinges.all$n==i]==TRUE)
	Corrects.hinges$Incorrect[Corrects.hinges$n==i]<-sum(hinges.all$correct[hinges.all$n==i]==FALSE)
	Corrects.hinges$termCorrect[Corrects.hinges$n==i]<-sum(hinges.all$termCorrect[hinges.all$n==i]==TRUE)
	Corrects.hinges$termIncorrect[Corrects.hinges$n==i]<-sum(hinges.all$termCorrect[hinges.all$n==i]==FALSE)

			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

Corrects.hinges$"NA"<-"0"
Corrects.hinges$noModel<-"0"
Corrects.hinges$TermnoModel<-"0"

Corrects.all<-data.frame(matrix(nrow=nrow(Corrects)+nrow(Corrects.hinges),ncol=7));colnames(Corrects.all)<-colnames(Corrects)
Corrects.all[1:12,]<-Corrects
Corrects.all[13:24,]<-Corrects.hinges


Corrects.all[,1]<-as.numeric(Corrects.all[,1])
Corrects.all[,2]<-as.numeric(Corrects.all[,2])
Corrects.all[,3]<-factor(Corrects.all[,3])
Corrects.all[,4]<-as.numeric(Corrects.all[,4])
Corrects.all[,5]<-as.numeric(Corrects.all[,5])
Corrects.all[,6]<-as.numeric(Corrects.all[,6])
Corrects.all[,7]<-as.numeric(Corrects.all[,7])



write.csv(Corrects.all,file="~/Dropbox/HingeModels/Simuation_duration_numbers_correct.csv",row.names=FALSE)

CA.n<-read.csv("~/Dropbox/HingeModels/Simuation_duration_numbers_correct.csv")

dur.melt<-melt(CA.n,id=c("n","model"))


dur.melt$value<-as.numeric(as.character(dur.melt$value))/100
dur.melt$param<-"Initiation"


dur.melt$variable<-"Correct"
dur.melt$variable[25:48]<-"Incorrect"
dur.melt$variable[73:96]<-"Incorrect"
dur.melt$variable[97:120]<-"noModel"
dur.melt$param<-"Initiation"
dur.melt$param[49:96]<-"Termination"



d<-ggplot(dur.melt,aes(x=n,y=value,colour=variable))+geom_point()+
facet_grid(dur.melt$param~model)+
#binomial_smooth(se=FALSE,alpha=.3)+
#geom_spline(alpha=.5,all.knots=FALSE,nknots=6)+ 
geom_line(stat="smooth",method = "loess",alpha = 0.3,lwd=.5)+
ylim(0,1)+
ylab("% of models")+
xlab("molt duration")+
 #xlim(c(0,100))+
theme_bw()

pdf("~/Dropbox/HingeModels/UZ_hinge_correct_param.pdf")
print(d)
dev.off()






pdf("~/Dropbox/HingeModels/all_params.pdf")
grid.arrange(p,v,d)
dev.off()

UZ.start<-data.frame(UZ.results$startDate,UZ.results$n);colnames(UZ.start)<-c("value","n");UZ.start$param<-"start";UZ.start$model="UZ"
UZ.end<-data.frame(UZ.results$startDate+UZ.results$duration,UZ.results$n);colnames(UZ.end)<-c("value","n");UZ.end$param<-"end";UZ.end$model="UZ"
hinge.start<-data.frame(hinges.all$startDate,hinges.all$n);colnames(hinge.start)<-c("value","n");hinge.start$param<-"start";hinge.start$model="hinge"
hinge.end<-data.frame(hinges.all$endDate,hinges.all$n);colnames(hinge.end)<-c("value","n");hinge.end$param<-"end";hinge.end$model="hinge"

params<-rbind(UZ.start,UZ.end,hinge.start,hinge.end)


q<-ggplot(params,aes(x=n,y=value))+geom_point(alpha=.05)+
facet_grid(param~model)+
theme_bw()

pdf("~/Dropbox/HingeModels/UZ_hinge_param_estimates.pdf")
print(q)
dev.off()

p.leg<-
p<-ggplot(n.melt,aes(x=n,y=value,colour=variable))+geom_point(size=.5)+
facet_grid(n.melt$param~model)+
binomial_smooth(se=FALSE)

library(cowplot)
legend <- get_legend(p.leg)


p<-ggplot(n.melt,aes(x=n,y=value,colour=variable))+geom_point(size=1)+
facet_grid(n.melt$param~model)+
binomial_smooth(se=FALSE,lwd=.5)+
#geom_spline(alpha=.5,all.knots=FALSE,nknots=6)+ 
#geom_line(stat="smooth",method = "loess",alpha = 0.3,lwd=.5)+
ylim(0,1)+
ylab("% of models")+
xlab("sample size")+
geom_rangeframe(col="black")+
theme_tufte()+
theme(legend.position = "none")




v<-ggplot(var.melt,aes(x=var,y=value,colour=variable))+geom_point(size=1)+
facet_grid(var.melt$param~model)+
binomial_smooth(se=FALSE,lwd=.5)+
#geom_spline(alpha=.5,all.knots=FALSE,nknots=6)+ 
#geom_line(stat="smooth",method = "loess",alpha = 0.3,lwd=.5)+
ylim(0,1)+
ylab("")+
xlab("initiation variance")+
geom_rangeframe(col="black")+
theme_tufte()+
theme(legend.position = "none")
 #xlim(c(0,100))+



d<-ggplot(dur.melt,aes(x=n,y=value,colour=variable))+geom_point(size=1)+
facet_grid(dur.melt$param~model)+
binomial_smooth(se=FALSE,lwd=.5)+
#geom_spline(alpha=.5,all.knots=FALSE,nknots=6)+ 
#geom_line(stat="smooth",method = "loess",alpha = 0.3,lwd=.5)+
ylim(0,1)+
ylab("")+
xlab("molt duration")+
geom_rangeframe(col="black")+
theme_tufte()+
theme(legend.position = "none")




 #xlim(c(0,100))+



g<-grid.arrange(p,v,d,legend,ncol=4)

pdf("~/Dropbox/HingeModels/all_sims.pdf",height=6,width=12)
grid.arrange(p,v,d,legend,ncol=4)
dev.off()

p.small<-ggplot(n.melt,aes(x=n,y=value,colour=variable))+geom_point(size=2)+
facet_grid(n.melt$param~model)+
binomial_smooth(se=FALSE,lwd=.5)+
geom_rangeframe(col="black")+
theme_tufte()+
#geom_spline(alpha=.5,all.knots=FALSE,nknots=6)+ 
#geom_line(stat="smooth",method = "loess",alpha = 0.3)+
ylim(0,1)+
xlim(5,50)+
ylab("% of models")+
xlab("sample size")+

pdf("~/Dropbox/HingeModels/small_sample.pdf",height=6,width=12)
print(p.small)
dev.off()



####Compare UZ to Hinge models
#FOR UZ  MODELS: convert SD to SE or (better) compute confidence intervals
##se = sd/squrt(n)

UZ.mat$SE<-UZ.mat$SD/sqrt(as.numeric(UZ.mat$n))


UZ.melt<-melt(UZ.mat)

#plot UZ and hinge estimate with confidence intervals and input parameters

##add raw data to top of grid, then four model plots below, and line up the four plots so the scales are the same

p1<-ggplot()+
geom_pointrange(data=UZ.mat,alpha=.5,aes(x=factor(as.numeric(n)),y=startDate,ymin=startDate-(SE*1.96),ymax=startDate+(SE*1.96)))+
theme_bw()+
geom_hline(yintercept=150,linetype="dashed")+
ggtitle("UZ estimate of molt initation")+
xlab("sample size")+ylab("Julian Day")+
theme(plot.title = element_text(hjust = 0.5))

p2<-ggplot()+
geom_pointrange(data=UZ.mat,alpha=.5,aes(x=factor(as.numeric(n)),y=startDate+duration,ymin=(startDate+duration)-(SE*1.96),ymax=(startDate+duration)+(SE*1.96)))+
theme_bw()+
geom_hline(yintercept=200,linetype="dashed")+
ggtitle("UZ estimate of molt completion")+
xlab("sample size")+ylab("Julian Day")+
theme(plot.title = element_text(hjust = 0.5))


p3<-ggplot()+
geom_pointrange(data=hinge.mat,alpha=.5,aes(x=factor(as.numeric(n)),y=startDate,ymin=(startmin),ymax=(startmax)))+
theme_bw()+
geom_hline(yintercept=150,linetype="dashed")+
ggtitle("Hinge estimate of molt initiation")+
xlab("sample size")+ylab("Julian Day")+
theme(plot.title = element_text(hjust = 0.5))


p4<-ggplot()+
geom_pointrange(data=hinge.mat,alpha=.5,aes(x=factor(as.numeric(n)),y=endDate,ymin=(endmin),ymax=(endmax)))+
theme_bw()+
geom_hline(yintercept=200,linetype="dashed")+
ggtitle("Hinge estimate of molt completion")+
xlab("sample size")+ylab("Julian Day")+
theme(plot.title = element_text(hjust = 0.5))

p5<-ggplot(simDataLong,aes(x=julian,y=pscore))+geom_point(size=.5,)+facet_wrap(~V3)+theme_bw()+
ggtitle("simulated molt data")+
geom_vline(xintercept=150,linetype="dashed")+
geom_vline(xintercept=200,linetype="dashed")+
xlab("Julian Day")+ylab("Molt completion %")+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
theme(plot.title = element_text(hjust = 0.5))

gs<-(c(as_grob(p5),as_grob(p2),as_grob(p4),as_grob(p1),as_grob(p3)))


p<-grid.arrange(arrangeGrob(p5),arrangeGrob(p2,p4,p1,p3,nrow=2),ncol=1)

pdf("/Users/rterrill/Dropbox/HingeModels/mainfig.pdf")
plot(p)
dev.off()


