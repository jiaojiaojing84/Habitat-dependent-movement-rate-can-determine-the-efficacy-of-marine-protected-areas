setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape)



library(ggplot2)



library(scales)



library(pheatmap)



JAP09<-function(t, inits,parameters) {
  
  
  with(as.list(c(inits, parameters)),{
    
    AF=R-(mu+F)*A11-DM*beta*A11+DM*A22   
    
    
    AM=R-mu*A22-DM*A22+DM*beta*A11
    
    
    list(c(AF,AM))
    
  })
  
}


###differential movement is defined as D1/D2

Timesteps=500

times <- seq(0, Timesteps, by = 1)

inits <- c(A11=1,A22=1)

DM<-seq(0,10,0.1)

F=0.25

A1.eq1<-c()
A2.eq1<-c()

A1.eq1.5<-c()
A2.eq1.5<-c()


A1.eq2<-c()
A2.eq2<-c()

A_bef<-c()

for(i in 1:length(DM))
{
  #random movement
  parameters <- c(R=2,beta=1,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP09, parms = parameters)
  
  A1.eq1[i]=out[Timesteps+1,2]    
  
  A2.eq1[i]=out[Timesteps+1,3]
  
  
  parameters <- c(R=2,beta=1.5,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP09, parms = parameters)
  
  A1.eq1.5[i]=out[Timesteps+1,2]    
  
  A2.eq1.5[i]=out[Timesteps+1,3]
  
  parameters <- c(R=2,beta=2,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP09, parms = parameters)
  
  A1.eq2[i]=out[Timesteps+1,2]    
  
  A2.eq2[i]=out[Timesteps+1,3]
  
  parameters_before <- c(R=2,beta=0,mu=0.5,F=0.25,DM=0) 
  out_bef= ode(y = inits, times = times, func = JAP09, parms = parameters_before)
  
  A_bef[i]<-out_bef[Timesteps+1,2]
}

###before MPA
###local effect: =1 since there is no MPA
loc_before=1
###regional abundance:eqnmean_ba[1,1,1]=2.666667
reg_before<-A_bef[1]*2
###fishing yield
fis_before<-F*reg_before

#local effect
loc1<-c()
loc1.5<-c()
loc2<-c()

for(i in 1:length(DM))
{
  loc1[i]<-A2.eq1[i]/A1.eq1[i]/loc_before
  loc1.5[i]<-A2.eq1.5[i]/A1.eq1.5[i]/loc_before
  loc2[i]<-A2.eq2[i]/A1.eq2[i]/loc_before
}

reg1<-c()
reg1.5<-c()
reg2<-c()

for(i in 1:length(DM))
{
  reg1[i]<-(A2.eq1[i]+A1.eq1[i])/reg_before
  reg1.5[i]<-(A2.eq1.5[i]+A1.eq1.5[i])/reg_before
  reg2[i]<-(A2.eq2[i]+A1.eq2[i])/reg_before
}

yie1<-c()
yie1.5<-c()
yie2<-c()

for(i in 1:length(DM))
{
  yie1[i]<-(A1.eq1[i]*F)/fis_before
  yie1.5[i]<-(A1.eq1.5[i]*F)/fis_before
  yie2[i]<-(A1.eq2[i]*F)/fis_before
}


tiff("Open_system_Fig.3.tiff", width=5,height=5, units='in',res=600)
par(mfrow=c(1,3))
par(mar=c(12,2,12,2))

#tiff("diff mov 9-11-17 Fig 3a.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
plot(DM,loc1,ylim=c(1,2),ylab="", xlab="",type="l",lwd=2)
points(DM,loc1.5,type="l",lty=2,lwd=2)
points(DM,loc2,type="l",lty=3,lwd=2)
#legend(4,0.95,lty=c(1,2,3),c())
#dev.off()

#tiff("diff mov 9-11-17 Fig 3b.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
plot(DM,reg1,ylim=c(1.20,1.3),ylab="", xlab="",type="l",lwd=2)
points(DM,reg1.5,type="l",lty=2,lwd=2)
points(DM,reg2,type="l",lty=3,lwd=2)
#dev.off()

#tiff("diff mov 9-11-17 Fig 3c.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
plot(DM,yie1,ylim=c(0.4,0.6),ylab="", xlab="",type="l",lwd=2)
points(DM,yie1.5,type="l",lty=2,lwd=2)
points(DM,yie2,type="l",lty=3,lwd=2)
dev.off()

