####two-patch model with constant larval rain

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

DM<-c(0.1,0.5,5)
beta<-seq(1,4,0.5)

F=0.25

A1.eq0.1<-c()
A2.eq0.1<-c()

A1.eq0.5<-c()
A2.eq0.5<-c()

A1.eq5<-c()
A2.eq5<-c()

A_bef<-c()

for(i in 1:length(beta))
{
  #random movement
  parameters0.1 <- c(R=2,beta=beta[i],mu=0.5,F=0.25,DM=DM[1])   
  
  out0.1= ode(y = inits, times = times, func = JAP09, parms = parameters0.1)
  
  A1.eq0.1[i]=out0.1[Timesteps+1,2]    
  
  A2.eq0.1[i]=out0.1[Timesteps+1,3]
  
  
  parameters0.5 <- c(R=2,beta=beta[i],mu=0.5,F=0.25,DM=DM[2])   
  
  out0.5= ode(y = inits, times = times, func = JAP09, parms = parameters0.5)
  
  A1.eq0.5[i]=out0.5[Timesteps+1,2]    
  
  A2.eq0.5[i]=out0.5[Timesteps+1,3]
  
  parameters5 <- c(R=2,beta=beta[i],mu=0.5,F=0.25,DM=DM[3])   
  
  out5= ode(y = inits, times = times, func = JAP09, parms = parameters5)
  
  A1.eq5[i]=out5[Timesteps+1,2]    
  
  A2.eq5[i]=out5[Timesteps+1,3]
  
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
loc0.1<-c()
loc0.5<-c()
loc5<-c()

for(i in 1:length(beta))
{
  loc0.1[i]<-A2.eq0.1[i]/A1.eq0.1[i]/(loc_before)
  loc0.5[i]<-A2.eq0.5[i]/A1.eq0.5[i]/(loc_before)
  loc5[i]<-A2.eq5[i]/A1.eq5[i]/(loc_before)
}

reg0.1<-c()
reg0.5<-c()
reg5<-c()


for(i in 1:length(beta))
{
  reg0.1[i]<-(A2.eq0.1[i]+A1.eq0.1[i])/reg_before
  reg0.5[i]<-(A2.eq0.5[i]+A1.eq0.5[i])/reg_before
  reg5[i]<-(A2.eq5[i]+A1.eq5[i])/reg_before
}

yie0.1<-c()
yie0.5<-c()
yie5<-c()

for(i in 1:length(beta))
{
  yie0.1[i]<-A1.eq0.1[i]*F/fis_before
  yie0.5[i]<-A1.eq0.5[i]*F/fis_before
  yie5[i]<-A1.eq5[i]*F/fis_before
}


tiff("Open_system_Fig.2.tiff", width=5,height=5, units='in',res=600)
par(mfrow=c(1,3))
par(mar=c(12,2,12,2))

plot(beta,loc0.1,ylim=c(min(loc5),max(loc5)),ylab="", xlab="",type="l",lwd=2)
points(beta,loc0.5,type="l",lty=2,lwd=2)
points(beta,loc5,type="l",lty=3,lwd=2)
#legend(4,0.95,lty=c(1,2,3),c())
#dev.off()


plot(beta,reg0.1,ylim=c(min(reg5),max(reg5)),ylab="", xlab="",type="l",lwd=2)
points(beta,reg0.5,type="l",lty=2,lwd=2)
points(beta,reg5,type="l",lty=3,lwd=2)
#dev.off()

#tiff("diff mov 9-11-17 Fig 3c.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
plot(beta,yie0.1,ylim=c(min(yie5),max(yie5)),ylab="", xlab="",type="l",lwd=2)
points(beta,yie0.5,type="l",lty=2,lwd=2)
points(beta,yie5,type="l",lty=3,lwd=2)
dev.off()

