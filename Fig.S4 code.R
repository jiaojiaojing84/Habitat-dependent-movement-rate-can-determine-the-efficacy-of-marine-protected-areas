#setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape)



library(ggplot2)



library(scales)



library(pheatmap)


####Open system
JAP07<-function(t, inits,parameters) {
  
  
  with(as.list(c(inits, parameters)),{
    
    AF=R-(mu+F)*A11-DM*beta*A11+DM*A22   
    
    
    AM=R-mu*A22-DM*A22+DM*beta*A11
    
    
    list(c(AF,AM))
    
  })
  
}


####close system
JAP08<-function(t, inits,parameters) {
  
  
  with(as.list(c(inits, parameters)),{
    
    AF=r*A11-e*A11^2-(mu+F)*A11-DM*beta*A11+DM*A22   
    
    
    AM=r*A22-e*A22^2-mu*A22-DM*A22+DM*beta*A11
    
    
    list(c(AF,AM))
    
  })
  
}


#semi-close system
JAP09<-function(t, inits,parameters) {
  
  
  with(as.list(c(inits, parameters)),{
    
    AF=r/2*(A11+A22)-e*A11^2-(mu+F)*A11-DM*beta*A11+DM*A22   
    
    
    AM=r/2*(A11+A22)-e*A22^2-mu*A22-DM*A22+DM*beta*A11
    
    
    list(c(AF,AM))
    
  })
  
}


###differential movement is defined as D1/D2

Timesteps=500

times <- seq(0, Timesteps, by = 1)

inits <- c(A11=1,A22=1)

DM<-seq(0,10,0.1)

F=0.25


###Open system

A1.eq1<-c()
A2.eq1<-c()

A1.eq1.5<-c()
A2.eq1.5<-c()


A1.eq2<-c()
A2.eq2<-c()

A_bef<-c()
####close system

B1.eq1<-c()
B2.eq1<-c()

B1.eq2<-c()
B2.eq2<-c()

B1.eq5<-c()
B2.eq5<-c()

####Semi-close system

C1.eq1<-c()
C2.eq1<-c()

C1.eq1.2<-c()
C2.eq1.2<-c()

C1.eq2<-c()
C2.eq2<-c()

C1.eq5<-c()
C2.eq5<-c()

for(i in 1:length(DM))
{
  ####open system
  parameters <- c(R=2,beta=1,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP07, parms = parameters)
  
  A1.eq1[i]=out[Timesteps+1,2]    
  
  A2.eq1[i]=out[Timesteps+1,3]
  
  
  parameters <- c(R=2,beta=1.5,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP07, parms = parameters)
  
  A1.eq1.5[i]=out[Timesteps+1,2]    
  
  A2.eq1.5[i]=out[Timesteps+1,3]
  
  parameters <- c(R=2,beta=2,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP07, parms = parameters)
  
  A1.eq2[i]=out[Timesteps+1,2]    
  
  A2.eq2[i]=out[Timesteps+1,3]
  
  parameters_before <- c(R=2,beta=0,mu=0.5,F=0.25,DM=0) 
  out_bef= ode(y = inits, times = times, func = JAP07, parms = parameters_before)
  
  A_bef[i]<-out_bef[Timesteps+1,2]

  
  ####close system
  parameters <- c(r=1,e=0.01,beta=1,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP08, parms = parameters)
  
  B1.eq1[i]=out[Timesteps+1,2]    
  
  B2.eq1[i]=out[Timesteps+1,3]
  
  
  parameters <- c(r=1,e=0.01,beta=2,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP08, parms = parameters)
  
  B1.eq2[i]=out[Timesteps+1,2]    
  
  B2.eq2[i]=out[Timesteps+1,3]
  
  parameters <- c(r=1,e=0.01,beta=5,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP08, parms = parameters)
  
  B1.eq5[i]=out[Timesteps+1,2]    
  
  B2.eq5[i]=out[Timesteps+1,3]
  
  ####Semi-close system
  
  parameters <- c(r=1,e=0.01,beta=1,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP09, parms = parameters)
  
  C1.eq1[i]=out[Timesteps+1,2]    
  
  C2.eq1[i]=out[Timesteps+1,3]
  
  
  parameters <- c(r=1,e=0.01,beta=1.2,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP09, parms = parameters)
  
  C1.eq1.2[i]=out[Timesteps+1,2]    
  
  C2.eq1.2[i]=out[Timesteps+1,3]
  
  parameters <- c(r=1,e=0.01,beta=1.5,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP09, parms = parameters)
  
  C1.eq2[i]=out[Timesteps+1,2]    
  
  C2.eq2[i]=out[Timesteps+1,3]
  
  parameters <- c(r=1,e=0.01,beta=2,mu=0.5,F=0.25,DM=DM[i])   
  
  out= ode(y = inits, times = times, func = JAP09, parms = parameters)
  
  C1.eq5[i]=out[Timesteps+1,2]    
  
  C2.eq5[i]=out[Timesteps+1,3]
}

###before indexes
loc_beforeA<-1
loc_beforeB<-1
loc_beforeC<-1

reg_beforeA<-A_bef[1]*2
##K_F=100, mu=0.5, F=0.25
reg_beforeB<-100*(1-0.5-0.25)*2
reg_beforeC<-reg_beforeB

fis_beforeA<-F*reg_beforeA
fis_beforeB<-F*reg_beforeB
fis_beforeC<-F*reg_beforeC
  
#local effect
locA1<-c()
locA1.5<-c()
locA2<-c()

locB1<-c()
locB2<-c()
locB5<-c()

locC1<-c()
locC1.2<-c()
locC2<-c()
locC5<-c()

for(i in 1:length(DM))
{
  locA1[i]<-A2.eq1[i]/A1.eq1[i]/loc_beforeA
  locA1.5[i]<-A2.eq1.5[i]/A1.eq1.5[i]/loc_beforeA
  locA2[i]<-A2.eq2[i]/A1.eq2[i]/loc_beforeA
  
  locB1[i]<-B2.eq1[i]/B1.eq1[i]/loc_beforeB
  locB2[i]<-B2.eq2[i]/B1.eq2[i]/loc_beforeB
  locB5[i]<-B2.eq5[i]/B1.eq5[i]/loc_beforeB
  
  locC1[i]<-C2.eq1[i]/C1.eq1[i]/loc_beforeC
  locC1.2[i]<-C2.eq1.2[i]/C1.eq1.2[i]/loc_beforeC
  locC2[i]<-C2.eq2[i]/C1.eq2[i]/loc_beforeC
  locC5[i]<-C2.eq5[i]/C1.eq5[i]/loc_beforeC
}

regA1<-c()
regA1.5<-c()
regA2<-c()

regB1<-c()
regB2<-c()
regB5<-c()

regC1<-c()
regC1.2<-c()
regC2<-c()
regC5<-c()

for(i in 1:length(DM))
{
  regA1[i]<-(A2.eq1[i]+A1.eq1[i])/reg_beforeA
  regA1.5[i]<-(A2.eq1.5[i]+A1.eq1.5[i])/reg_beforeA
  regA2[i]<-(A2.eq2[i]+A1.eq2[i])/reg_beforeA
  
  regB1[i]<-(B2.eq1[i]+B1.eq1[i])/reg_beforeB
  regB2[i]<-(B2.eq2[i]+B1.eq2[i])/reg_beforeB
  regB5[i]<-(B2.eq5[i]+B1.eq5[i])/reg_beforeB
  
  regC1[i]<-(C2.eq1[i]+C1.eq1[i])/reg_beforeC
  regC1.2[i]<-(C2.eq1.2[i]+C1.eq1.2[i])/reg_beforeC
  regC2[i]<-(C2.eq2[i]+C1.eq2[i])/reg_beforeC
  regC5[i]<-(C2.eq5[i]+C1.eq5[i])/reg_beforeC
  
}

yieA1<-c()
yieA1.5<-c()
yieA2<-c()

yieB1<-c()
yieB2<-c()
yieB5<-c()

yieC1<-c()
yieC1.2<-c()
yieC2<-c()
yieC5<-c()

for(i in 1:length(DM))
{
  yieA1[i]<-(A1.eq1[i]*F)/fis_beforeA
  yieA1.5[i]<-A1.eq1.5[i]*F/fis_beforeA
  yieA2[i]<-A1.eq2[i]*F/fis_beforeA
  
  yieB1[i]<-B1.eq1[i]*F/fis_beforeB
  yieB2[i]<-B1.eq2[i]*F/fis_beforeB
  yieB5[i]<-B1.eq5[i]*F/fis_beforeB
  
  yieC1[i]<-C1.eq1[i]*F/fis_beforeC
  yieC1.2[i]<-C1.eq1.2[i]*F/fis_beforeC
  yieC2[i]<-C1.eq2[i]*F/fis_beforeC
  yieC5[i]<-C1.eq5[i]*F/fis_beforeC
  
}


tiff("Fig.4_Appendix.tiff", width=5,height=5, units='in',res=600)
par(mfrow=c(3,3))
par(mar=c(1.9,2,1,2))

#plot(DM,locA1,ylim=c(min(locA1,locA1.5,locA2),max(locA1,locA1.5,locA2)),ylab="", xlab="",type="l",lwd=2)
plot(DM,locA1,ylim=c(1,2),ylab="", xlab="",type="l",lwd=2)
points(DM,locA1.5,type="l",lty=2,lwd=2)
points(DM,locA2,type="l",lty=3,lwd=2)
#dev.off()

#tiff("diff mov 9-11-17 Fig 3b.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
#plot(DM,regA1,ylim=c(min(regA1,regA1.5,regA2),max(regA1,regA1.5,regA2)),ylab="", xlab="",type="l",lwd=2)
plot(DM,regA1,ylim=c(1.2,1.3),ylab="", xlab="",type="l",lwd=2)
points(DM,regA1.5,type="l",lty=2,lwd=2)
points(DM,regA2,type="l",lty=3,lwd=2)
#dev.off()

#tiff("diff mov 9-11-17 Fig 3c.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
#plot(DM,yieA1,ylim=c(min(yieA1,yieA1.5,yieA2),max(yieA1,yieA1.5,yieA2)),ylab="", xlab="",type="l",lwd=2)
plot(DM,yieA1,ylim=c(0.4,0.6),ylab="", xlab="",type="l",lwd=2)
points(DM,yieA1.5,type="l",lty=2,lwd=2)
points(DM,yieA2,type="l",lty=3,lwd=2)


#####close system
#plot(DM,locB1,ylim=c(min(locB1,locB2,locB5),max(locB1,locB2,locB5)),ylab="", xlab="",type="l",lwd=2)
plot(DM,locB1,ylim=c(-1,5),ylab="", xlab="",type="l",lwd=2)
points(DM,locB2,type="l",lty=2,lwd=2)
points(DM,locB5,type="l",lty=3,lwd=2)
#dev.off()

#tiff("diff mov 9-11-17 Fig 3b.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
#plot(DM,regB1,ylim=c(min(regB1,regB2,regB5),max(regB1,regB2,regB5)),ylab="", xlab="",type="l",lwd=2)
plot(DM,regB1,ylim=c(1.2,1.8),ylab="", xlab="",type="l",lwd=2)
points(DM,regB2,type="l",lty=2,lwd=2)
points(DM,regB5,type="l",lty=3,lwd=2)
#dev.off()

#tiff("diff mov 9-11-17 Fig 3c.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
#plot(DM,yieB1,ylim=c(min(yieB1,yieB2,yieB5),max(yieB1,yieB2,yieB5)),ylab="", xlab="",type="l",lwd=2)
plot(DM,yieB1,ylim=c(0.2,0.8),ylab="", xlab="",type="l",lwd=2)
points(DM,yieB2,type="l",lty=2,lwd=2)
points(DM,yieB5,type="l",lty=3,lwd=2)

####semi_close system
#plot(DM,locC1,ylim=c(min(locC1,locC1.2,locC2,locC5),max(locC1,locC1.2,locC2,locC5)),ylab="", xlab="",type="l",lwd=2)
plot(DM,locC1,ylim=c(0.2,2),ylab="", xlab="",type="l",lwd=2)
points(DM,locC1.2,type="l",lty=2,lwd=2)
points(DM,locC2,type="l",lty=3,lwd=2)
points(DM,locC5,type="l",lty=5,lwd=2)
#dev.off()

#tiff("diff mov 9-11-17 Fig 3b.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
#plot(DM,regC1,ylim=c(min(regC1,regC1.2,regC2,regC5),max(regC1,regC1.2,regC2,regC5)),ylab="", xlab="",type="l",lwd=2)
plot(DM,regC1,ylim=c(1.5,1.59),ylab="", xlab="",type="l",lwd=2)
points(DM,regC1.2,type="l",lty=2,lwd=2)
points(DM,regC2,type="l",lty=3,lwd=2)
points(DM,regC5,type="l",lty=5,lwd=2)
#dev.off()

#tiff("diff mov 9-11-17 Fig 3c.tiff", width=5,height=5, units='in',res=600)
#par(mar=c(3,6,3,6))
#plot(DM,yieC1,ylim=c(min(yieC1,yieC1.2,yieC2,yieC5),max(yieC1,yieC1.2,yieC2,yieC5)),ylab="", xlab="",type="l",lwd=2)
plot(DM,yieC1,ylim=c(0.5,0.9),ylab="", xlab="",type="l",lwd=2)
points(DM,yieC1.2,type="l",lty=2,lwd=2)
points(DM,yieC2,type="l",lty=3,lwd=2)
points(DM,yieC5,type="l",lty=5,lwd=2)
dev.off()

