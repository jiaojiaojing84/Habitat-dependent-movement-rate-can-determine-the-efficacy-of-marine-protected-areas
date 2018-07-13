##take off exponential term +1 and interference competition 
#setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape)



library(ggplot2)



library(scales)



library(pheatmap)




################test
#####L is fishing ground while T is MPAs, total area is S
JAP08<-function(t, inits,parameters) {
  
  
  with(as.list(c(inits, parameters)),{
    
    x<-inits[1:L]
    y<-inits[(L+1):(L+T)]
    
    A<-array(NA,dim=c(1,(L+T)))
    
    #calculate larvae redistribution
    #R<-(sum(x[1:L])+sum(y[1:T]))/(T+L)
    
    if(L==1)
    {
      A[1]<-r*(sum(x[1:L])+sum(y[1:T]))/(T+L)-e*x[1]^2-(mu+F)*x[1]-D1*x[1]+D2/2*y[1]+D2/2*y[T]
    }
    else
    {
      A[1]<-r*(sum(x[1:L])+sum(y[1:T]))/(T+L)-e*x[1]^2-(mu+F)*x[1]-D1*x[1]+D1/2*x[2]+D2/2*y[1]
      
      A[L]<-r*(sum(x[1:L])+sum(y[1:T]))/(T+L)-e*x[L]^2-(mu+F)*x[L]-D1*x[L]+D1/2*x[L-1]+D2/2*y[T]
      
    }
    if(L-1>=2)
    {
      for(i in 2:(L-1))
      {
        A[i]<-r*(sum(x[1:L])+sum(y[1:T]))/(T+L)-e*x[i]^2-(mu+F)*x[i]-D1*x[i]+D1/2*(x[i-1]+x[i+1])
      }
    }
    if(T==1)
    {
      A[(L+1)]<-r*(sum(x[1:L])+sum(y[1:T]))/(T+L)-e*y[1]^2-mu*y[1]-D2*y[1]+D1/2*x[L]+D1/2*x[1]	
    }
    else
    {
      A[(L+1)]<-r*(sum(x[1:L])+sum(y[1:T]))/(T+L)-e*y[1]^2-mu*y[1]-D2*y[1]+D2/2*y[2]+D1/2*x[1]
      A[(T+L)]<-r*(sum(x[1:L])+sum(y[1:T]))/(T+L)-e*y[T]^2-mu*y[T]-D2*y[T]+D2/2*y[T-1]+D1/2*x[L]
    }
    if(T-1>=2)
    {
      for(i in (L+2):(T+L-1))
      {
        A[i]<-r*(sum(x[1:L])+sum(y[1:T]))/(T+L)-e*y[i-L]^2-mu*y[i-L]-D2*y[i-L]+D2/2*(y[i-L-1]+y[i-L+1])
      }
    }
    list(c(A))
    
  }) 
}



#################

Timesteps=500

times <- seq(0, Timesteps, by = 1)

S=10

inits <- rep(1,S)


r<-c(1,1.2,2,5)

m<-0.5

F=0.25

##h is for MPA size, r is for differential movement
eqn1<-array(NA,dim=c(length(r),length(m),S))
eqn5<-array(NA,dim=c(length(r),length(m),S))
eqn9<-array(NA,dim=c(length(r),length(m),S))



for(i in 1:length(r))
{
  for(j in 1:length(m))
  {
    for(z in 1:S)
    {
      parameters <- c(T=1,L=9,r=1,e=0.01,mu=0.5,D1=m[j]*r[i],D2=m[j],F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn1[i,j,z]<-out[Timesteps+1,z+1]
      
      parameters <- c(T=5,L=5,r=1,e=0.01,mu=0.5,D1=m[j]*r[i],D2=m[j],F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn5[i,j,z]<-out[Timesteps+1,z+1]
      
      parameters <- c(T=9,L=1,r=1,e=0.01,mu=0.5,D1=m[j]*r[i],D2=m[j],F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn9[i,j,z]<-out[Timesteps+1,z+1]
    }
  }
}

###before density KF=100, mu=0.5, and F=0.25
den_bef<-100*(1-0.5-0.25)

##local effect
eqnmean1<-array(NA,dim=c(length(r),length(m),2))
eqnmean5<-array(NA,dim=c(length(r),length(m),2))
eqnmean9<-array(NA,dim=c(length(r),length(m),2))

loceff1.1<-array(NA,dim=c(length(r),length(m)))
loceff5.1<-array(NA,dim=c(length(r),length(m)))
loceff9.1<-array(NA,dim=c(length(r),length(m)))

for(i in 1:length(r))
{
  for(j in 1:length(m))
  {
    eqnmean1[i,j,1]<-mean(eqn1[i,j,1:9])
    eqnmean1[i,j,2]<-mean(eqn1[i,j,10])
    loceff1.1[i,j]<-eqnmean1[i,j,2]/eqnmean1[i,j,1]
    
    eqnmean5[i,j,1]<-mean(eqn5[i,j,1:5])
    eqnmean5[i,j,2]<-mean(eqn5[i,j,6:10])
    loceff5.1[i,j]<-eqnmean5[i,j,2]/eqnmean5[i,j,1]
    
    eqnmean9[i,j,1]<-mean(eqn9[i,j,1])
    eqnmean9[i,j,2]<-mean(eqn9[i,j,2:10])
    loceff9.1[i,j]<-eqnmean9[i,j,2]/eqnmean9[i,j,1]
    
  }
} 

###combine three datasets at D_M=0.5
bind1<-rbind(eqn1[,1,],eqn5[,1,],eqn9[,1,])
bind<-rbind(eqn1[,1,],eqn5[,1,],eqn9[,1,])/den_bef


###local effect calculation
L1high<-mean(eqn1[1,10,10])/mean(eqn1[1,10,1:9])
L2high<-mean(eqn1[2,10,10])/mean(eqn1[2,10,1:9])
L3high<-mean(eqn1[3,10,10])/mean(eqn1[3,10,1:9])

####local effect calculation
L1media<-mean(eqn5[1,10,6:10])/mean(eqn5[1,10,1:5])
L2media<-mean(eqn5[2,10,6:10])/mean(eqn5[2,10,1:5])
L3media<-mean(eqn5[3,10,6:10])/mean(eqn5[3,10,1:5])

##local effect calculation
L1small<-mean(eqn9[1,10,2:10])/mean(eqn9[1,10,1])
L2small<-mean(eqn9[2,10,2:10])/mean(eqn9[2,10,1])
L3small<-mean(eqn9[3,10,2:10])/mean(eqn9[3,10,1])


tiff("Semi_close_Fig.S3.tiff", width=5,height=5, units='in',res=600)
par(mfrow=c(3,4))
par(mar=c(2,2,1.5,1.6))

###at small MPA size SM=1 
yval6<-eqn9[1,1,]/den_bef
plot(seq(1,10,1),eqn9[1,1,]/den_bef,lwd=2,xlab="",ylab="",ylim=c(min(bind),max(bind)),col="#31A9B8",pch=16)
points(seq(1,10,1),yval6,lwd=1,type="p",col="black")
points(1,yval6[1],col="#CF3721",pch=16)
abline(v=1.5,lwd=2,col="black",lty=5)
text(4,1.5,paste("LE=",round(mean(eqn9[1,1,2:10])/mean(eqn9[1,1,1]),digits=2)),cex=0.8)

yval7<-eqn9[2,1,]/den_bef
plot(seq(1,10,1),eqn9[2,1,]/den_bef,lwd=2,type="p",lty=2,ylim=c(min(bind),max(bind)),col="#31A9B8",pch=16)
points(seq(1,10,1),yval7,lwd=1,type="p",col="black")
points(1,yval7[1],col="#CF3721",pch=16)
abline(v=1.5,lwd=2,col="black",lty=5)
text(4,1.5,paste("LE=",round(mean(eqn9[2,1,2:10])/mean(eqn9[2,1,1]),digits=2)),cex=0.8)

yval8<-eqn9[3,1,]/den_bef
plot(seq(1,10,1),eqn9[3,1,]/den_bef,lwd=2,type="p",lty=3,ylim=c(min(bind),max(bind)),col="#31A9B8",pch=16)
points(seq(1,10,1),yval8,lwd=1,type="p",col="black")
points(1,yval8[1],col="#CF3721",pch=16)
abline(v=1.5,lwd=2,col="black",lty=5)
text(4,1.5,paste("LE=",round(mean(eqn9[3,1,2:10])/mean(eqn9[3,1,1]),digits=2)),cex=0.8)

yval80<-eqn9[4,1,]/den_bef
plot(seq(1,10,1),eqn9[4,1,]/den_bef,lwd=2,type="p",lty=3,ylim=c(min(bind),max(bind)),col="#31A9B8",pch=16)
points(seq(1,10,1),yval80,lwd=1,type="p",col="black")
points(1,yval80[1],col="#CF3721",pch=16)
abline(v=1.5,lwd=2,col="black",lty=5)
text(4,1.5,paste("LE=",round(mean(eqn9[4,1,2:10])/mean(eqn9[4,1,1]),digits=2)),cex=0.8)

###at medium MPA size SM=5 
yval3<-eqn5[1,1,]/den_bef
plot(seq(1,10,1),yval3,lwd=2,type="p",pch=16,xlab="",ylab="",ylim=c(min(bind),max(bind)),col="#31A9B8")
points(seq(1,10,1),eqn5[1,1,]/den_bef,lwd=1,type="p",col="black")
points(1:5,yval3[1:5],type="p",pch=16,lwd=2,col="#CF3721")
abline(v=5.5,lwd=2,col="black",lty=5)
text(3.1,0.9,paste("LE=",round(mean(eqn5[1,1,6:10])/mean(eqn5[1,1,1:5]),digits=2)),cex=0.8)

yval4<-eqn5[2,1,]/den_bef
plot(seq(1,10,1),eqn5[2,1,]/den_bef,lwd=2,type="p",lty=2,ylim=c(min(bind),max(bind)),pch=16,col="#31A9B8")
points(seq(1,10,1),yval4,lwd=1,type="p",col="black")
points(1:5,yval4[1:5],pch=16,lwd=2,col="#CF3721")
abline(v=5.5,lwd=2,col="black",lty=5)
text(3.1,0.9,paste("LE=",round(mean(eqn5[2,1,6:10])/mean(eqn5[2,1,1:5]),digits=2)),cex=0.8)

yval5<-eqn5[3,1,]/den_bef
plot(seq(1,10,1),eqn5[3,1,]/den_bef,lwd=2,type="p",lty=3,ylim=c(min(bind),max(bind)),col="#31A9B8",pch=16)
points(seq(1,10,1),yval5,lwd=1,type="p",col="black")
points(1:5,yval5[1:5],lwd=2,pch=16,col="#CF3721")
abline(v=5.5,lwd=2,col="black",lty=5)
text(3.1,0.9,paste("LE=",round(mean(eqn5[3,1,6:10])/mean(eqn5[3,1,1:5]),digits=2)),cex=0.8)

yval50<-eqn5[4,1,]/den_bef
plot(seq(1,10,1),eqn5[4,1,]/den_bef,lwd=2,type="p",lty=3,ylim=c(min(bind),max(bind)),col="#31A9B8",pch=16)
points(seq(1,10,1),yval50,lwd=1,type="p",col="black")
points(1:5,yval50[1:5],lwd=2,pch=16,col="#CF3721")
abline(v=5.5,lwd=2,col="black",lty=5)
text(3.1,0.9,paste("LE=",round(mean(eqn5[4,1,6:10])/mean(eqn5[4,1,1:5]),digits=2)),cex=0.8)

###at large MPA size SM=9 
yval<-eqn1[1,1,]/den_bef
plot(seq(1,10,1),yval,lwd=2,xlab="",ylab="",ylim=c(min(bind),max(bind)),pch=16,col="#31A9B8")
points(seq(1,10,1),eqn1[1,1,]/den_bef,lwd=1,type="p",col="black")
points(1:9,yval[1:9],type="p",pch=16,lwd=2,col="#CF3721")
abline(v=9.5,lwd=2,col="black",lty=5)
text(3.5,1.5,paste("LE=", round(mean(eqn1[1,1,10])/mean(eqn1[1,1,1:9]),digits=2)),cex=0.8)
#plot(seq(1,10,1),eqn1[2,10,]/mean(eqn5[2,10,]),lwd=2,type="p",lty=2,ylim=c(min(bind),max(bind)))
yval1<-eqn1[2,1,]/den_bef
plot(1:10,yval1,lwd=2,type="p",lty=2,pch=16,ylim=c(min(bind),max(bind)),col="#31A9B8")
points(seq(1,10,1),eqn1[2,1,]/den_bef,lwd=1,type="p",col="black")
points(1:9,yval1[1:9],lwd=2,type="p",pch=16,col="#CF3721")
abline(v=9.5,lwd=2,col="black",lty=5)
text(3.5,1.5,paste("LE=",round(mean(eqn1[2,1,10])/mean(eqn1[2,1,1:9]),digits=2)),cex=0.8)

yval2<-eqn1[3,1,]/den_bef
plot(seq(1,10,1),yval2,lwd=2,type="p",pch=16,lty=3,ylim=c(min(bind),max(bind)),col="#31A9B8")
points(seq(1,10,1),eqn1[3,1,]/den_bef,lwd=1,type="p",col="black")
points(1:9,yval2[1:9],lwd=2,type="p",pch=16,col="#CF3721")
abline(v=9.5,lwd=2,col="black",lty=5)
text(3.5,1.5,paste("LE=",round(mean(eqn1[3,1,10])/mean(eqn1[3,1,1:9]),digits=2)),cex=0.8)

yval20<-eqn1[4,1,]/den_bef
plot(seq(1,10,1),yval20,lwd=2,type="p",pch=16,lty=3,ylim=c(min(bind),max(bind)),col="#31A9B8")
points(seq(1,10,1),eqn1[4,1,]/den_bef,lwd=1,type="p",col="black")
points(1:9,yval20[1:9],lwd=2,type="p",pch=16,col="#CF3721")
abline(v=9.5,lwd=2,col="black",lty=5)
text(3.5,1.5,paste("LE=",round(mean(eqn1[4,1,10])/mean(eqn1[4,1,1:9]),digits=2)),cex=0.8)

dev.off()

