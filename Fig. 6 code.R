##take off exponential term +1 and interference competition 
setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

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
    
    if(L==1)
    {
      A[1]<-R-(mu+F)*x[1]-D1*x[1]+D2/2*y[1]+D2/2*y[T]
    }
    else
    {
      A[1]<-R-(mu+F)*x[1]-D1*x[1]+D1/2*x[2]+D2/2*y[1]
      
      A[L]<-R-(mu+F)*x[L]-D1*x[L]+D1/2*x[L-1]+D2/2*y[T]
      
    }
    if(L-1>=2)
    {
      for(i in 2:(L-1))
      {
        A[i]<-R-(mu+F)*x[i]-D1*x[i]+D1/2*(x[i-1]+x[i+1])
      }
    }
    if(T==1)
    {
      A[(L+1)]<-R-mu*y[1]-D2*y[1]+D1/2*x[L]+D1/2*x[1]	
    }
    else
    {
      A[(L+1)]<-R-mu*y[1]-D2*y[1]+D2/2*y[2]+D1/2*x[1]
      A[(T+L)]<-R-mu*y[T]-D2*y[T]+D2/2*y[T-1]+D1/2*x[L]
    }
    if(T-1>=2)
    {
      for(i in (L+2):(T+L-1))
      {
        A[i]<-R-mu*y[i-L]-D2*y[i-L]+D2/2*(y[i-L-1]+y[i-L+1])
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


r<-c(1,1.5,2)

m<-seq(0,10,0.5)

F=0.25

##h is for MPA size, r is for differential movement
eqn1<-array(NA,dim=c(length(r),length(m),S))
eqn5<-array(NA,dim=c(length(r),length(m),S))
eqn9<-array(NA,dim=c(length(r),length(m),S))
eqn_before<-array(NA,dim=c(length(r),length(m),S))


for(i in 1:length(r))
{
  for(j in 1:length(m))
  {
    for(z in 1:S)
    {
      parameters <- c(T=1,L=9,R=2,mu=0.5,D1=m[j]*r[i],D2=m[j],F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn1[i,j,z]<-out[Timesteps+1,z+1]
      
      parameters <- c(T=5,L=5,R=2,mu=0.5,D1=m[j]*r[i],D2=m[j],F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn5[i,j,z]<-out[Timesteps+1,z+1]
      
      parameters <- c(T=9,L=1,R=2,mu=0.5,D1=m[j]*r[i],D2=m[j],F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn9[i,j,z]<-out[Timesteps+1,z+1]
      
      parameters_before <- c(T=5,L=5,R=2,mu=0.5,D1=0,D2=0,F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn_before[i,j,z]<-out[Timesteps+1,z+1]
    }
  }
}

###before MPA
###cell density
den_bef<-eqn_before[1,1,1]
###local effect: =1 since there is no MPA
loc_before=1
###regional abundance:eqnmean_ba[1,1,1]=2.666667
reg_before<-eqn_before[1,1,1]*10
###fishing yield
fis_before<-F*reg_before

bind1<-rbind(eqn5[1,10,]/den_bef,eqn5[2,10,]/den_bef,eqn5[3,10,]/den_bef)
bind<-bind1

tiff("Open_system_Fig.6.tiff", width=5,height=5, units='in',res=600)
par(mfrow=c(1,3))
par(mar=c(12,2,12,2))
###at large MPA size SM=9 
#plot(seq(1,10,1),eqn1[1,10,],lwd=2,xlab="",ylab="",ylim=c(min(bind),max(bind)))
#abline(v=9.5,lwd=2,col="red")
#text(3.5,4.5,paste("LE=", round(mean(eqn1[1,10,10])/mean(eqn1[1,10,1:9]),digits=2)))
yvalue<-eqn5[1,10,]/den_bef
plot(seq(1,10,1),eqn5[1,10,]/den_bef,lwd=2,type="p",lty=2,ylim=c(min(bind)-0.03,max(bind)+0.03),xlab="",ylab="",pch=16,col="#31A9B8")
points(seq(1,10,1),eqn5[1,10,]/den_bef,lwd=1,type="p",col="black")
points(1:5,yvalue[1:5],col="red",pch=16,type="p")
abline(v=5.5,lwd=2,lty=5)
text(3.3,1.4,paste("LE=",round(mean(eqn5[1,10,6:10])/mean(eqn5[1,10,1:5]),digits=2)))

yvalue1<-eqn5[2,10,]/den_bef
plot(seq(1,10,1),eqn5[2,10,]/den_bef,lwd=2,type="p",lty=2,ylim=c(min(bind)-0.03,max(bind)+0.03),xlab="",ylab="",pch=16,col="#31A9B8")
points(seq(1,10,1),eqn5[2,10,]/den_bef,lwd=1,type="p",col="black")
points(1:5,yvalue1[1:5],col="red",pch=16,type="p")
abline(v=5.5,lwd=2,lty=5)
text(3.3,1.4,paste("LE=",round(mean(eqn5[2,10,6:10])/mean(eqn5[2,10,1:5]),digits=2)))

yvalue2<-eqn5[3,10,]/den_bef
plot(seq(1,10,1),eqn5[3,10,]/den_bef,lwd=2,type="p",lty=2,ylim=c(min(bind)-0.03,max(bind)+0.03),xlab="",ylab="",pch=16,col="#31A9B8")
points(seq(1,10,1),eqn5[3,10,]/den_bef,lwd=1,type="p",col="black")
points(1:5,yvalue2[1:5],col="red",pch=16,type="p")
abline(v=5.5,lwd=2,lty=5)
text(3.3,1.4,paste("LE=",round(mean(eqn5[3,10,6:10])/mean(eqn5[3,10,1:5]),digits=2)))
dev.off()

