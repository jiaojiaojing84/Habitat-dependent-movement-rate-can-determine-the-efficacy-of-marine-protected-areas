#updated by JJ on 6-7 on heatmap direction and take off cell number
####mu_F=0.1
rm(list=ls())

#setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape2)



library(ggplot2)



library(scales)



library(pheatmap)



library(tidyverse)



################test
#####L is fishing ground while T is MPAs, total area is S
JAP08<-function(t, inits,parameters) {
  
  
  with(as.list(c(inits, parameters)),{
    
    x<-inits[1:L]
    y<-inits[(L+1):(L+T)]
    
    A<-array(NA,dim=c(1,(L+T)))
    
    if(L==1)
    {
      A[1]<-r*x[1]-e*x[1]^2-(mu+F)*x[1]-D1*x[1]+D2/2*y[1]+D2/2*y[T]
    }
    else
    {
      A[1]<-r*x[1]-e*x[1]^2-(mu+F)*x[1]-D1*x[1]+D1/2*x[2]+D2/2*y[1]
      
      A[L]<-r*x[1]-e*x[1]^2-(mu+F)*x[L]-D1*x[L]+D1/2*x[L-1]+D2/2*y[T]
      
    }
    if(L-1>=2)
    {
      for(i in 2:(L-1))
      {
        A[i]<-r*x[i]-e*x[i]^2-(mu+F)*x[i]-D1*x[i]+D1/2*(x[i-1]+x[i+1])
      }
    }
    if(T==1)
    {
      A[(L+1)]<-r1*y[1]-e1*y[1]^2-mu*y[1]-D2*y[1]+D1/2*x[L]+D1/2*x[1]	
    }
    else
    {
      A[(L+1)]<-r1*y[1]-e1*y[1]^2-mu*y[1]-D2*y[1]+D2/2*y[2]+D1/2*x[1]
      A[(T+L)]<-r1*y[T]-e1*y[T]^2-mu*y[T]-D2*y[T]+D2/2*y[T-1]+D1/2*x[L]
    }
    if(T-1>=2)
    {
      for(i in (L+2):(T+L-1))
      {
        A[i]<-r1*y[i-L]-e1*y[i-L]^2-mu*y[i-L]-D2*y[i-L]+D2/2*(y[i-L-1]+y[i-L+1])
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


beta<-seq(1,8,0.5)

h<-seq(1,9,1)

#F=seq(0.01,2,0.2)
F=0.1

##h is for MPA size, r is for differential movement
#eqn<-array(NA,dim=c(length(beta),length(h),length(F),S))

eqn<-array(NA,dim=c(length(beta),length(h),S))


for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    
    for(z in 1:S)
    {
      
      parameters <- c(T=h[j],L=S-h[j],r=1,e=0.01,e1=0.01,r1=1,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.1) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn[i,j,z]<-out[Timesteps+1,z+1]
    }
  }
}

###before status
###local effect: =1 since there is no MPA
loc_before1=1
###regional abundance:
reg_before1<-40*10
#reg_before1<-1
###fishing yield
fis_before1<-F*40*10
#fis_before1<-1

#local effect
yield<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    yield[i,j]<-sum(eqn[i,j,1:(S-h[j])])*F
  }
} 
##local effect
eqnmean<-array(NA,dim=c(length(beta),length(h),2))
loceff<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    eqnmean[i,j,1]<-mean(eqn[i,j,1:(S-h[j])])
    eqnmean[i,j,2]<-mean(eqn[i,j,(S-h[j]+1):S])
    loceff[i,j]<-eqnmean[i,j,2]/eqnmean[i,j,1]
    
  }
} 
##regional effect
regeff<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    regeff[i,j]<-sum(eqn[i,j,1:S]) 
  }
} 

yield1<-melt(yield)
loceff1<-melt(loceff)
regeff1<-melt(regeff)



theme_set(theme_bw(20))


yield2<-yield1
yield2$value<-yield1$value/fis_before1

loceff2<-loceff1
loceff2$value<-loceff1$value/loc_before1

regeff2<-regeff1
regeff2$value<-regeff1$value/reg_before1

tiff("Fig.S6_yield_close_F_0.1.tiff", width=5,height=5, units='in',res=600)


p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) )+guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S6_local-effect_close_F_0.1.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S6_regional_close_F_0.1.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(regeff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


#####mu_F=0.25
rm(list=ls())

#setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape2)



library(ggplot2)



library(scales)



library(pheatmap)



library(tidyverse)



################test
#####L is fishing ground while T is MPAs, total area is S
JAP08<-function(t, inits,parameters) {
  
  
  with(as.list(c(inits, parameters)),{
    
    x<-inits[1:L]
    y<-inits[(L+1):(L+T)]
    
    A<-array(NA,dim=c(1,(L+T)))
    
    if(L==1)
    {
      A[1]<-r*x[1]-e*x[1]^2-(mu+F)*x[1]-D1*x[1]+D2/2*y[1]+D2/2*y[T]
    }
    else
    {
      A[1]<-r*x[1]-e*x[1]^2-(mu+F)*x[1]-D1*x[1]+D1/2*x[2]+D2/2*y[1]
      
      A[L]<-r*x[1]-e*x[1]^2-(mu+F)*x[L]-D1*x[L]+D1/2*x[L-1]+D2/2*y[T]
      
    }
    if(L-1>=2)
    {
      for(i in 2:(L-1))
      {
        A[i]<-r*x[i]-e*x[i]^2-(mu+F)*x[i]-D1*x[i]+D1/2*(x[i-1]+x[i+1])
      }
    }
    if(T==1)
    {
      A[(L+1)]<-r1*y[1]-e1*y[1]^2-mu*y[1]-D2*y[1]+D1/2*x[L]+D1/2*x[1]	
    }
    else
    {
      A[(L+1)]<-r1*y[1]-e1*y[1]^2-mu*y[1]-D2*y[1]+D2/2*y[2]+D1/2*x[1]
      A[(T+L)]<-r1*y[T]-e1*y[T]^2-mu*y[T]-D2*y[T]+D2/2*y[T-1]+D1/2*x[L]
    }
    if(T-1>=2)
    {
      for(i in (L+2):(T+L-1))
      {
        A[i]<-r1*y[i-L]-e1*y[i-L]^2-mu*y[i-L]-D2*y[i-L]+D2/2*(y[i-L-1]+y[i-L+1])
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


beta<-seq(1,8,0.5)

h<-seq(1,9,1)

#F=seq(0.01,2,0.2)
F=0.25

##h is for MPA size, r is for differential movement
#eqn<-array(NA,dim=c(length(beta),length(h),length(F),S))

eqn<-array(NA,dim=c(length(beta),length(h),S))


for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    
    for(z in 1:S)
    {
      
      parameters <- c(T=h[j],L=S-h[j],r=1,e=0.01,e1=0.01,r1=1,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn[i,j,z]<-out[Timesteps+1,z+1]
    }
  }
}

###before status
###local effect: =1 since there is no MPA
loc_before1=1
###regional abundance:
reg_before1<-25*10
#reg_before1<-1
###fishing yield
fis_before1<-F*25*10
#fis_before1<-1

#local effect
yield<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    yield[i,j]<-sum(eqn[i,j,1:(S-h[j])])*F
  }
} 
##local effect
eqnmean<-array(NA,dim=c(length(beta),length(h),2))
loceff<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    eqnmean[i,j,1]<-mean(eqn[i,j,1:(S-h[j])])
    eqnmean[i,j,2]<-mean(eqn[i,j,(S-h[j]+1):S])
    loceff[i,j]<-eqnmean[i,j,2]/eqnmean[i,j,1]
    
  }
} 
##regional effect
regeff<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    regeff[i,j]<-sum(eqn[i,j,1:S]) 
  }
} 

yield1<-melt(yield)
loceff1<-melt(loceff)
regeff1<-melt(regeff)


##at 50% MPA and FG
#rho<-loceff[,5]
#plot(rho,lwd=2,xlim=c(0,100),ylim=c(0,100),ylab="rho",xlab="beta",type="l")
#y<-function(x) {y=x}
#curve(y,0,100,lwd=2,add=T,col="red",type="l")
#y1<-function(x1) {y1=0.5*x1}
#curve(y1,0,100,lwd=2,add=T,col="red",lty=2)

theme_set(theme_bw(20))

#tiff("density in FG and MPA 5-2.tiff", width=5,height=5, units='in',res=600)
#par(mfrow=c(2,1))
#par(mar=c(1.98,4,1,0.8))
#plot(eqnmean[,5,1],type="l",lwd=2,xlab="beta",ylab="density in FG")
#plot(eqnmean[,5,2],type="l",lwd=2,xlab="beta",ylab="density in MPA")
#abline(h=100,col="red",lwd=2)
#dev.off()
yield2<-yield1
yield2$value<-yield1$value/fis_before1

loceff2<-loceff1
loceff2$value<-loceff1$value/loc_before1

regeff2<-regeff1
regeff2$value<-regeff1$value/reg_before1

tiff("Fig.S6_yield_close_F_0.25.tiff", width=5,height=5, units='in',res=600)

#p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) + geom_text(aes(label=round(value,2)),position = position_dodge(width=0.1),  size=1.8)+guides(fill = guide_colorbar(title=""))

p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) )+guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S6_local-effect_close_F_0.25.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S6_regional_close_F_0.25.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(regeff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()




####mu_F=0.45
rm(list=ls())

#setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape2)



library(ggplot2)



library(scales)



library(pheatmap)



library(tidyverse)



################test
#####L is fishing ground while T is MPAs, total area is S
JAP08<-function(t, inits,parameters) {
  
  
  with(as.list(c(inits, parameters)),{
    
    x<-inits[1:L]
    y<-inits[(L+1):(L+T)]
    
    A<-array(NA,dim=c(1,(L+T)))
    
    if(L==1)
    {
      A[1]<-r*x[1]-e*x[1]^2-(mu+F)*x[1]-D1*x[1]+D2/2*y[1]+D2/2*y[T]
    }
    else
    {
      A[1]<-r*x[1]-e*x[1]^2-(mu+F)*x[1]-D1*x[1]+D1/2*x[2]+D2/2*y[1]
      
      A[L]<-r*x[1]-e*x[1]^2-(mu+F)*x[L]-D1*x[L]+D1/2*x[L-1]+D2/2*y[T]
      
    }
    if(L-1>=2)
    {
      for(i in 2:(L-1))
      {
        A[i]<-r*x[i]-e*x[i]^2-(mu+F)*x[i]-D1*x[i]+D1/2*(x[i-1]+x[i+1])
      }
    }
    if(T==1)
    {
      A[(L+1)]<-r1*y[1]-e1*y[1]^2-mu*y[1]-D2*y[1]+D1/2*x[L]+D1/2*x[1]	
    }
    else
    {
      A[(L+1)]<-r1*y[1]-e1*y[1]^2-mu*y[1]-D2*y[1]+D2/2*y[2]+D1/2*x[1]
      A[(T+L)]<-r1*y[T]-e1*y[T]^2-mu*y[T]-D2*y[T]+D2/2*y[T-1]+D1/2*x[L]
    }
    if(T-1>=2)
    {
      for(i in (L+2):(T+L-1))
      {
        A[i]<-r1*y[i-L]-e1*y[i-L]^2-mu*y[i-L]-D2*y[i-L]+D2/2*(y[i-L-1]+y[i-L+1])
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


beta<-seq(1,8,0.5)

h<-seq(1,9,1)

#F=seq(0.01,2,0.2)
F=0.45

##h is for MPA size, r is for differential movement
#eqn<-array(NA,dim=c(length(beta),length(h),length(F),S))

eqn<-array(NA,dim=c(length(beta),length(h),S))


for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    
    for(z in 1:S)
    {
      
      parameters <- c(T=h[j],L=S-h[j],r=1,e=0.01,e1=0.01,r1=1,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.45) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn[i,j,z]<-out[Timesteps+1,z+1]
    }
  }
}

###before status
###local effect: =1 since there is no MPA
loc_before1=1
###regional abundance:
reg_before1<-5*10
#reg_before1<-1
###fishing yield
fis_before1<-F*5*10
#fis_before1<-1

#local effect
yield<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    yield[i,j]<-sum(eqn[i,j,1:(S-h[j])])*F
  }
} 
##local effect
eqnmean<-array(NA,dim=c(length(beta),length(h),2))
loceff<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    eqnmean[i,j,1]<-mean(eqn[i,j,1:(S-h[j])])
    eqnmean[i,j,2]<-mean(eqn[i,j,(S-h[j]+1):S])
    loceff[i,j]<-eqnmean[i,j,2]/eqnmean[i,j,1]
    
  }
} 
##regional effect
regeff<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    regeff[i,j]<-sum(eqn[i,j,1:S]) 
  }
} 

yield1<-melt(yield)
loceff1<-melt(loceff)
regeff1<-melt(regeff)


##at 50% MPA and FG
#rho<-loceff[,5]
#plot(rho,lwd=2,xlim=c(0,100),ylim=c(0,100),ylab="rho",xlab="beta",type="l")
#y<-function(x) {y=x}
#curve(y,0,100,lwd=2,add=T,col="red",type="l")
#y1<-function(x1) {y1=0.5*x1}
#curve(y1,0,100,lwd=2,add=T,col="red",lty=2)

theme_set(theme_bw(20))

#tiff("density in FG and MPA 5-2.tiff", width=5,height=5, units='in',res=600)
#par(mfrow=c(2,1))
#par(mar=c(1.98,4,1,0.8))
#plot(eqnmean[,5,1],type="l",lwd=2,xlab="beta",ylab="density in FG")
#plot(eqnmean[,5,2],type="l",lwd=2,xlab="beta",ylab="density in MPA")
#abline(h=100,col="red",lwd=2)
#dev.off()
yield2<-yield1
yield2$value<-yield1$value/fis_before1

loceff2<-loceff1
loceff2$value<-loceff1$value/loc_before1

regeff2<-regeff1
regeff2$value<-regeff1$value/reg_before1

tiff("6-13-18-Fig.S6_yield_close_F_0.45.tiff", width=5,height=5, units='in',res=600)

#p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) + geom_text(aes(label=round(value,2)),position = position_dodge(width=0.1),  size=1.8)+guides(fill = guide_colorbar(title=""))

p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) )+guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("6-13-18-Fig.S6_local-effect_close_F_0.45.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("6-13-18-Fig.S6_regional_close_F_0.45.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(regeff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()




