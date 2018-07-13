###For mu_F=0.1
rm(list=ls())

#setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape)



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


beta<-seq(1,4.1,0.25)

h<-seq(1,9,1)

##F here is mu_F in the model
F=0.1

##h is for MPA size, r is for differential movement
eqn<-array(NA,dim=c(length(beta),length(h),S))
eqn_before<-array(NA,dim=c(length(beta),length(h),S))



for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    for(z in 1:S)
    {
      parameters <- c(T=h[j],L=S-h[j],R=2,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.1) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn[i,j,z]<-out[Timesteps+1,z+1]
      
      ###before data for each cell
      parameters_before <- c(T=h[j],L=S-h[j],R=2,mu=0.5,D1=0*beta[i],D2=0,F=0.1) 
      out_before= ode(y = inits, times = times, func = JAP08, parms = parameters_before)
      eqn_before[i,j,z]<-out_before[Timesteps+1,z+1]
      
    }
  }
}
###before-after
###using after density at each cell in MPA / before density at each cell

####after mean/before mean indicates the before-after effect
eqnmean_ba<-array(NA,dim=c(length(beta),length(h),2))
before_after<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    eqnmean_ba[i,j,1]<-mean(eqn_before[i,j,1:(S-h[j])]) # mean density before MPA
    eqnmean_ba[i,j,2]<-mean(eqn[i,j,(S-h[j]+1):S]) # mean density in MPA 
    before_after[i,j]<-eqnmean_ba[i,j,2]/eqnmean_ba[i,j,1]
    
    
  }
} 

###before status
###local effect: =1 since there is no MPA
loc_before=1
###regional abundance:eqnmean_ba[1,1,1]=2.666667
reg_before<-eqnmean_ba[1,1,1]*10
###fishing yield
fis_before<-F*reg_before

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

before_after1<-melt(before_after)
yield1<-melt(yield)
loceff1<-melt(loceff)
regeff1<-melt(regeff)


theme_set(theme_bw(20))


before_after2<-before_after1
#before_after2$value<-before_after1$value/before_after1$value[61]
before_after2$value<-before_after1$value/1

yield2<-yield1
#yield2$value<-yield1$value/yield1$value[61]
yield2$value<-yield1$value/fis_before


loceff2<-loceff1
#loceff2$value<-loceff1$value/loceff1$value[61]
loceff2$value<-loceff1$value/1

regeff2<-regeff1
#regeff2$value<-regeff1$value/regeff1$value[61]
regeff2$value<-regeff1$value/reg_before



tiff("Fig.S5_fishing yield-0.1.tiff", width=5,height=5, units='in',res=600)


p1 <- ggplot(yield2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value))+guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S5_local-effect-0.1.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value))+guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S5_regional-effect-0.1.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(regeff2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value)) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue',
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))
dev.off()



###For mu_F=0.25
rm(list=ls())

#setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape)



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


beta<-seq(1,4.1,0.25)

h<-seq(1,9,1)

##F here is mu_F in the model
F=0.25

##h is for MPA size, r is for differential movement
eqn<-array(NA,dim=c(length(beta),length(h),S))
eqn_before<-array(NA,dim=c(length(beta),length(h),S))



for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
  	for(z in 1:S)
  	{
  	parameters <- c(T=h[j],L=S-h[j],R=2,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.25) 
  	out= ode(y = inits, times = times, func = JAP08, parms = parameters)
  	eqn[i,j,z]<-out[Timesteps+1,z+1]
  	
  	###before data for each cell
  	parameters_before <- c(T=h[j],L=S-h[j],R=2,mu=0.5,D1=0*beta[i],D2=0,F=0.25) 
  	out_before= ode(y = inits, times = times, func = JAP08, parms = parameters_before)
  	eqn_before[i,j,z]<-out_before[Timesteps+1,z+1]
  	
     }
   }
}
###before-after
###using after density at each cell in MPA / before density at each cell

####after mean/before mean indicates the before-after effect
eqnmean_ba<-array(NA,dim=c(length(beta),length(h),2))
before_after<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    eqnmean_ba[i,j,1]<-mean(eqn_before[i,j,1:(S-h[j])]) # mean density before MPA
    eqnmean_ba[i,j,2]<-mean(eqn[i,j,(S-h[j]+1):S]) # mean density in MPA 
    before_after[i,j]<-eqnmean_ba[i,j,2]/eqnmean_ba[i,j,1]
    
    
  }
} 

###before status
###local effect: =1 since there is no MPA
loc_before=1
###regional abundance:eqnmean_ba[1,1,1]=2.666667
reg_before<-eqnmean_ba[1,1,1]*10
###fishing yield
fis_before<-F*reg_before

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

before_after1<-melt(before_after)
yield1<-melt(yield)
loceff1<-melt(loceff)
regeff1<-melt(regeff)


theme_set(theme_bw(20))


before_after2<-before_after1
#before_after2$value<-before_after1$value/before_after1$value[61]
before_after2$value<-before_after1$value/1

yield2<-yield1
#yield2$value<-yield1$value/yield1$value[61]
yield2$value<-yield1$value/fis_before


loceff2<-loceff1
#loceff2$value<-loceff1$value/loceff1$value[61]
loceff2$value<-loceff1$value/1

regeff2<-regeff1
#regeff2$value<-regeff1$value/regeff1$value[61]
regeff2$value<-regeff1$value/reg_before



tiff("Fig.S5_fishing yield-0.25.tiff", width=5,height=5, units='in',res=600)


p1 <- ggplot(yield2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value))+guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S5_local-effect-0.25.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value))+guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S5_regional-effect-0.25.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(regeff2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value)) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue',
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))
dev.off()



####For mu_F=0.45
rm(list=ls())

#setwd("C:\\Users\\jiaojin1\\Downloads\\PhD work")

library(deSolve)



library(reshape)



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


beta<-seq(1,4.1,0.25)

h<-seq(1,9,1)

##F here is mu_F in the model
F=0.45

##h is for MPA size, r is for differential movement
eqn<-array(NA,dim=c(length(beta),length(h),S))
eqn_before<-array(NA,dim=c(length(beta),length(h),S))

for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    for(z in 1:S)
    {
      parameters <- c(T=h[j],L=S-h[j],R=2,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.45) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn[i,j,z]<-out[Timesteps+1,z+1]
      
      ###before data for each cell
      parameters_before <- c(T=h[j],L=S-h[j],R=2,mu=0.5,D1=0*beta[i],D2=0,F=0.45) 
      out_before= ode(y = inits, times = times, func = JAP08, parms = parameters_before)
      eqn_before[i,j,z]<-out_before[Timesteps+1,z+1]
      
    }
  }
}
###before-after
###using after density at each cell in MPA / before density at each cell

####after mean/before mean indicates the before-after effect
eqnmean_ba<-array(NA,dim=c(length(beta),length(h),2))
before_after<-array(NA,dim=c(length(beta),length(h)))
for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    eqnmean_ba[i,j,1]<-mean(eqn_before[i,j,1:(S-h[j])]) # mean density before MPA
    eqnmean_ba[i,j,2]<-mean(eqn[i,j,(S-h[j]+1):S]) # mean density in MPA 
    before_after[i,j]<-eqnmean_ba[i,j,2]/eqnmean_ba[i,j,1]
    
    
  }
} 

###before status
###local effect: =1 since there is no MPA
loc_before=1
###regional abundance:eqnmean_ba[1,1,1]=2.666667
reg_before<-eqnmean_ba[1,1,1]*10
###fishing yield
fis_before<-F*reg_before

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

before_after1<-melt(before_after)
yield1<-melt(yield)
loceff1<-melt(loceff)
regeff1<-melt(regeff)


theme_set(theme_bw(20))


before_after2<-before_after1
#before_after2$value<-before_after1$value/before_after1$value[61]
before_after2$value<-before_after1$value/1

yield2<-yield1
#yield2$value<-yield1$value/yield1$value[61]
yield2$value<-yield1$value/fis_before


loceff2<-loceff1
#loceff2$value<-loceff1$value/loceff1$value[61]
loceff2$value<-loceff1$value/1

regeff2<-regeff1
#regeff2$value<-regeff1$value/regeff1$value[61]
regeff2$value<-regeff1$value/reg_before



tiff("Fig.S5_fishing yield-0.45.tiff", width=5,height=5, units='in',res=600)


p1 <- ggplot(yield2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value))+guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S5_local-effect-0.45.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value))+guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S5_regional-effect-0.45.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(regeff2, aes(beta[X1], h[X2]/10)) + geom_tile(aes(fill = value)) +guides(fill = guide_colorbar(title=""))



p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue',
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))
dev.off()





