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


####For mu_F=0.1

#############

Timesteps=500

times <- seq(0, Timesteps, by = 1)

S=10

inits <- rep(1,S)


beta<-seq(1,8,0.5)

h<-seq(1,9,1)

##here F is mu_F
F=0.1


eqn<-array(NA,dim=c(length(beta),length(h),S))


for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    
    for(z in 1:S)
    {
      
      parameters <- c(T=h[j],L=S-h[j],r=1,e=0.01,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.1) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn[i,j,z]<-out[Timesteps+1,z+1]
      ###before data
      
    }
  }
}



###before status
###local effect: =1 since there is no MPA
loc_before1=1
###regional abundance: density in each cell* cell number (10 here)
##the density in each cell = K_F*(1-(mu_N+mu_F)/r_F)=100*(1-(0.5+0.1)/1)=40
reg_before1<-40*10

###fishing yield
fis_before1<-F*40*10


#yield
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


tiff("Fig.S7_yield_F_0.1.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S7_local-effect_F_0.1.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(high = 'red', mid='white',low = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S7_regional-effect_F_0.1.tiff", width=5,height=5, units='in',res=600)


p1 <- ggplot(regeff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) )+guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(high = 'red', mid='white',low = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))
dev.off()



#####For mu_F=0.25

rm(list=ls())

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



Timesteps=500

times <- seq(0, Timesteps, by = 1)

S=10

inits <- rep(1,S)


beta<-seq(1,8,0.5)

h<-seq(1,9,1)

###F here is mu_F
F=0.25

eqn<-array(NA,dim=c(length(beta),length(h),S))


for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    
    for(z in 1:S)
    {
      
      parameters <- c(T=h[j],L=S-h[j],r=1,e=0.01,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.25) 
      out= ode(y = inits, times = times, func = JAP08, parms = parameters)
      eqn[i,j,z]<-out[Timesteps+1,z+1]
      ###before data
      
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

#yield
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


tiff("Fig.S7_yield_F_0.25.tiff", width=5,height=5, units='in',res=600)


p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(high = 'red', mid='white',low = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S7_local-effect_F_0.25.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(high = 'red', mid='white',low = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S7_regional-effect_F_0.25.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(regeff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) )+guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(high = 'red', mid='white',low = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))
dev.off()




####for mu_F=0.45

rm(list=ls())


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


#beta<-seq(1,4.0,0.25)
beta<-seq(1,8,0.5)

h<-seq(1,9,1)


F=0.45

##h is for MPA size, r is for differential movement
#eqn<-array(NA,dim=c(length(beta),length(h),length(F),S))

eqn<-array(NA,dim=c(length(beta),length(h),S))
#eqn_before<-array(NA,dim=c(length(beta),length(h),S))


for(i in 1:length(beta))
{
  for(j in 1:length(h))
  {
    
  	  for(z in 1:S)
  	 {
  	  
  	parameters <- c(T=h[j],L=S-h[j],r=1,e=0.01,mu=0.5,D1=0.5*beta[i],D2=0.5,F=0.45) 
  	out= ode(y = inits, times = times, func = JAP08, parms = parameters)
  	eqn[i,j,z]<-out[Timesteps+1,z+1]
  	###before data
  	
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

#yield
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
#yield2$value<-yield1$value/yield1$value[61]
yield2$value<-yield1$value/fis_before1

loceff2<-loceff1
loceff2$value<-loceff1$value/loc_before1

regeff2<-regeff1
regeff2$value<-regeff1$value/reg_before1


tiff("Fig.S7_yield_F_0.45.tiff", width=5,height=5, units='in',res=600)

#p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) + geom_text(aes(label=round(value,1)),position = position_dodge(width=0.1),  size=1.8)+guides(fill = guide_colorbar(title=""))

p1 <- ggplot(yield2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S7_local-effect_F_0.45.tiff", width=5,height=5, units='in',res=600)

p1 <- ggplot(loceff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) +guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))

dev.off()


tiff("Fig.S7_regional-effect_F_0.45.tiff", width=5,height=5, units='in',res=600)

#p1 <- ggplot(regeff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) ) + geom_text(aes(label=round(value,1)),position = position_dodge(width=0.05),  size=1.8)+guides(fill = guide_colorbar(title=""))

p1 <- ggplot(regeff2, aes(beta[Var1], h[Var2]/10)) + geom_tile(aes(fill = value) )+guides(fill = guide_colorbar(title=""))


p1+xlab("")+ylab('')+ scale_fill_gradient2(low = 'red', mid='white',high = 'steelblue' ,
                                           
                                           midpoint=1, space = "rgb", na.value = "grey50", guide = 
                                             "colourbar")+scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand 
                                                                                                                    = c(0, 0))
dev.off()




