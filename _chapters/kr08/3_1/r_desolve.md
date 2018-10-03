---
interact_link: notebooks/kr08/3_1/r_desolve.ipynb
title: 'R using deSolve'
permalink: 'chapters/kr08/3_1/r_desolve'
previouschapter:
  url: chapters/kr08/3_1/intro
  title: 'Program 3.1: SIS with risk groups'
nextchapter:
  url: chapters/kr08/3_1/julia
  title: 'Julia'
redirect_from:
  - 'chapters/kr08/3-1/r-desolve'
---

## Problem 3.1: SIS model with risk groups

Author: Emma Accorsi @emmaaccorsi

Date: 2018-10-01

Import libraries.


{:.input_area}
```R
library(deSolve)
library(ggplot2)
library(reshape2)
```

Specify SIS model function.


{:.input_area}
```R
sis_ode <- function(times,x,parms){
  with(as.list(c(parms,x)),{
    # ODEs
    SH<-nH-x[[1]] #Calculate SH as nH-IH
    SL<-(1-nH)-x[[2]] #Calculate SL as nL-IL
    dIH <-+(betaHH*IH+betaHL*IL)*SH-gamma*IH
    dIL <-+(betaLH*IH+betaLL*IL)*SL-gamma*IL
    
    der<-c(dIH,dIL)
    list(der)
  })
}
```

Specify parameter values and run SIS model


{:.input_area}
```R
parms <- c(betaHH=10,betaHL=0.1,betaLH=0.1,betaLL=1,gamma=1,nH=0.2)
x <- c(IH=0.00001,IL=0.001)
times <-seq(0,15,1)
sis_out <- as.data.frame(lsoda(x,times,sis_ode,parms))
```

Create visualization with ggplot2


{:.input_area}
```R
sis_out_long <- melt(sis_out,"time") #Collapse dataset from "wide" to "long" format for plotting
ggplot(sis_out_long,aes(x=time,y=value,colour=variable,group=variable))+
  # Add line
  geom_line(lwd=2)+
  #Add labels
  labs(x="Time (Years)",y="Proportion of Population",color="Risk Group")
```




![png](../../../images/chapters/kr08/3_1/r_desolve_8_1.png)

