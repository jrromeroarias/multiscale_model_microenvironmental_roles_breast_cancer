## Microenvironmental factors in cell segregation and heterogeneity in breast cancer development
## By Roberto Romero, IIMAS UNAM CU MEXICO

# remove all
rm(list=ls())
# load lybraries
library(deSolve)
library(phaseR)
source('Grind.r')

# List of parameters
#***********************************************************************************************
alpha1 = -1.0;  beta1 = 9.0;  mu1 = 0.2;  v1 = 1.5;  gamma1 = 6;  K1 = 0.5; eta1 = 0.5;
alpha2 = -1.0;  beta2 = 0.2;  mu2 = 0.2;  v2 = 1;  gamma2 = 1;  K2 = 0.5; eta2 = 0.5;
alpha3 = -1.0;  beta3 = 0.2;  mu3 = 0.2;  v3 = 1;  gamma3 = 1;  K3 = 0.5; eta3 = 0.5;
alpha4 = -1.0;  beta4 = 0.2;  mu4 = 0.2;  v4 = 1;  gamma4 = 3;  K4 = 0.5; eta4 = 0.5;
alpha5 = 1.0;  beta5 = 20.0;  mu5 = 0.2;  v5 = 1;  gamma5 = 1;  K5 = 0.5; eta5 = 0.5;
alpha6 = 1.0;  beta6 = 2.0;  mu6 = 0.2;  v6 = 1.2;  gamma6 = 4;  K6 = 0.5; eta6 = 0.5;
alpha7 = -1.0;  beta7 = 0.0;  mu7 = 0.2;  v7 = 1;  gamma7 = 3;  K7 = 0.5; eta7 = 0.5;
alpha8 = 1.0;  beta8 = 0.2;  mu8 = 0.2;  v8 = 1;  gamma8 = 1;  K8 = 0.5; eta8 = 0.5;
alpha9 = 1.0;  beta9 = 3.0;  mu9 = 0.2;  v9 = 0.5;  gamma9 = 2;  K9 = 0.5; eta9 = 0.5;
alpha10 = 1.0; beta10 = 2.0; mu10 = 0.2; v10 = 1; gamma10 = 2; K10 = 0.5; eta10 = 0.5;

vy5 = 0.5; gammay5 = 1; Ky5 = 0.5;
lambda0 = 0.4; delta0 = 1;
lambda1 = 0.4; delta1 = 0.1;
#*************************************************************************************************
#input stimuli
y1 = 0.5; y2 = 0.5;
y6 = y1; y7 = y2;
y3 = y2; y8 = y1;
y4 = y2; y9 = y1;
y5 = y2; y10 = y2;

#initial values
x1 = 0; x6 = 0;
x2 = 0; x7 = 0;
x3 = 0; x8 = 0;
x4 = 0; x9 = 0;
x5 = 0; x10 = 0;

# initial conditions an time evolution
tspan <- seq(from = 0, to = 10, by = 0.01)
stepm = 0.001;


#******************************************************
# Dynamical system
#******************************************************
# for bifurcations with Grind.r
model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
   
    dx1  = (  alpha1*(y1/eta1 - 1)*(1 - x1)  + beta2*x1*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  + beta4*x1*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4)) - beta5*x1*((v5*x5)^gamma5/(K5 + (v5*x5)^gamma5)) - beta6*x1*((v6*x6)^gamma6/(K6 + (v6*x6)^gamma6)) +  beta8*x1*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta9*x1*((v9*x9)^gamma9/(K9 + (v9*x9)^gamma9)) - mu1*x1);
    dx2  = (  alpha2*(y2/eta2 - 1)*(1 - x2)  + beta4*x4*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4))  - mu2*x2);
    dx3  = (  alpha3*(y3/eta3 - 1)*(1 - x3)  + beta1*x3*((v1*x1)^gamma1/(K1 + (v1*x1)^gamma1))  - mu3*x3);
    dx4  = (  alpha4*(y4/eta4 - 1)*(1 - x4)  + beta2*x4*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  - beta8*x4*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta10*x4*((v10*x10)^gamma10/(K10 + (v10*x10)^gamma10)) - mu4*x4);
    dx5  = (  alpha5*(y5/eta5 - 1)*(1 - x5)  + lambda0*(x5/(delta0 + x5))*((vy5*y5)^gammay5/(Ky5     + (vy5*y5)^gammay5)) - mu5*x5);
    dx6  = (  alpha6*(y6/eta6 - 1)*(1 - x6)  + beta1*x6*((v1*x1)^gamma1/(K1 + (v1*x1)^gamma1))  + beta3*x6*((v3*x3)^gamma3/(K3 + (v3*x3)^gamma3)) + beta8*x6*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta2*x6*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2)) - mu6*x6);
    dx7  = (  alpha7*(y7/eta7 - 1)*(1 - x7)  + beta2*x7*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  + beta4*x7*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4)) - beta8*x7*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - mu7*x7); 
    dx8  = (  alpha8*(y8/eta8 - 1)*(1 - x8)  + beta9*x8*((v9*x9)^gamma9/(K9 + (v9*x9)^gamma9))  - mu8*x8);
    dx9  = (  alpha9*(y9/eta9 - 1)*(1 - x9)  - beta1*x9*((v1*x1)^gamma1/(K1 + (v1*x1)^gamma1))  + beta5*x9*((v5*x5)^gamma5/(K5 + (v5*x5)^gamma5)) - mu9*x9);
    dx10 = (  alpha10*(y10/eta10 - 1)*(1 - x10) + beta9*x10*((v9*x9)^gamma9/(K9 + (v9*x9)^gamma9)) + lambda1*(delta1/(delta1 + x10)) - mu10*x10); 
  
     return(list(c(dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx10)))  
  })
}
#*****************************************************************************************************

#*****************************************************
#STATES
######################
#Cancer
#high levels
# HER2   MDM2    AKT1     P21       CDK2     
x5 = 1; x6 = 1;  x8 = 1; x9 = 1;  x10 = 1;
#low levels
# TP53   ATR     ATM     BRCA1   CHEK1  
x1 = 0; x2 = 0; x3 = 0; x4 = 0;  x7 = 0; 
#Premalignant
#high levels
s <- c(x1=0,x2=0,x3=0,x4=0,x5=1,x6=1,x7=0,x8=1,x9=1,x10=1)  #cancer
#low levels
#Normal
#low levels
# HER2   MDM2    AKT1     P21       CDK2     
x5 = 0; x6 = 0;  x8 = 0; x9 = 0;  x10 = 0;
#high levels
# TP53   ATR     ATM     BRCA1   CHEK1  
x1 = 1; x2 = 1; x3 = 1; x4 = 1;  x7 = 1; 
#*******************************************************************

s <- c(x1=1,x2=1,x3=1,x4=1,x5=0,x6=0,x7=1,x8=0,x9=0,x10=0)  #normal

# Vector parameters
p <- c(alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beya7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,v1=v1,v2=v2,v3=v3,v4=v4,v5=v5,v6=v6,v7=v7,v8=v8,v9=v9,v10=v10,vy5=vy5,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gammay5=gammay5,K1=K1,K2=K2,K3=K3,K4=K4,K5=K5,K6=K6,K7=K7,K8=K8,K9=K9,K10=K10,Ky5=Ky5,lambda0=lambda0,lambda1=lambda1,delta0=delta0,delta1=delta1,y1=y1,y2=y2,y3=y3,y4=y4,y5=y5,y6=y6,y7=y7,y8=y8,y9=y9,y10=y10,eta1=eta1,eta2=eta2,eta3=eta3,eta4=eta4,eta5=eta5,eta6=eta6,eta7=eta7,eta8=eta8,eta9=eta9,eta10=eta10,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9,x10=x10)

#Vector states
s <- c(x1=0,x2=0,x3=0,x4=0,x5=0,x6=2,x7=0,x8=1,x9=1,x10=1)


#******************************************************
# Bifurcacion diagrams
#******************************************************

#*********************************
#System x1/x6
#*********************************

model2 <- function(t, state, parms){  
  with(as.list(c(state,parms)), {

    dx1  = (  alpha1*(y1/eta1 - 1)*(1 - x1) + beta2*x1*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  + beta4*x1*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4)) - beta5*x1*((v5*x5)^gamma5/(K5 + (v5*x5)^gamma5)) - beta6*x1*((v6*x6)^gamma6/(K6 + (v6*x6)^gamma6)) +  beta8*x1*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta9*x1*((v9*x9)^gamma9/(K9 + (v9*x9)^gamma9)) - mu1*x1);
    dx6  = (  alpha6*(y6/eta6 - 1)*(1 - x6) + beta1*x6*((v1*x1)^gamma1/(K1 + (v1*x1)^gamma1))  + beta3*x6*((v3*x3)^gamma3/(K3 + (v3*x3)^gamma3)) + beta8*x6*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta2*x6*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2)) - mu6*x6);

    x1 = max(min(x1,1),0);
    x6 = max(min(x6,1),0);
    
    return(list(c(dx1,dx6)))
  })
}

##################
x1 = 0; x6 = 0;
x2 = 0; x7 = 0;
x3 = 0; x8 = 0;
x4 = 0; x9 = 0;
x5 = 0; x10 = 0;

#precancer
#[1 0 1 0 0 0 1 0 0 1];
#normal
#[1 1 1 1 0 0 1 0 0 0];
y1=0.2; y6=y1; eta1=0.5; eta6=eta1; x2=1;x3=1;x4=0.6;x5=0;x7=1;x8=0.0;x9=0;x10=0;beta5=20.0;
s <- c(x1=1,x6=0.1)
p1 <- c(alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beya7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,v1=v1,v2=v2,v3=v3,v4=v4,v5=v5,v6=v6,v7=v7,v8=v8,v9=v9,v10=v10,vy5=vy5,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gammay5=gammay5,K1=K1,K2=K2,K3=K3,K4=K4,K5=K5,K6=K6,K7=K7,K8=K8,K9=K9,K10=K10,Ky5=Ky5,lambda0=lambda0,lambda1=lambda1,delta0=delta0,delta1=delta1,y1=y1,y2=y2,y3=y3,y4=y4,y5=y5,y6=y6,y7=y7,y8=y8,y9=y9,y10=y10,eta1=eta1,eta2=eta2,eta3=eta3,eta4=eta4,eta5=eta5,eta6=eta6,eta7=eta7,eta8=eta8,eta9=eta9,eta10=eta10,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9,x10=x10);
plane(parms=p1,odes=model2,state=s,xmax=1.0,ymax=2,ymin=0,x=1,y=2 )

mid <- newton(odes=model2,parms=p1,c(x1=0.5,x6=0.5),plot=T,)
hig <- newton(odes=model2,parms=p1,c(x1=0.8,x6=0.1),plot=T)
low <- newton(odes=model2,parms=p1,c(x1=0.0,x6=1.3),plot=T)


pmin = 0.0;
pmax = 30.0;
pymin = 0;
pymax = 1.5;

continue(state=hig, parms=p1, odes=model2, x="beta1", step=stepm, xmin=pmin, xmax=pmax,y="x1", ymin=pymin, ymax=pymax,main = "Genetic States x1 vs alpha6") # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p1, odes=model2, x="beta1", step=stepm, xmin=pmin, xmax=pmax,y="x1", ymin=pymin, ymax=pymax, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p1, odes=model2, x="beta1", step=stepm, xmin=pmin, xmax=pmax,y="x1", ymin=pymin, ymax=pymax, log="", time=0, positive=TRUE, add=TRUE)

continue(state=hig, parms=p1, odes=model2, x="beta1", step=stepm, xmin=pmin, xmax=pmax,y="x6", ymin=pymin, ymax=pymax,main = "Genetic States. x6 vs alpha6") # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p1, odes=model2, x="beta1", step=stepm, xmin=pmin, xmax=pmax,y="x6", ymin=pymin, ymax=pymax, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p1, odes=model2, x="beta1", step=stepm, xmin=pmin, xmax=pmax,y="x6", ymin=pymin, ymax=pymax, log="", time=0, positive=TRUE, add=TRUE)




#*********************************
#System x1/x9
#*********************************

model3 <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
   
    dx1  = (  alpha1*(y1/eta1 - 1)*(1 - x1)  + beta2*x1*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  + beta4*x1*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4)) - beta5*x1*((v5*x5)^gamma5/(K5 + (v5*x5)^gamma5)) - beta6*x1*((v6*x6)^gamma6/(K6 + (v6*x6)^gamma6)) +  beta8*x1*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta9*x1*((v9*x9)^gamma9/(K9 + (v9*x9)^gamma9)) - mu1*x1);
    dx9  = (  alpha9*(y9/eta9 - 1)*(1 - x9)  - beta1*x9*((v1*x1)^gamma1/(K1 + (v1*x1)^gamma1))  + beta5*x9*((v5*x5)^gamma5/(K5 + (v5*x5)^gamma5)) - mu9*x9);

    return(list(c(dx1,dx9)))
  })
}

##################
x1 = 0; x6 = 0;
x2 = 0; x7 = 0;
x3 = 0; x8 = 0;
x4 = 0; x9 = 0;
x5 = 0; x10 = 0;

#precancer
#[1 0 1 0 0 0 1 0 0 1];
#normal
#[1 1 1 1 0 0 1 0 0 0];
#cancer
#[0 0 0 0 0 1 0 1 1 1];

y1=0.2; y9=y1; eta1=0.5; eta9=eta1; x2=1;x3=1;x4=1;x5=0;x6=0;x7=1;x8=0;x10=0;beta5=20.0;
s <- c(x1=1,x9=0.1)
p1 <- c(alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beya7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,v1=v1,v2=v2,v3=v3,v4=v4,v5=v5,v6=v6,v7=v7,v8=v8,v9=v9,v10=v10,vy5=vy5,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gammay5=gammay5,K1=K1,K2=K2,K3=K3,K4=K4,K5=K5,K6=K6,K7=K7,K8=K8,K9=K9,K10=K10,Ky5=Ky5,lambda0=lambda0,lambda1=lambda1,delta0=delta0,delta1=delta1,y1=y1,y2=y2,y3=y3,y4=y4,y5=y5,y6=y6,y7=y7,y8=y8,y9=y9,y10=y10,eta1=eta1,eta2=eta2,eta3=eta3,eta4=eta4,eta5=eta5,eta6=eta6,eta7=eta7,eta8=eta8,eta9=eta9,eta10=eta10,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9,x10=x10);
plane(parms=p1,odes=model3,state=s,xmax=1.0,ymax=3.0,ymin=0,x=1,y=2 )


mid <- newton(odes=model3,parms=p1,c(x1=0.2,x9=.5),plot=T)
low <- newton(odes=model3,parms=p1,c(x1=1.0,x9=0.3),plot=T)
hig <- newton(odes=model3,parms=p1,c(x1=0.2,x9=1.5),plot=T)
low<- c(x1=1.0,x9=0.0)

pmin = 0;
pmax = 2.0;
pymin = -0.2;
pymax = 2.0;

continue(state=low,tstep = 1e-5, parms=p1, odes=model3, x="alpha9", step=stepm, xmin=pmin, xmax=pmax,y="x1", ymin=pymin, ymax=pymax,main = "") # log="", time=0, positive=TRUE, add=TRUE)
continue(state=hig,tstep = 1e-6, parms=p1, odes=model3, x="alpha9", step=stepm, xmin=pmin, xmax=pmax,y="x1", ymin=pymin, ymax=pymax, log="", time=0, positive=TRUE, add=TRUE)


continue(state=low,tstep = 1e-5, parms=p1, odes=model3, x="alpha9", step=stepm, xmin=pmin, xmax=pmax,y="x9", ymin=pymin, ymax=pymax,main = "") # log="", time=0, positive=TRUE, add=TRUE)
continue(state=hig,tstep = 1e-6, parms=p1, odes=model3, x="alpha9", step=stepm, xmin=pmin, xmax=pmax,y="x9", ymin=pymin, ymax=pymax, log="", time=0, positive=TRUE, add=TRUE)


#*********************************
#System x2/x4
#*********************************

model4 <- function(t, state, parms){  
  with(as.list(c(state,parms)), {

    dx2  = (  alpha2*(y2/eta2 - 1)*(1 - x2)  + beta4*x4*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4))  - mu2*x2);
    dx4  = (  alpha4*(y4/eta4 - 1)*(1 - x4)  + beta2*x4*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  - beta8*x4*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta10*x4*((v10*x10)^gamma10/(K10 + (v10*x10)^gamma10)) - mu4*x4);
    
    return(list(c(dx2,dx4)))
  })
}

##################
x1 = 0; x6 = 0;
x2 = 0; x7 = 0;
x3 = 0; x8 = 0;
x4 = 0; x9 = 0;
x5 = 0; x10 = 0;

#precancer
#[1 0 1 0 0 0 1 0 0 1];
#normal
#[1 1 1 1 0 0 1 0 0 0];
#cancer
#[0 0 0 0 0 1 0 1 1 1];

y2=0.2; y4=y2; eta2=0.25;  x1=0;x3=0;x5=0;x6=1;x7=0;x8=1;x9=1;x10=1;beta5=20.0;
s <- c(x2=1,x4=0.1)
p1 <- c(alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4,alpha5=alpha5,alpha6=alpha6,alpha7=alpha7,alpha8=alpha8,alpha9=alpha9,alpha10=alpha10,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,beta5=beta5,beta6=beta6,beya7=beta7,beta8=beta8,beta9=beta9,beta10=beta10,v1=v1,v2=v2,v3=v3,v4=v4,v5=v5,v6=v6,v7=v7,v8=v8,v9=v9,v10=v10,vy5=vy5,gamma1=gamma1,gamma2=gamma2,gamma3=gamma3,gamma4=gamma4,gamma5=gamma5,gamma6=gamma6,gamma7=gamma7,gamma8=gamma8,gamma9=gamma9,gamma10=gamma10,gammay5=gammay5,K1=K1,K2=K2,K3=K3,K4=K4,K5=K5,K6=K6,K7=K7,K8=K8,K9=K9,K10=K10,Ky5=Ky5,lambda0=lambda0,lambda1=lambda1,delta0=delta0,delta1=delta1,y1=y1,y2=y2,y3=y3,y4=y4,y5=y5,y6=y6,y7=y7,y8=y8,y9=y9,y10=y10,eta1=eta1,eta2=eta2,eta3=eta3,eta4=eta4,eta5=eta5,eta6=eta6,eta7=eta7,eta8=eta8,eta9=eta9,eta10=eta10,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9,x10=x10);
plane(parms=p1,odes=model4,state=s,xmax=1.0,ymax=2,ymin=0,x=1,y=2 )

mid <- newton(odes=model4,parms=p1,c(x2=0.8,x4=.5),plot=T)
low <- newton(odes=model4,parms=p1,c(x2=0.0,x4=0.3),plot=T)
hig <- newton(odes=model4,parms=p1,c(x2=0.5,x4=1.5),plot=T)

pmin = -2.0;
pmax =1.0;
pymin = 0;
pymax = 2;

continue(state=hig, parms=p1, odes=model4, x="alpha4", step=stepm, xmin=pmin, xmax=pmax,y="x2", ymin=pymin, ymax=pymax,main = "Genetic States. X1->X6") # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p1, odes=model4, x="alpha4", step=stepm, xmin=pmin, xmax=pmax,y="x2", ymin=pymin, ymax=pymax, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p1, odes=model4, x="alpha4", step=stepm, xmin=pmin, xmax=pmax,y="x2", ymin=pymin, ymax=pymax, log="", time=0, positive=TRUE, add=TRUE)

