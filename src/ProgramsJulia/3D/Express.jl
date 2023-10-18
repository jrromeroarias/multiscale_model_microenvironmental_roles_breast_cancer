function expression(cell,parms,i,j,k,im,jm,km) #Y0,η0,G0)
 
#    %******** List of parameter
alpha1 = -1;  beta1 = 9.0;  mu1 = 0.2;  v1 = 1.5;  gamma1 = 6;  K1 = 0.5;
alpha2 = -1;  beta2 = 0.2;   mu2 = 0.2;  v2 = 1;    gamma2 = 1;  K2 = 0.5;
alpha3 = -1;  beta3 = 0.2;   mu3 = 0.2;  v3 = 1;    gamma3 = 1;  K3 = 0.5;
alpha4 = -1;  beta4 = 0.2;   mu4 = 0.2;  v4 = 1;    gamma4 = 3;  K4 = 0.5;
alpha5 = 1;  beta5 = 20.0;  mu5 = 0.2;  v5 = 1;    gamma5 = 1;  K5 = 0.5;
alpha6 = 1;  beta6 = 4.0;   mu6 = 0.2;  v6 = 1.2;  gamma6 = 4;  K6 = 0.5;
alpha7 = -1;  beta7 = 0.0;   mu7 = 0.2;  v7 = 1;    gamma7 = 3;  K7 = 0.5;
alpha8 = 1;  beta8 = 0.2;   mu8 = 0.2;  v8 = 1;    gamma8 = 1;  K8 = 0.5;
alpha9 = 1;  beta9 = 3.0;   mu9 = 0.2;  v9 = 0.5;  gamma9 = 2;  K9 = 0.5;
alpha10 = 1; beta10 = 4.0;  mu10 = 0.2; v10 = 1;   gamma10 = 2; K10 = 0.5;
    
vy5 = 0.5; gammay5 = 1;  Ky5 = 0.5;
vy10 = 1;  gammay10 = 1; Ky10 = 5;

lambda0 = 0.4; delta0 = 1;
lambda1 = 0.4; delta1 = 0.1;
    
#******** 
y1  = cell[i,j,k].M;   # TP53,   stress, oxygen                   # oxygen
y2  = cell[i,j,k].E/(cell[i,j,k].E + parms.var17);   # ATR,                                     # estrogens
y3  = cell[i,j,k].E/(cell[i,j,k].E + parms.var17);   # ATM,                                     # estrogens
y4  = cell[i,j,k].E/(cell[i,j,k].E + parms.var17);   # BRCA1, Methylation --->                  # estrogens
y5  = cell[i,j,k].E/(cell[i,j,k].E + parms.var17);   # HER2, estrogens, growth factor, glucose  # estrogens
y6  = cell[i,j,k].M;   # MDM2, phosphorylation rates---->         # oxygen
y7  = cell[i,j,k].E/(cell[i,j,k].E + parms.var17);   # CHEK1                                    # estrogens
y8  = cell[i,j,k].M;   # AKT1                                     # oxygen
y9  = cell[i,j,k].M;   # P21                                      # oxygen
y10 = cell[i,j,k].E/(cell[i,j,k].E + parms.var17);   # CDK2                                     # estrogens
#***********

    DT = 3600;
    dt = 0.001;
    

x1 = rand();#cell[i,j].geneexpression[1]; 
x2 = rand();#cell[i,j].geneexpression[2]; 
x3 = rand();#cell[i,j].geneexpression[3];
x4 = rand();#cell[i,j].geneexpression[4];
x5 = rand();#cell[i,j].geneexpression[5]; 
x6 = rand();#cell[i,j].geneexpression[6];
x7 = rand();#cell[i,j].geneexpression[7];
x8 = rand();#cell[i,j].geneexpression[8];
x9 = rand();#cell[i,j].geneexpression[9];
x10 = rand();#cell[i,j].geneexpression[10];
ftp = cell[i,j,k].geneexpression[11];
    
#eta1 = max.(parms.var13 - cell[i,j].thresholdexpression[1], 0.05); 
#eta2 = max.(parms.var13 - cell[i,j].thresholdexpression[2], 0.05); 
#eta3 = max.(parms.var13 - cell[i,j].thresholdexpression[3], 0.05); 
#eta4 = max.(parms.var13 - cell[i,j].thresholdexpression[4], 0.05); 
#eta5 = max.(parms.var13 - cell[i,j].thresholdexpression[5], 0.05); 
#eta6 = max.(parms.var13 - cell[i,j].thresholdexpression[6], 0.05); 
#eta7 = max.(parms.var13 - cell[i,j].thresholdexpression[7], 0.05); 
#eta8 = max.(parms.var13 - cell[i,j].thresholdexpression[8], 0.05); 
#eta9 = max.(parms.var13 - cell[i,j].thresholdexpression[9], 0.05); 
#eta10 = max.(parms.var13 - cell[i,j].thresholdexpression[10], 0.05);     

eta1 = cell[i,j,k].thresholdexpression[1]; 
eta2 = cell[i,j,k].thresholdexpression[2]; 
eta3 = cell[i,j,k].thresholdexpression[3]; 
eta4 = cell[i,j,k].thresholdexpression[4]; 
eta5 = cell[i,j,k].thresholdexpression[5]; 
eta6 = cell[i,j,k].thresholdexpression[6]; 
eta7 = cell[i,j,k].thresholdexpression[7]; 
eta8 = cell[i,j,k].thresholdexpression[8]; 
eta9 = cell[i,j,k].thresholdexpression[9]; 
eta10 = cell[i,j,k].thresholdexpression[10];




   #***********    
     
    for i=1:DT-1
    %#0  = dt*i;    
    x1  = x1  + dt*(  alpha1*(y1/eta1 - 1)*(1 - x1)  + beta2*x1*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  + beta4*x1*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4)) - beta5*x1*((v5*x5)^gamma5/(K5 + (v5*x5)^gamma5)) - beta6*x1*((v6*x6)^gamma6/(K6 + (v6*x6)^gamma6)) +  beta8*x1*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta9*x1*((v9*x9)^gamma9/(K9 + (v9*x9)^gamma9)) - mu1*x1);
    x2  = x2  + dt*(  alpha2*(y2/eta2 - 1)*(1 - x2)  + beta4*x4*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4))  - mu2*x2);
    x3  = x3  + dt*(  alpha3*(y3/eta3 - 1)*(1 - x3)  + beta1*x3*((v1*x1)^gamma1/(K1 + (v1*x1)^gamma1))  - mu3*x3);
    x4  = x4  + dt*(  alpha4*(y4/eta4 - 1)*(1 - x4)  + beta2*x4*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  - beta8*x4*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta10*x4*((v10*x10)^gamma10/(K10 + (v10*x10)^gamma10)) - mu4*x4);
    x5  = x5  + dt*(  alpha5*(y5/eta5 - 1)*(1 - x5)  + lambda0*(x5/(delta0 + x5))*((vy5*y5)^gammay5/(Ky5     + (vy5*y5)^gammay5)) - mu5*x5);
    x6  = x6  + dt*(  alpha6*(y6/eta6 - 1)*(1 - x6)  + beta1*x6*((v1*x1)^gamma1/(K1 + (v1*x1)^gamma1))  + beta3*x6*((v3*x3)^gamma3/(K3 + (v3*x3)^gamma3)) + beta8*x6*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - beta2*x6*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2)) - mu6*x6);
    x7  = x7  + dt*(  alpha7*(y7/eta7 - 1)*(1 - x7)  + beta2*x7*((v2*x2)^gamma2/(K2 + (v2*x2)^gamma2))  + beta4*x7*((v4*x4)^gamma4/(K4 + (v4*x4)^gamma4)) - beta8*x7*((v8*x8)^gamma8/(K8 + (v8*x8)^gamma8)) - mu7*x7); 
    x8  = x8  + dt*(  alpha8*(y8/eta8 - 1)*(1 - x8)  + beta9*x8*((v9*x9)^gamma9/(K9 + (v9*x9)^gamma9))  - mu8*x8);
    x9  = x9  + dt*(  alpha9*(y9/eta9 - 1)*(1 - x9)  - beta1*x9*((v1*x1)^gamma1/(K1 + (v1*x1)^gamma1))  + beta5*x9*((v5*x5)^gamma5/(K5 + (v5*x5)^gamma5)) - mu9*x9);
    x10 = x10 + dt*( alpha10*(y10/eta10 - 1)*(1 - x10) + beta9*x10*((v9*x9)^gamma9/(K9 + (v9*x9)^gamma9)) + lambda1*(delta1/(delta1 + x10)) - mu10*x10); 
    
    x1  = max.(x1,0);  x1  = min.(x1,1.0);
    x2  = max.(x2,0);  x2  = min.(x2,1.0);
    x3  = max.(x3,0);  x3  = min.(x3,1.0);
    x4  = max.(x4,0);  x4  = min.(x4,1.0);
    x5  = max.(x5,0);  x5  = min.(x5,1.0);    
    x6  = max.(x6,0);  x6  = min.(x6,1.0);
    x7  = max.(x7,0);  x7  = min.(x7,1.0);
    x8  = max.(x8,0);  x8  = min.(x8,1.0);
    x9  = max.(x9,0);  x9  = min.(x9,1.0);
    x10 = max.(x10,0); x10 = min.(x10,1.0);
       
    end
  
    ## This is for fenotype
    umb = 0.8; # umbral 
    #ftp = 8;  #precancer
    if (x1>=umb && x2>=umb && x3>=umb && x4>=umb && x5<umb && x6<umb && x7>=umb && x8<umb && x9<umb && x10<umb)
        ftp = 2; #normal
    elseif  (x1<umb && x2<umb && x3<umb && x4<umb && x5>=umb && x6>=umb && x7<umb && x8>=umb && x9>=umb && x10>=umb)
        ftp = 10; # Cancer
    elseif (x1<umb && x2<umb && x3<umb && x4<umb && x5<umb && x6<umb && x7<umb && x8<umb && x9<umb && x10<umb)
        ftp = 8;  #precancer
    end
#Ecell = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 ftp];

    cell[i+im,j+jm,k+km].geneexpression = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,ftp];
 
    return cell  

end

