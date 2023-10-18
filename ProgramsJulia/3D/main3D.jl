#!/usr/local/bin/julia
####
# Microenvironmental factors in cell segregation and heterogeneity in breast cancer development
## Propagation of Breast Cancer
### by Roberto Romero - IIMAS/UNAM
#####
using Plots,PyPlot
using Compat, Distributions
using ProgressMeter, Random
using CSV, DataFrames 
####

include("Initial.jl"); #Initiation of variables
include("Neib.jl");    #Check the number of neigbors with cancer and update the celular clock
include("Repro.jl");   #Calcules proliferation. Solve Eqs.(3)-(4) with Eq. (12)
include("Dead.jl");    #Proliferation function. Solve Eq. (2)
include("DiffG.jl");   #Solve global reaction-difusion equations. Solve Eq. (1)
include("DiffL.jl");   #Solve local reaction-difusion equations.
include("Express.jl"); #Calculated gene expression. Solve Eq.(13)
include("Muta.jl");    #Calculated gene mutations. Solve Eqs. (7)-(11)
include("Segre.jl");   #Calculated the segregation index, solve Eq.(5) and Eq. (17)
include("Graph.jl");   #Export graphics
include("Boxfrac.jl"); #Calculated fractal dimension
include("Fravol.jl");  #Estimated volume fraction of tumor 


mutable struct ATTR
    N::Float64
    M::Float64
    E::Float64
    R::Float64
    Type::Vector{Int64}
    segregation::Vector{Float64}
    genemutation::Vector{Int64}
    geneexpression::Vector{Float64}
    thresholdexpression::Vector{Float64}
end

mutable struct parameters
    var1::Float64 
    var2::Float64
    var3::Float64
    var4::Float64
    var5::Float64
    var6::Float64
    var7::Float64
    var8::Float64
    var9::Float64
    var10::Int64
    var11::Int64
    var12::Int64
    var13::Float64
    var14::Float64
    var15::Float64
    var16::Int64
    var17::Float64
end
##### Parameters

Nx = 220;   # lattice size
Ny = 220;
Nz = Ny;
tmax = 1500;
genes = 10;  #number of genes
alphaN = 8.0/500; ##nutrient consumption rate
alphaM = 4.0/500; ##oxygen consumption rate
alphaE = 16.0/500; ##estrogens consumption rate
lambdaN = 100;  ##excess consumption of N by cancer cells
lambdaM = 50;  #excess consumption of M by cancer cells
lambdaE = 50;  #excess consumption of E by cancer cells
thetadiv = 0.3;  #shape parameter; division response to nutrients
thetadie = 0.01; #shape parameter; death response to oxygen
beta = 1.0;  # segregation importance
eta = 0.5;   # Threshold for gene expression
chi = 4.0; # Threshold susceptibility for mutation and microenviroment
xi = 2.0;     # Threshold susceptibility heritage
radio = 20 # scale of neighbourhood
kap = 0.1; # estrogen charge capacity 

parms = parameters(alphaN,lambdaN,alphaM,lambdaM,alphaE,lambdaE,thetadiv,thetadie,beta,Nx,Ny,genes,eta,chi,xi,radio,kap);

open("C-parmsp.txt","w") do file
    write(file,"........alphaN.lambdaN.alphaM.lambdaM.alphaE.lambdaE.tdiv.tdie.beta...Nx...Ny...genes.eta..chi..xi..radio..kap \n")
    write(file,string(parms));
end

cell0 = InitialV(parms);

#############################
#############################
#for ixi in (0.5,1.0,2.0,4.0,6.0,8.0,10.0,12.0)
    
    for ichi in (2.0)
    
#    xi = ixi;    
    chi = ichi;
    
parms = parameters(alphaN,lambdaN,alphaM,lambdaM,alphaE,lambdaE,thetadiv,thetadie,beta,Nx,Ny,genes,eta,chi,xi,radio,kap);

cellT = deepcopy(cell0);

MsureT = zeros(1,12);


for t = 1:tmax    
 
    cellT,ik,jk,kk = neibs(cellT,parms);
    for mm = 1: length(ik)
        
        act = rand(1:2);
        
        if (act == 1)
            if (cellT[ik[mm],jk[mm],kk[mm]].R >35.0 )
                cellT = reproduction(cellT,parms,ik[mm],jk[mm],kk[mm]);
            end
        elseif (act == 10 && length(ik)>50 && rand()<=exp( -( cellT[ik[mm],jk[mm],kk[mm]].M / parms.var8 )^2 )) #death
                cellT = death(cellT,parms,ik[mm],jk[mm],kk[mm]);
        end
                cellT = Difflocal(cellT,parms,ik[mm],jk[mm],kk[mm],4);
    end

        cellT = Diffglobal(cellT,parms,1); 

    if mod(t,40)==0
         #count = Int64(count + 1);
         MsureT = grafs(cellT,parms,MsureT,t);
		 #CSV.write("A-beta-"*string(parms.var9)*"-kappa-"*string(parms.var17)*"-chi-"*string(parms.var14)*"-xi-"*string(parms.var15)*".csv", DataFrame(MsureT),header = false);
    	 
    end
#println(t)
end


    end
#end


    

     

 
