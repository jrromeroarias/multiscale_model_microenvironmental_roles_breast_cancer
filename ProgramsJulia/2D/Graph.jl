function grafs(cell,parms,Msure,tim)

TcellG = zeros(Nx,Ny);
TcellS = zeros(Nx,Ny);
TcellT = zeros(Nx,Ny);
TcellE = zeros(Nx,Ny);
TcellU = zeros(Nx,Ny);
    for ii =2:Nx-1
        for jj = 2:Ny-1                    
            TcellG[ii,jj] = sum(cell[ii,jj].genemutation[:]);
            TcellS[ii,jj] = sum(cell[ii,jj].segregation[:]);
            TcellT[ii,jj] = cell[ii,jj].geneexpression[11];
            TcellE[ii,jj] = cell[ii,jj].E; 
            TcellU[ii,jj] = sum(cell[ii,jj].thresholdexpression[:])/10;

        end
    end

matshow(TcellG)
colorbar()
PyPlot.savefig("Muta-beta-"*string(parms.var9)*"-kappa-"*string(parms.var17)*"-chi-"*string(parms.var14)*"-xi-"*string(parms.var15)*"-alphaN-"*string(parms.var1)*"-lambdaN-"*string(parms.var2)*"-alphaE-"*string(parms.var5)*"-lambdaE-"*string(parms.var6)*"-t-"*string(tim)*".pdf")    
matshow(TcellS)
colorbar()
PyPlot.savefig("Segre-beta-"*string(parms.var9)*"-kappa-"*string(parms.var17)*"-chi-"*string(parms.var14)*"-xi-"*string(parms.var15)*"-alphaN-"*string(parms.var1)*"-lambdaN-"*string(parms.var2)*"-alphaE-"*string(parms.var5)*"-lambdaE-"*string(parms.var6)*"-t-"*string(tim)*".pdf") 
matshow(TcellT)
colorbar()
PyPlot.savefig("Feno-beta-"*string(parms.var9)*"-kappa-"*string(parms.var17)*"-chi-"*string(parms.var14)*"-xi-"*string(parms.var15)*"-alphaN-"*string(parms.var1)*"-lambdaN-"*string(parms.var2)*"-alphaE-"*string(parms.var5)*"-lambdaE-"*string(parms.var6)*"-t-"*string(tim)*".pdf") 


###################################################
#### Estimation of Shannon Index ##################
#####   and Fractal dimension #####################

inda = reshape(TcellG,parms.var10*parms.var11);
indt = length(findall(inda.>0));
pdi = zeros(Int64(maximum(TcellG)));

for ip = 1:length(pdi)
    pdi[ip] = length(findall(inda.==ip))/indt;
end
indpdi = filter!(x->x>0,pdi);
shindex = -1.0*sum(indpdi.*log.(indpdi));
bxt = mean(filter!(x->x>0,boxcount(TcellT)));       
#########################################
########################################
SI = reshape((TcellS/parms.var12),parms.var10*parms.var11);
EM = reshape((TcellE/parms.var12),parms.var10*parms.var11);

meansi = sum(TcellS)/indt;     
meanE = mean(EM);
varsi = var(SI[findall(inda.>0)]);

#########################################
tcellcan = sum(sign.(TcellG));
indnormal = reshape(TcellT,parms.var10*parms.var11);
indtnormal = length(findall(indnormal.==2.0))/tcellcan;
indtprecancer = length(findall(indnormal.==8))/tcellcan;
indtcancer = length(findall(indnormal.==10))/tcellcan;
#######################################
typomean = filter!(x->x>0,reshape(TcellT,parms.var10*parms.var11));
mutmean = filter!(x->x>0,reshape(TcellG,parms.var10*parms.var11));
thersmean = filter!(x->x>0,reshape(TcellU,parms.var10*parms.var11));

#

Msure = vcat(Msure, [tim shindex bxt meansi varsi meanE indtnormal indtprecancer indtcancer mean(typomean) mean(mutmean) mean(thersmean)]);
#  time, shannon idex, fractal dimension, mean of SI, var of SI, 
#...Feno normal, Feno Precancer,  Feno cancer 

return Msure

    
end