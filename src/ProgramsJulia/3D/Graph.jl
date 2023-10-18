function grafs(cell,parms,Msure,tim)

TcellG = zeros(Nx,Ny,Nz);
TcellS = zeros(Nx,Ny,Nz);
TcellT = zeros(Nx,Ny,Nz);
TcellE = zeros(Nx,Ny,Nz);
TcellU = zeros(Nx,Ny,Nz);
TcellM = zeros(Nx,Ny,Nz);

    for ii =2:Nx-1
        for jj = 2:Ny-1 
            for kk = 2:Ny-1                   
                TcellG[ii,jj,kk] = sum(cell[ii,jj,kk].genemutation[:]);
                TcellS[ii,jj,kk] = sum(cell[ii,jj,kk].segregation[:]);
                TcellT[ii,jj,kk] = cell[ii,jj,kk].geneexpression[11];
                TcellE[ii,jj,kk] = cell[ii,jj,kk].E; 
                TcellU[ii,jj,kk] = sum(cell[ii,jj,kk].thresholdexpression[:])/10;
                TcellM[ii,jj,kk] = cell[ii,jj,kk].Type[2];
            end
        end
    end

    plotly()
    cell, ik,jk,kk = neibs(cell,parms);
    lix = length(ik);  
    #scar = zeros(lix,1);
    #scarx = scary =scarz = cmm =zeros(lix,1) 
    alphax = zeros(lix,1);
    cm = cgrad(:viridis)
    colors = cm;        
#           
    fig=figure()
#        
    for k=minimum(kk):maximum(kk)
#                  
    #    cc = max.(Int64(round(256*(TcellG[ik[k],jk[k],kk[k]] +1)/(maximum(TcellG) +1) )),1);
      
    #    alphax[k] = (TcellG[ik[k],jk[k],kk[k]]+1)/(maximum(TcellG)+1); 
       
    #    scatter3D(ik[k],jk[k],kk[k]*TcellM[ik[k],jk[k],kk[k]],color=(red(colors[cc]), green(colors[cc]), blue(colors[cc])) ) #,alpha = alphax[k] )

    matshow(TcellG[:,:,k])
    colorbar()
    PyPlot.savefig("A05-"*string(parms.var9)*"-kappa-"*string(parms.var17)*"-chi-"*string(parms.var14)*"-xi-"*string(parms.var15)*"-t-"*string(tim)*"-K"*string(k)*".pdf")
    
    end
   

############
 
for k=minimum(kk):maximum(kk)

    #cc = max.(Int64(round( 256*(TcellT[ik[k],jk[k],kk[k]]+ 1 )/(maximum(TcellT)+1) )),1);
    #alphax[k] = (TcellT[ik[k],jk[k],kk[k]] +1)/(maximum(TcellT)+1); 
    #scatter3D(ik[k],jk[k],kk[k]*TcellM[ik[k],jk[k],kk[k]],color=(red(colors[cc]), green(colors[cc]), blue(colors[cc])),alpha = alphax[k] )

    matshow(TcellT[:,:,k])
    colorbar()
    PyPlot.savefig("B05-"*string(parms.var9)*"-kappa-"*string(parms.var17)*"-chi-"*string(parms.var14)*"-xi-"*string(parms.var15)*"-t-"*string(tim)*"-K"*string(k)*".pdf")
         
end
 



###################################################
#### Estimation of Shannon Index ##################
#####   and Fractal dimension #####################

inda = reshape(TcellG,parms.var10*parms.var11*parms.var11);
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
SI = reshape((TcellS/parms.var12),parms.var10*parms.var11*parms.var11);
EM = reshape((TcellE/parms.var12),parms.var10*parms.var11*parms.var11);

meansi = sum(TcellS)/indt;     
meanE = mean(EM);
varsi = var(SI[findall(inda.>0)]);

#########################################
tcellcan = sum(sign.(TcellG));
indnormal = reshape(TcellT,parms.var10*parms.var11*parms.var11);
indtnormal = length(findall(indnormal.==2.0))/tcellcan;
indtprecancer = length(findall(indnormal.==8))/tcellcan;
indtcancer = length(findall(indnormal.==10))/tcellcan;
#######################################
typomean = filter!(x->x>0,reshape(TcellT,parms.var10*parms.var11*parms.var11));
mutmean = filter!(x->x>0,reshape(TcellG,parms.var10*parms.var11*parms.var11));
thersmean = filter!(x->x>0,reshape(TcellU,parms.var10*parms.var11*parms.var11));

#

Msure = vcat(Msure, [tim shindex bxt meansi varsi meanE indtnormal indtprecancer indtcancer mean(typomean) mean(mutmean) mean(thersmean)]);
#Msure[coum,:] = [tim shindex bxt meansi varsi meanE indtnormal indtprecancer indtcancer mean(typomean) mean(mutmean) mean(thersmean)];

#  time, shannon idex, fractal dimension, mean of SI, var of SI, 
#...Feno normal, Feno Precancer,  Feno cancer 

return Msure

    
end
