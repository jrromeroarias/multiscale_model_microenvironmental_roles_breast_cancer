function InitialV(parms)
   
## Initial values for microenvironment
Naux = zeros(parms.var10,parms.var11,parms.var11);
Maux = zeros(parms.var10,parms.var11,parms.var11);
Eaux = zeros(parms.var10,parms.var11,parms.var11);
A = exp(parms.var1*parms.var10)/(exp(parms.var1*parms.var10) + exp(-parms.var1*parms.var10));
B = 1 - (exp(parms.var1*parms.var10)/(exp(parms.var1*parms.var10) + exp(-parms.var1*parms.var10)));
x = range(1, length=parms.var10, stop=parms.var10);
Con = A*exp.(-parms.var1*x) + B*exp.(parms.var1*x);
Con2 = ones(parms.var10,parms.var11,parms.var11);
#for is = 1:parms.var10
#    for js = 1:parms.var11
        for ks = 1:parms.var11
            Naux[:,:,ks] = Con[ks].*Con2[:,:,ks] .+ 0.01*rand(-25:25);  # nutrients
            Maux[:,:,ks] = Con[ks].*Con2[:,:,ks] .+ 0.01*rand(-25:25);  # oxigen
            Eaux[:,:,ks] = Con[ks].*Con2[:,:,ks] .+ 0.01*rand(-25:25);  # Estrogens
        end
#    end
#end

Naux[:,:,:] = max.(Naux[:,:,:],0);
Maux[:,:,:] = max.(Maux[:,:,:],0);
Eaux[:,:,:] = max.(Eaux[:,:,:],0);
Naux[:,:,:] = min.(Naux[:,:,:],1);
Maux[:,:,:] = min.(Maux[:,:,:],1);
Eaux[:,:,:] = min.(Eaux[:,:,:],1);
    
## Inititation
AB1 = zeros(parms.var12,1);
AB2 = zeros(parms.var12+1,1);
cellInt = Array{ATTR}(undef, parms.var10*parms.var10*parms.var10);
for i = 1:parms.var10^3
   cellInt[i] = ATTR(0,0,0,0,[0,0,0],AB1[:],AB1[:],AB2[:],AB1[:]); # for type [normal,cancer,dead]
end
cellInt = reshape(cellInt,parms.var10,parms.var10,parms.var10);    

    
for i = 1:parms.var10
    for j = 1:parms.var11
        for k =1:parms.var11
            cellInt[i,j,k].N = Naux[i,j,k];  # Glucose
            cellInt[i,j,k].M = Maux[i,j,k];  # Oxygen
            cellInt[i,j,k].E = Eaux[i,j,k];  # Estrogens
            cellInt[i,j,k].R = 20*rand();  # Clock
        end
    end
end

# Initial cancer cell
pos_i = Int64(rand((parms.var10/2)-20:(parms.var10/2)+20));
pos_j = Int64(rand((parms.var11/2)-20:(parms.var11/2)+20));
pos_k = Int64(rand((parms.var11/2)-20:(parms.var11/2)+20));
cellInt[pos_i,pos_j,pos_k].Type[1] = 0;
cellInt[pos_i,pos_j,pos_k].Type[2] = 1;
cellInt[pos_i,pos_j,pos_k].genemutation[:] = [1 1 1 1 1 1 1 1 1 1];
cellInt[pos_i,pos_j,pos_k].geneexpression[:] = [0 0 0 0 1 1 0 1 1 1 10];
cellInt[pos_i,pos_j,pos_k].R = 36;

#cellInt = reshape(cellInt,parms.var10*parms.var10*parms.var10);

#CSV.write("ABcell.csv", DataFrame(cellInt),header = false);
        
    return cellInt
    
end