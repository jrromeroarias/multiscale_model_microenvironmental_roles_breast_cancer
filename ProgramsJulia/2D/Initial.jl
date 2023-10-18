function InitialV(parms)
   
## Initial values for microenvironment
Naux = zeros(parms.var10,parms.var11);
Maux = zeros(parms.var10,parms.var11);
Eaux = zeros(parms.var10,parms.var11);
A = exp(parms.var1*parms.var10)/(exp(parms.var1*parms.var10) + exp(-parms.var1*parms.var10));
B = 1 - (exp(parms.var1*parms.var10)/(exp(parms.var1*parms.var10) + exp(-parms.var1*parms.var10)));
x = range(1, length=parms.var10, stop=parms.var10);
Con = A*exp.(-parms.var1*x) + B*exp.(parms.var1*x);
Naux = Con.*ones(parms.var10,parms.var11) + 0.01*rand(-25:25,parms.var10,parms.var11);  # nutrients
Maux = Con.*ones(parms.var10,parms.var11) + 0.01*rand(-25:25,parms.var10,parms.var11);  # oxigen
Eaux = Con.*ones(parms.var10,parms.var11) + 0.01*rand(-25:25,parms.var10,parms.var11);  # Estrogens
Naux[:,:] = max.(Naux[:,:],0);
Maux[:,:] = max.(Maux[:,:],0);
Eaux[:,:] = max.(Eaux[:,:],0);
Naux[:,:] = min.(Naux[:,:],1);
Maux[:,:] = min.(Maux[:,:],1);
Eaux[:,:] = min.(Eaux[:,:],1);
    
## Inititation
AB1 = zeros(parms.var12,1);
AB2 = zeros(parms.var12+1,1);
cellInt = Array{ATTR}(undef, parms.var10*parms.var11);
for i = 1:parms.var10^2
   cellInt[i] = ATTR(0,0,0,0,[0,0,0],AB1[:],AB1[:],AB2[:],AB1[:]); # for type [normal,cancer,dead]
end
cellInt = reshape(cellInt,parms.var10,parms.var10);    

    
for i = 1:parms.var10
    for j = 1:parms.var11
        cellInt[i,j].N = Naux[i,j];  # Glucose
        cellInt[i,j].M = Maux[i,j];  # Oxygen
        cellInt[i,j].E = Eaux[i,j];  # Estrogens
        cellInt[i,j].R = 20*rand();  # Clock
    end
end

# Initial cancer cell
pos_i = Int64(rand((parms.var10/2)-20:(parms.var10/2)+20));
pos_j = Int64(rand((parms.var11/2)-20:(parms.var11/2)+20));
cellInt[pos_i,pos_j].Type[1] = 0;
cellInt[pos_i,pos_j].Type[2] = 1;
cellInt[pos_i,pos_j].genemutation[:] = [1 1 1 1 1 1 1 1 1 1];
cellInt[pos_i,pos_j].geneexpression[:] = [0 0 0 0 1 1 0 1 1 1 10];
cellInt[pos_i,pos_j].R = 36;
    
    
    return cellInt
    
end