function  mutations(cell,parms,i,j,im,jm)
#This function returns the new probability of mutation and increse the gene mutations    
#Genes that will mutate

p = exp( -( cell[i,j].N / parms.var7 )^2 );
p = max(0.001,p); #savers
p = min(0.999,p); 

cell, numa = segregationlocal(cell,parms,i,j);
    

for il=1:parms.var12   # all genes can change at the same time
 
    pp = ( p + cell[i,j].segregation[il] - p* cell[i,j].segregation[il] );
    gene_mut = min.(rand( NegativeBinomial(numa[1],1-p),parms.var12),parms.var12);
    
     num_p = 0;
     num_g = 0;
     
    if (  isempty(gene_mut[il])==0 && gene_mut[il]>=1 && pp>rand() ) #p>probm # nut<θdiv*sqrt(-log(probm))
        # Tau-Leaping algorithm
        tau = -log(rand());
        num_p = rand(Poisson(tau));
        num_g = rand(Geometric(p));    
    end
        
    cell[i+im,j+jm].genemutation[il] = cell[i,j].genemutation[il] + num_p*num_g;
    cell[i+im,j+jm].thresholdexpression[il] = cell[i,j].thresholdexpression[il];               

end        

    cell = expression(cell,parms,i,j,im,jm);

   
return cell        
        
        
end    