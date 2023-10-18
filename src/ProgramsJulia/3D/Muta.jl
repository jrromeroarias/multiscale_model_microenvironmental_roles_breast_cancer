function  mutations(cell,parms,i,j,k,im,jm,km)
#This function returns the new probability of mutation and increse the gene mutations    
#Genes that will mutate

p = exp( -( cell[i,j,k].N / parms.var7 )^2 );
#probm = probm/β1;
p = max(0.001,p); #savers
p = min(0.999,p); 
#gene_mut = min(unique(nbinrnd(3.1189,0.6560,[length(celula(1,1).genemutation),1])),length(celula(1,1).genemutation));

cell, numa = segregationlocal(cell,parms,i,j,k);
    
#numa = rand(1:10);
for il=1:parms.var12   # all genes can change at the same time
 
    pp = ( p + cell[i,j,k].segregation[il] - p* cell[i,j,k].segregation[il] );
    gene_mut = min.(rand( NegativeBinomial(numa[1],1-p),parms.var12),parms.var12);
    
     num_p = 0;
     num_g = 0;
     #isempty(gene_mut(il))==0 &&
    if (  isempty(gene_mut[il])==0 && gene_mut[il]>=1 && pp>rand() ) #p>probm # nut<θdiv*sqrt(-log(probm))
        # Tau-Leaping algorithm
        tau = -log(rand());
        num_p = rand(Poisson(tau));
        num_g = rand(Geometric(p));    
    end
        
    cell[i+im,j+jm,k+km].genemutation[il] = cell[i,j,k].genemutation[il] + num_p*num_g;
    cell[i+im,j+jm,k+km].thresholdexpression[il] = cell[i,j,k].thresholdexpression[il];               

end        
#los umbrales eta dependeran  de las etas y las mutaciones de la hija
#cell[i+im,j+jm].geneexpression = expression(cell,i,j,im,jm);
   cell = expression(cell,parms,i,j,k,im,jm,km);
#;.+ 0*rand(Poisson.(0.0005),gens).*rand(Geometric(0.5),gens); # + is driver 
   
return cell        
        
        
end    