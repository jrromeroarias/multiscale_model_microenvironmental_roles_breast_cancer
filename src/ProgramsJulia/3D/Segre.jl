function  segregationlocal(cell,parms,i,j,k)

#Local changes for segregation index
h1, h2, h3 = size(cell);
rdio = 10;#parms.var16;
row1 = max.(i-parms.var16,1);  #take an area around the change
row2 = min.(i+parms.var16,h1);
col1 = max.(j-parms.var16,1);
col2 = min.(j+parms.var16,h2);
hig1 = max.(k-parms.var16,1);
hig2 = min.(k+parms.var16,h3);
Lrow = row2-row1+1;
Lcol = col2-col1+1;
Lhig = hig2-hig1+1;

num_a = zeros(parms.var12);
num_aa = zeros(parms.var12);

for il = 1:parms.var12
    
    for ii = row1+1:row2-1
        for jj = col1+1:col2-1
            for kk = hig1+1:hig2-1     
                num_a[il] = num_a[il] + parms.var9*((cell[ii,jj,kk].E/( cell[ii,jj,kk].E + 0.5 ))*sign.( abs.( cell[ii,jj,kk].genemutation[il] - cell[i,j,k].genemutation[il] ) )); 
                num_aa[il] = num_aa[il] + cell[ii,jj,kk].genemutation[il]; #sum para saber cuantas mutaciones por gen 
            end
        end
    end
    
    cell[i,j,k].segregation[il] = num_a[il]/Lrow/Lcol/Lhig; 
    
    cell[i,j,k].thresholdexpression[il] = max.(parms.var13 - parms.var14*cell[i,j,k].segregation[il] - parms.var15*cell[i,j,k].thresholdexpression[il], 0.05);     

end

num_b = Tuple.(findall(num_a.==maximum(num_a)));
num_c = num_a[first.(num_b)];

    if (num_c[1]==0)
        num_c[1] =1;
    end
    
    return cell, num_c

end