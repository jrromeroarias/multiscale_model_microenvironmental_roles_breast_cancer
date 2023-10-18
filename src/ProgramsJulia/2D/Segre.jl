function  segregationlocal(cell,parms,i,j)

#Local changes for segregation index
h1, h2 = size(cell);
rdio = 10;#parms.var16;
row1 = max.(i-parms.var16,1);  #take an area around the change
row2 = min.(i+parms.var16,h1);
col1 = max.(j-parms.var16,1);
col2 = min.(j+parms.var16,h2);
Lrow = row2-row1+1;
Lcol = col2-col1+1;

num_a = zeros(parms.var12);
num_aa = zeros(parms.var12);

for il = 1:parms.var12
    
    for ii = row1+1:row2-1
        for jj = col1+1:col2-1     
            
            num_a[il] = num_a[il] + parms.var9*((cell[ii,jj].E/( cell[ii,jj].E + 0.5 ))*sign.( abs.( cell[ii,jj].genemutation[il] - cell[i,j].genemutation[il] ) )); 
            
            num_aa[il] = num_aa[il] + cell[ii,jj].genemutation[il]; #sum para saber cuantas mutaciones por gen 
                
        end
    end
    
    cell[i,j].segregation[il] = num_a[il]/Lrow/Lcol; 
    
    
    cell[i,j].thresholdexpression[il] = max.(parms.var13 - parms.var14*cell[i,j].segregation[il] - parms.var15*cell[i,j].thresholdexpression[il], 0.05);     

end

num_b = Tuple.(findall(num_a.==maximum(num_a)));
num_c = num_a[first.(num_b)];

    if (num_c[1]==0)
        num_c[1] =1;
    end
    
    return cell, num_c

end