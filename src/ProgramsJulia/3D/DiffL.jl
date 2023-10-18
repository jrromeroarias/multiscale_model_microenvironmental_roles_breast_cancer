function  Difflocal(cell,parms,i,j,k,timm)

#Local changes to nutrient concentration: Solve steady-state diffusion equations
max_it = timm;  #number of iterations
h1, h2, h3 = size(cell);

row1 = max.(i-10,1);  #take an area around the change
row2 = min.(i+10,h1);
col1 = max.(j-10,1);
col2 = min.(j+10,h2);
hig1 = max.(k-10,1);
hig2 = min.(k+10,h3); 
Lrow = row2-row1+1;
Lcol = col2-col1+1;
Lhig = hig2-hig1+1;

h1 = Lrow;
h2 = Lcol;
h3 = Lhig;

for it = 1:max_it
    
    for ii = row1+1:row2-1
        for jj = col1+1:col2-1
            for kk = hig1+1:hig2-1
                    #glusose
                    cell[ii,jj,kk].N = ( (cell[ii+1,jj,kk].N + cell[ii-1,jj,kk].N) + (cell[ii,jj+1,kk].N + cell[ii,jj-1,kk].N) + (cell[ii,jj,kk+1].N + cell[ii,jj,kk-1].N))./6.0 - parms.var1^2*cell[ii,jj,kk].N*cell[ii,jj,kk].Type[1]-parms.var2*parms.var1^2*cell[ii,jj,kk].N*cell[ii,jj,kk].Type[2];
                    cell[ii,jj,kk].N = max.(cell[ii,jj,kk].N,0); cell[ii,jj,kk].N = min.(cell[ii,jj,kk].N,1);
                    #oxygen
                    cell[ii,jj,kk].M = ( (cell[ii+1,jj,kk].M + cell[ii-1,jj,kk].M) + (cell[ii,jj+1,kk].M + cell[ii,jj-1,kk].M) + (cell[ii,jj,kk+1].M + cell[ii,jj,kk-1].M))./6.0 - parms.var3^2*cell[ii,jj,kk].M*cell[ii,jj,kk].Type[1]-parms.var4*parms.var3^2*cell[ii,jj,kk].M*cell[ii,jj,kk].Type[2];
                    cell[ii,jj,kk].M = max.(cell[ii,jj,kk].M,0); cell[ii,jj,kk].M = min.(cell[ii,jj,kk].M,1);
                    #estrogens
                    cell[ii,jj,kk].E = ( (cell[ii+1,jj,kk].E + cell[ii-1,jj,kk].E) + (cell[ii,jj+1,kk].E + cell[ii,jj-1,kk].E) + (cell[ii,jj,kk+1].E + cell[ii,jj,kk-1].E))./6.0 - parms.var5^2*cell[ii,jj,kk].E*cell[ii,jj,kk].Type[1]-parms.var6*parms.var5^2*cell[ii,jj,kk].E*cell[ii,jj,kk].Type[2];
                    cell[ii,jj,kk].E = max.(cell[ii,jj,kk].E,0); cell[ii,jj,kk].E = min.(cell[ii,jj,kk].E,1);
            end
        end
    end
end
    
return cell
                                            
end

