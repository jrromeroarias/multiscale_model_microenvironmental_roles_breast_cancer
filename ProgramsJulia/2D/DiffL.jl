function  Difflocal(cell,parms,i,j,timm)

#Local changes to nutrient concentration: Solve steady-state diffusion equations
max_it = timm;  #number of iterations
h1, h2 = size(cell);

row1 = max.(i-10,1);  #take an area around the change
row2 = min.(i+10,h1);
col1 = max.(j-10,1);
col2 = min.(j+10,h2);
Lrow = row2-row1+1;
Lcol = col2-col1+1;

h1 = Lrow;
h2 = Lcol;

for it = 1:max_it
    
    for ii = row1+1:row2-1
        for jj = col1+1:col2-1
                 #glusose
                cell[ii,jj].N = (h2^2*(cell[ii+1,jj].N + cell[ii-1,jj].N) + h1^2*(cell[ii,jj+1].N + cell[ii,jj-1].N))./(2*(h1^2+h2^2)) - parms.var1^2*cell[ii,jj].N*cell[ii,jj].Type[1]-parms.var2*parms.var1^2*cell[ii,jj].N*cell[ii,jj].Type[2];
                cell[ii,jj].N = max.(cell[ii,jj].N,0); cell[ii,jj].N = min.(cell[ii,jj].N,1);
                #oxygen
                cell[ii,jj].M = (h2^2*(cell[ii+1,jj].M + cell[ii-1,jj].M) + h1^2*(cell[ii,jj+1].M + cell[ii,jj-1].M))./(2*(h1^2+h2^2)) - parms.var3^2*cell[ii,jj].M*cell[ii,jj].Type[1]-parms.var4*parms.var3^2*cell[ii,jj].M*cell[ii,jj].Type[2];
                cell[ii,jj].M = max.(cell[ii,jj].M,0); cell[ii,jj].M = min.(cell[ii,jj].M,1);
                #estrogens
                cell[ii,jj].E = (h2^2*(cell[ii+1,jj].E + cell[ii-1,jj].E) + h1^2*(cell[ii,jj+1].E + cell[ii,jj-1].E))./(2*(h1^2+h2^2)) - parms.var5^2*cell[ii,jj].E*cell[ii,jj].Type[1]-parms.var6*parms.var5^2*cell[ii,jj].E*cell[ii,jj].Type[2];
                cell[ii,jj].E = max.(cell[ii,jj].E,0); cell[ii,jj].E = min.(cell[ii,jj].E,1);
        end
    end
end
    
return cell
                                            
end

