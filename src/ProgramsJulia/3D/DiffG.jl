function Diffglobal(cell,parms,timm)

    h1,h2,h3 = size(cell);
   FV = 0.3;#Fracvolume(cell,parms); 
    # Boundary conditions
    for ii = 1:h3
        for jj = 1:h2
            for kk = 1:h3
                cell[1,jj,kk].N = 0.0;
                cell[h1,jj,kk].N = 0.0;
                cell[1,jj,kk].M = 0.0;
                cell[h1,jj,kk].M = 0.0;
                cell[1,jj,kk].E = 0.0;
                cell[h1,jj,kk].E = 0.0;
           
                cell[ii,1,kk].N = 0.0;
                cell[ii,h2,kk].N = 0.0;
                cell[ii,1,kk].M = 0.0;
                cell[11,h2,kk].M = 0.0;
                cell[ii,1,kk].E = 0.0;
                cell[ii,h2,kk].E = 0.0;

                cell[ii,jj,1].N = FV;
                cell[ii,jj,h3].N = 0.0;
                cell[ii,jj,1].M = FV;
                cell[ii,jj,h3].M = 0.0;
                cell[ii,jj,1].E = FV;
                cell[ii,jj,h3].E = 0.0;
            end
        end
    end
     
    for it = 1:timm
        
        for ii = 2:h1 -1
            for jj = 2:h2 -1
                for kk = 2:h3-1
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