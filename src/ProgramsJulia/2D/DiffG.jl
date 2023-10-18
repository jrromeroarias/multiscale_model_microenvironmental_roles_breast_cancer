function Diffglobal(cell,parms,timm)

    h1,h2 = size(cell);
   FV = Fracvolume(cell,parms); 
    # Boundary conditions
    for ii = 1:h1
        cell[ii,1].N = 0.0;
        cell[ii,h2].N = 0.0;
        cell[ii,1].M = 0.0;
        cell[ii,h2].M = 0.0;
        cell[ii,1].E = 0.0;
        cell[ii,h2].E = 0.0;
    end
    for jj = 1:h2
        cell[1,jj].N = FV;
        cell[h1,jj].N = 0.0;
        cell[1,jj].M = FV;
        cell[h1,jj].M = 0.0;
        cell[1,jj].E = FV;
        cell[h1,jj].E = 0.0;
    end
     
    for it = 1:timm
        
        for ii = 2:h1 -1
            for jj = 2:h2 -1
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