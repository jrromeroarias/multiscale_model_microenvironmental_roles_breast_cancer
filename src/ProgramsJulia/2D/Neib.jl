function neibs(cell,parms)
    
h1,h2 = size(cell);

NNcell = zeros(h1,h2);
Tcell = zeros(h1,h2);
RRa = zeros(h1,h2);

gamma = 2*pi/parms.var5^2/0.5/(72/6);  # cell cycle normal cells 72, mean E and 60 div of cancer cell by 1 normal cell #t~1hr
gammax = gamma*parms.var5^2;
    for ii =2:h1-1
        for jj = 2:h2-1             
            NNcell[ii,jj]= cell[ii,jj].Type[2]*( cell[ii+1,jj].Type[2] +  cell[ii-1,jj].Type[2] + cell[ii,jj+1].Type[2] + cell[ii,jj-1].Type[2] );          
            Tcell[ii,jj] = cell[ii,jj].Type[2];
            cell[ii,jj].R =  min.(cell[ii,jj].R + gammax*( parms.var6*cell[ii,jj].Type[2] + cell[ii,jj].Type[1])*cell[ii,jj].E, 36.0*cell[ii,jj].Type[2] + 72.0*cell[ii,jj].Type[1] );
            RRa[ii,jj] = cell[ii,jj].R;            
        end
    end

ind=  Tuple.(findall( (Tcell.==1) .& (NNcell.<4)  ));  
ip =first.(ind); jp =last.(ind);    

    
ipl = length(ip);
jpl = randperm(ipl);
ipp = ip[jpl];
jpp = jp[jpl];
    
  

return cell, ipp, jpp                                         
                                            
end
