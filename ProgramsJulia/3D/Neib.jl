function neibs(cell,parms)
    
h1,h2,h3 = size(cell);

NNcell = zeros(h1,h2,h3);
Tcell = zeros(h1,h2,h3);
RRa = zeros(h1,h2,h3);

gamma = 2*pi/parms.var5^2/0.5/(72/6);  # cell cycle normal cells 72,mean E and 60 div of cancer cell by 1 normal cell #t~1hr
gammax = gamma*parms.var5^2;
    for ii =2:h1-1
        for jj = 2:h2-1
            for kk = 2:h3-1             
                NNcell[ii,jj,kk]= cell[ii,jj,kk].Type[2]*( cell[ii+1,jj,kk].Type[2] +  cell[ii-1,jj,kk].Type[2] + cell[ii,jj+1,kk].Type[2] + cell[ii,jj-1,kk].Type[2] + cell[ii,jj,kk+1].Type[2] + cell[ii,jj,kk-1].Type[2]);          
                Tcell[ii,jj,kk] = cell[ii,jj,kk].Type[2];
                cell[ii,jj,kk].R =  min.(cell[ii,jj,kk].R + gammax*( parms.var6*cell[ii,jj,kk].Type[2] + cell[ii,jj,kk].Type[1])*cell[ii,jj,kk].E, 36.0*cell[ii,jj,kk].Type[2] + 72.0*cell[ii,jj,kk].Type[1] );
                RRa[ii,jj,kk] = cell[ii,jj,kk].R; 
            end           
        end
    end

#ind=  Tuple.(findall( (Tcell.==1) .& (NNcell.<=4) .& (RRa.> 35) )); 
ind=  Tuple.(findall( (Tcell.==1) .& (NNcell.<6)  ));  
#ind=  Tuple.(findall( (Tcell.==1) )); 
ip = getindex.(ind, 1) ;  jp = getindex.(ind, 2); kp = getindex.(ind, 3);   

ipl = length(ip);
jpl = randperm(ipl);
ipp = ip[jpl];
jpp = jp[jpl];
kpp = kp[jpl];    
  
#return ipp, jpp, kpp     
return cell, ipp, jpp, kpp                                         
                                            
end
