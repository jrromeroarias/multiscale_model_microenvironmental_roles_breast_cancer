function  reproduction(cell,parms,i,j)

na = rand(1:5);
aux = 0; 
                
    if (cell[i+1,j].Type[2]>=1 && cell[i-1,j].Type[2]>=1 && cell[i,j+1].Type[2]>=1 && cell[i,j-1].Type[2]>=1 && na ==5 && aux ==0)
        cell[i,j].Type = [0,1,0];  #[normal cancer death]
        cell[i,j].R = 0;
        cell = mutations(cell,parms,i,j,0,0);
        aux = 1;
    elseif (cell[i+1,j].Type[2]==0 && na ==1 && aux ==0 )
        cell[i+1,j].Type = [0,1,0];
        cell[i,j].R = 0;
        cell[i+1,j].R = 0;
        cell = mutations(cell,parms,i,j,1,0);
        aux = 1;
    elseif (cell[i-1,j].Type[2]==0 && na ==2 && aux ==0 )
        cell[i-1,j].Type = [0,1,0];
        cell[i,j].R = 0;
        cell[i-1,j].R = 0;
        cell = mutations(cell,parms,i,j,-1,0);
        aux = 1;
    elseif (cell[i,j+1].Type[2]==0 && na ==3 && aux ==0 )
        cell[i,j+1].Type = [0,1,0];
        cell[i,j].R = 0;
        cell[i,j+1].R = 0;
        cell = mutations(cell,parms,i,j,0,1);
        aux = 1;
    elseif (cell[i,j-1].Type[2]==0 && na ==4 && aux ==0 )
        cell[i,j-1].Type = [0,1,0];
        cell[i,j].R = 0;
        cell[i,j-1].R = 0;
        cell = mutations(cell,parms,i,j,0,-1);
        aux = 1;
    end
                                                        
    
return cell                                                
    
end

