function death(cell,parms,i,j,k)

    cell[i,j,k].Type = [0,0,0];
    cell[i,j,k].R = 0;
    cell[i,j,k].geneexpression[:] = zeros(parms.var12+1);
    cell[i,j,k].genemutation[:] = zeros(parms.var12);
    cell[i,j,k].thresholdexpression[:] = zeros(parms.var12);
    
return cell
end