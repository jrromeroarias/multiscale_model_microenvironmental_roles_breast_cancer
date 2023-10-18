function death(cell,parms,i,j)

    cell[i,j].Type = [0,0,0];
    cell[i,j].R = 0;
    cell[i,j].geneexpression[:] = zeros(parms.var12+1);
    cell[i,j].genemutation[:] = zeros(parms.var12);
    cell[i,j].thresholdexpression[:] = zeros(parms.var12);
    
return cell
end