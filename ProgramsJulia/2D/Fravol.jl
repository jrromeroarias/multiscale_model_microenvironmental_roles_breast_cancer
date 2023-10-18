function Fracvolume(cell,parms)

    h1,h2 = size(cell);

TcellNp = zeros(Nx,Ny);
        for ii = 2:h1 -1
            for jj = 2:h2 -1
                #glusose
                TcellNp[ii,jj] = cell[ii,jj].Type[2]; 
            end
        end

    fracvol = 1.0 - (sum(TcellNp)/(Nx*Ny));
    return fracvol

end