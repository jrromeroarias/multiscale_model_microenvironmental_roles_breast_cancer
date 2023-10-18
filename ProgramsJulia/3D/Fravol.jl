function Fracvolume(cell,parms)

    h1,h2,h3 = size(cell);

TcellNp = zeros(Nx,Ny,Nz);
        for ii = 2:h1 -1
            for jj = 2:h2 -1
                for kk = 2:h3-1
                    #glusose
                    TcellNp[ii,jj,kk] = cell[ii,jj,kk].Type[2];
                end
            end
        end

    fracvol = 1.0 - (sum(TcellNp)/(Nx*Ny*Nz));
    return fracvol

end