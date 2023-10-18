function boxcount(c)   
    
    if ndims(c) == 3
        if size(c,3)==3 && size(c,1)>= 8 && size(c,1)>=8
            c = sum(c,3);
        end
    end
    
    c = Int64.(BitArray(sign.(c)))
    
    dim = ndims(c);
    
    if ndims(c) > 3 
        println("Error: Maximum dimensions is 3")
        return 0,0;
    end
    
    
    if length(c) == size(c,1)
        dim = 1;
        if size(c,1) !=1   
            c = c';
        end   
    end
    
    width = maximum(size(c));
    p     = log(width)/log(2);
            
    if p != round(p) || any(size(c) != width)
        p = Int32(ceil(p));
        width = 2^p; 
        if dim == 1
            mz = zeros(Int32,Int32(width))
            mz[1:length(c)] = c
            c = mz;
        elseif dim == 2
            mz = zeros(Int32,width, width)
            mz[1:size(c,1), 1:size(c,2)] = c;
            c = mz;
        elseif dim == 3
            mz = zeros(Int32,width,width,width);
            mz[1:size(c,1),1:size(c,2),1:size(c,3)] = c;
            c = mz;
        end
    end
    
    n = zeros(1,p+1);

    if dim == 1
        n[p+1] = sum(c)
        for g = (p-1):-1:0
            siz = 2^(p-g);
            siz2 = round(Int32,siz/2);
            for i=1:siz:(width-siz+1)
                c[i] = Int32(Bool(c[i]) || Bool(c[i]));
            end
            n[g+1] = sum(c[1:siz:(width-siz+1)]);
        end
        
    elseif dim == 2
        n[p+1] = sum(c);
        for g=(p-1):-1:0
            siz = 2^(p-g);
            siz2 = round(Int32,siz/2);
            for i=1:siz:(width-siz+1)
                for j=1:siz:(width-siz+1)
                    c[i,j] = Int32( Bool(c[i,j]) || Bool(c[i+siz2,j]) || Bool(c[i,j+siz2]) || Bool(c[i+siz2,j+siz2]) );
                end
            end
            n[g+1] = sum(c[1:siz:(width-siz+1),1:siz:(width-siz+1)]);
        end
    elseif dim == 3
        n[p+1] = sum(c);
        for g=(p-1):-1:0
            siz = 2^(p-g);
            siz2 = round(Int32,siz/2);
            for i=1:siz:(width-siz+1)
                for j=1:siz:(width-siz+1)
                    for k=1:siz:(width-siz+1)
                        c[i,j,k] = Int32( Bool(c[i,j,k]) || Bool(c[i+siz2,j,k]) || Bool(c[i,j+siz2,k]) || Bool(c[i+siz2,j+siz2,k]) || Bool(c[i,j,k+siz2]) || Bool(c[i+siz2,j,k+siz2]) || Bool(c[i,j+siz2,k+siz2])  || Bool(c[i+siz2,j+siz2,k+siz2]));
                    end
                end
            end
            n[g+1] = sum(c[1:siz:(width-siz+1),1:siz:(width-siz+1),1:siz:(width-siz+1)]);
        end
    end
    
    n = n[end:-1:1];
    r = 2.0.^(0:p); 

    m1 = log.(n)
    m2 = log.(r)
    
    m  = -(m1[2:end]-m1[1:end-1]) ./ (m2[2:end]-m2[1:end-1])  
    
    return m
end