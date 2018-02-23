%...  The MatMol Group (2016)
    function dfxdx = dflux_dx(n,ne,x)

%... Computes the flux gradient 

    global v     
    dfxdx(1:n,1,1) = v;
    dfxdx(1:n,1,2) = 0;
    dfxdx(1:n,2,1) = 0;
    dfxdx(1:n,2,2) = 0;

    end 