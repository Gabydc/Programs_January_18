%... The Matmol group (2016)
    function uxy = mixed_second_order_derivatives_2D(u,nx,ny,D1x,D1y)

%... use finite differences in 1D along x and along y
    v = reshape(u,nx,[]);
    uxy = reshape((D1x*v)*D1y',nx*ny,1);
