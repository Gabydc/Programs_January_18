%... The Matmol group (2016)
    function [ux uy] = first_order_derivatives_2D(u,nx,ny,D1x,D1y)
%... use finite differences in 1D along x and along y
    v = reshape(u,nx,[]);
    ux = reshape(D1x*v,nx*ny,1);
    uy = reshape(v*D1y',nx*ny,1);
