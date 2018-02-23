%...  The MatMol Group (2016)
    function [uxx uyy] = second_order_derivatives_2D(u,nx,ny,D2x,D2y)
%... use finite differences in 1D along x and along y
    v = reshape(u,nx,[]);
    uxx = reshape(D2x*v,nx*ny,1);
    uyy = reshape(v*D2y',nx*ny,1);
