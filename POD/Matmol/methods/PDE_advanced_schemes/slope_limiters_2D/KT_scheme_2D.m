%...  The Matmol group (2016)
     function [fx fy] = KT_scheme_2D(u,x,y,nx,ny,t)
%...
    fkurgx = zeros(nx,ny);
    fkurgy = zeros(ny,nx);
    usquarex = reshape(u,nx,[]);
    for i = 1:ny
        fkurgx(:,i) = KT_centered_slope_limiter_fz(1,nx,x,t,usquarex(:,i),@flux,@dflux);
    end
    fx = reshape(fkurgx,nx*ny,1);
    usquarey = usquarex';
    for i = 1:nx
        fkurgy(:,i) = KT_centered_slope_limiter_fz(1,ny,y,t,usquarey(:,i),@fluy,@dfluy);
    end
    fy = reshape(fkurgy',nx*ny,1);
