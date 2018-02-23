%... The MatMol Group (2016)
     function ut = generalized_burgers_pde(t,u)
%...
%... set global variables
    global eps nu 
    global nx ny nv x y D1x_up D1x_c3 D1x_c5 D1y_up D1y_c3 D1y_c5
    t
%...
%... boundary conditions
    u(1:nx,1) = 0;          % BC along bottom boundary
    u(nv-nx+1:nv,1) = 0;    % BC along upper boundary
    u(1:nx:nv-nx+1,1) = 0;  % BC along left boundary
    u(nx:nx:nv,1) = 0;      % BC along right boundary
%...
%... evolution of nu
%...
    for i = 1:nv
        if abs(u(i))> 0.5
            nu(i) = 1;
        else
            nu(i) = 0;
        end
    end
%...
%... spatial derivatives
%...
%...   computation of (u^2)x and (u^2)y :
%...
    [fx fy] = KT_scheme_2D(u,x,y,nx,ny,t);
%    [fx fy] = first_order_derivatives_2D(u.*u,nx,ny,D1x_c3,D1y_c3);
%...
%...   computation of (nu*ux)x and (nu*uy)y :
%...
    [ux uy] = first_order_derivatives_2D(u,nx,ny,D1x_c3,D1y_c3);
    nuux = nu.*ux;
    nuuy = nu.*uy;
    [nu_uxx nonutil] = first_order_derivatives_2D(nuux,nx,ny,D1x_c3,D1y_c3);
    [nonutil nu_uyy] = first_order_derivatives_2D(nuuy,nx,ny,D1x_c3,D1y_c3);
%...
%... temporal derivatives
    ut = - fx - fy + eps*(nu_uxx+nu_uyy);
