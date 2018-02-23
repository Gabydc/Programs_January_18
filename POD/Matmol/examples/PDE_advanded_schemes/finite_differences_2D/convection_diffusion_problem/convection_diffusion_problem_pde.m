%... The MatMol Group (2016)
     function ut = convection_diffusion_problem_pde(t,u)
%...
%... set global variables
    global nu b1 b2 f 
    global nx ny nv D1x_up D1y_up D1x_c3 D1y_c3 D1x_c5 D1y_c5 D2x_c3 D2y_c3 D2x_c5 D2y_c5 
%...
%... boundary conditions
    u(1:nx,1) = 0;          % BC along bottom boundary
    u(nv-nx+1:nv,1) = 0;    % BC along upper boundary
    u(1:nx:nv-nx+1,1) = 0;  % BC along left boundary
    u(nx:nx:nv,1) = 0;      % BC along right boundary
%...
%... spatial derivatives
%...
    [ux uy] = first_order_derivatives_2D(u,nx,ny,D1x_up,D1y_up);
    [uxx uyy] = second_order_derivatives_2D(u,nx,ny,D2x_c5,D2y_c5);
%...
%... temporal derivatives
    ut = - b1*ux - b2*uy + nu*(uxx + uyy) + f;
