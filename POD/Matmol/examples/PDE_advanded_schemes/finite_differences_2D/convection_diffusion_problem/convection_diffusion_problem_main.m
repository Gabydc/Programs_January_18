%... The MatMol Group (2016)
%...   Convection-diffusion problem
%...
%...   ut = -b1*ux - b2*uy + nu*(uxx + uyy) + f(x,y)
%...
%...   0 < x < 1     0 < y < 1
%...
%...   0 < t 
%...
%...   b1 = cos(theta)/2  b2 = sin(theta)/2 : - pi/2 < theta  < + pi/2
%...
%...   nu = 0.001
%...
%...   f(x,y) = 1
%...
%...   ICs :  u(x,y,0) = 0
%...
%...   BCs :  u(x,y,t) = 0 on the boundary
%...
%... the following code computes a solution to the convection-diffusion problem
%...
    close all
    clear all
%...
%... start a stopwatch timer
    tic
%...
%... set global variables
    global nu b1 b2 f 
    global nx ny nv D1x_up D1y_up D1x_c3 D1y_c3 D1x_c5 D1y_c5 D2x_c3 D2y_c3 D2x_c5 D2y_c5 
%...
%... spatial grid
    x0 = 0; y0 = 0;
    xL = 1; yL = 1;
    nx = 201; ny = 201;
    dx = (xL - x0)/(nx-1); dy = (yL - y0)/(ny-1);
    x = [x0:dx:xL]; y = [y0:dy:yL];
%...
%... model parameters
    nu = 1e-3;
%    theta = 0; 
%    theta = pi/2;
%    theta = pi/4;
    theta = - pi/4;
%    theta = -pi/3.25;
    b1 = cos(theta)/2;
    b2 = sin(theta)/2;
%...
%... initial conditions
    u0 = zeros(nx*ny,1);
    nv = length(u0);
    f = ones(nv,1);
%...
%... finite difference (FD) approximation of the spatial derivatives
%...
    D1x_up = two_point_upwind_D1(x,b1);%three_point_upwind_D1(x,b1);%four_point_biased_upwind_D1(x,b1);%five_point_biased_upwind_D1(x,b1)
    D1x_c3 = three_point_centered_D1(x);
    D1x_c5 = five_point_centered_D1(x);
    D2x_c3 = three_point_centered_D2(x);
    D2x_c5 = five_point_centered_D2(x);
%...
    D1y_up = two_point_upwind_D1(y,b2);%three_point_upwind_D1(y,b2);%four_point_biased_upwind_D1(y,b2);%five_point_biased_upwind_D1(y,b2)
    D1y_c3 = three_point_centered_D1(y);
    D1y_c5 = five_point_centered_D1(y);
    D2y_c3 = three_point_centered_D2(y);
    D2y_c5 = five_point_centered_D2(y);
%...  
%... call to ODE solver 
%...
    t = [0:.25:1];
    [tout,uout] = ode45(@convection_diffusion_problem_pde,t,u0);
%...  
%... impose BCs 
%...
    uout(:,1:nx) = 0;              % BC along y = 0
    uout(:,nv-nx+1:nv) = 0;        % BC along y = 1
    uout(:,1:nx:nv-nx+1) = 0;      % BC along x = 0
    uout(:,nx:nx:nv) = 0;          % BC along x = 1
%...  
%... plot results 
%...
    [X Y] = meshgrid(x,y(ny:-1:1));
%...
    figure
    for j = 2:length(tout)
%...    
    subplot(2,2,j-1)
%     figure
%...
    u(ny:-1:1,:) = reshape(uout(j,1:nv),nx,[])';
%...
    surf(X,Y,u,'FaceColor',[.95 .95 .95],'EdgeColor','none')
    hold
    axis([0 1.1 0 1 -.1 1.3])
    xlabel('x'); ylabel('y'); zlabel('u(x,y)');
    view([32 18]);
    camlight right ;
    lighting gouraud ;
end
%...     
%... read the stopwatch timer
    toc