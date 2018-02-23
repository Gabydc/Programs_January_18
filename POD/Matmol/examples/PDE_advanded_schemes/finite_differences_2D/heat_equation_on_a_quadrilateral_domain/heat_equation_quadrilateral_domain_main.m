%... The MatMol Group (2016)
%...   Heat equation on a convex quadrilateral domain
%...
%... The equation in 2D is given by
%...
%...   Tt = k[Txx + Tyy] 
%...          
%... It is defined on a quadrilater defined by 4 corners. Their positions
%... can be freely chosen and are stored (in whatever order) in the matrix
%... corners(4,2). The function grid will order these points anticlockwise
%... starting with corner S1, which has the smallest x-coordinate. If two
%... points are candidates (they have the same small x-coordinate), then the
%... one with the smallest y-coordinate is selected.
%...
%... with the conditions:
%...
%...   IC = T(x,y,0) = 0
%...
%...   BCs : Dirichlet : T = 0 along S1-S2 and S1-S4
%...
%...                     T = TM along S2-S3 and S3-S4
%...
%... The problem paramters are
%...   k = 10
%...
%...   TM = 100
%...
%... the following code computes a solution to the 2D heat equation
%...
    close all
    clear all
%...
%... start a stopwatch timer
    tic
%...
%... set global variables
    global k TM 
    global nv nksi neta ksi eta D1_ksi D1_eta D2_ksi D2_eta 
    global ksi_x ksi_y eta_x eta_y ksi_xx ksi_yy eta_xx eta_yy
%...
%... model parameters
    TM = 100;
    k = 10;
%...
%... spatial domain 
    corners = [ -3 -2 ; 2 -5 ; 1 4 ; -2 1 ];
%...
%... coordinates ksi and eta in a square domain [0 1]x[0 1]
    nksi = 101; neta = 101%
    nv = nksi*neta;
    dksi = 1/(nksi-1); ksi = [0:dksi:1]';
    deta = 1/(neta-1); eta = [0:deta:1]';
%...
%... spatial grid on the actual domain: 
%...      * nksi points are used on the edge S1-S2
%...      * neta points are used on the edge S1-S4
%...      * (x,y) are the actual coordinates of the grid points
%...      * rho and eta are the angles of the edge S1-S2, and S1-S4,
%...        respectively
%...      * ksi_x, ksi_y and eta_x, eta_y are the transformations to convert
%...        du/dksi and du/deta into du/dx and du/dy
%...      * the same holds for the second-order derivatives
%...
    [x y ksi_x ksi_y eta_x eta_y ksi_xx ksi_yy eta_xx eta_yy] = actual_grid(corners,nksi,neta,ksi,eta);
%...
%... initial conditions
    T0(1:nv,1) = zeros(nv,1);
%...
%... finite difference (FD) approximation of the spatial derivatives
%...
    D1_ksi = three_point_centered_D1(ksi);
    D1_eta = three_point_centered_D1(eta);
    D2_ksi = five_point_centered_D2(ksi);
    D2_eta = five_point_centered_D2(eta);
%...  
%... call to ODE solver 
%...
    t = [0:.0005:.01];
%...
    [tout,Tout] = ode45(@heat_equation_quadrilateral_domain_pde,t,T0);
%...  
%... impose BCs 
%...
    for k = 1:length(tout)
        Tout(k,1:nksi) = 0;
        Tout(k,1:nksi:nksi*(neta-1)+1) = 0;
        Tout(k,nksi:nksi:nksi*neta) = TM;
        Tout(k,nksi*(neta-1)+1:nksi*neta)= TM;
    end
%...  
%... plot results 
%...
    for j = 1:length(tout)
        figure
        u(neta:-1:1,:) = reshape(Tout(j,1:nv),nksi,[])';
        surf(x,y,u,'EdgeColor','none')
        view([25 32]);
        axis([-3 2 -5 4 0 100])
%        camlight right ;
%        lighting gouraud ;
        xlabel('x'); ylabel('y'); zlabel('T(x,y)')
    end
%...     
%... read the stopwatch timer
    toc