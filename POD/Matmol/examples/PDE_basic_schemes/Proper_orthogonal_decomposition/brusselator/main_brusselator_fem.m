%... The MatMol Group (2016)
%... Solution of brusselator problem using the FEM
%...
%... Brusselator
%...
%...           2
%... u  = A + u v - (B+1)u + alpha*u
%...  t                             zz
%...
%...            2
%... v  = Bu - u v + alpha*v
%...  t                     zz
%...
%...
%... where
%...
%...    A = 1
%...
%...    B = 3
%...  
%...    alpha = 1/50
%...
%...    0 < z <1
%...
%...    t > 0
%...
%... 
%... with
%...
%...    C.L. :  u(0,t) = u(1,t) = 1
%...
%...            v(0,t) = v(1,t) = 3
%...
%...    C.I. :  u(z,0) = 1 + sin(2*pi*z)
%...
%...            v(z,0) = 3

    clear all
    clc

%... Spatial coordinates
    nd = 71;                 %... Discretization points
    z  = linspace(0,1,nd)';   

%... Parameters
    A     = 1;
    B     = 3;
    alpha = 0.02;

%... FEM matrices
    [MM, DM, BM] = matfem (z, 'dir', 'dir', 'MM', 'DM', 'BM');
    iMM      = eye(size(MM))/MM;
    lapop    = -iMM*(alpha*DM);

%... Initial conditions
    u0 = 1+sin(2*pi*z);
    v0 = 3*ones(nd,1);
    x0 = [u0; v0];

%... Time span
    tlist = 0:0.1:20;

%... ODE integration
    options       = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [tout, x_num] = ode15s(@pde_brusselator_fem, tlist, x0, options,...
                       lapop, nd, A, B, BM(1,1), iMM);   

%... Field separation
    u = x_num(:,1:nd)';
    v = x_num(:,1+nd:2*nd)';

%... Plot the results
    mesh(z, tout, u')
    figure
    mesh(z, tout, v')