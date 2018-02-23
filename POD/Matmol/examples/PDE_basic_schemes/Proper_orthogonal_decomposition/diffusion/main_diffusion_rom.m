%... The MatMol Group (2016)
%... The diffusion problem solved using the POD

    clear all
    clc

%... Spatial coordinates
    nd = 101;                 %... Discretization points
    z  = linspace(0,1,nd)';   

%... Parameters
    k = 0.1;

%... FEM matrices
    [MM, DM] = matfem (z, 'neu', 'neu', 'MM', 'DM');
    iMM      = eye(size(MM))/MM;
    lapop    = -iMM*DM;

%... Load the POD data (previously computed using file podcomputation.m)
    load podbasis_101
    neig = 8;
    phi  = pods.phi(:,1:neig);

%... Initial conditions
    x0 = 5*(z.^2/2 - z.^4/4) + 1;
    mx0 = phi'*MM*x0;

%... Projection of the laplacian operator
    mlapop = phi'*MM*lapop*phi;

%... Time span
    tlist = 0:0.1:4;

%... ODE integration
    options        = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [tout, mx_num] = ode15s(@pde_dif_pod, tlist, mx0, options, mlapop, k);   

%... Recovery the field
    x_num = phi*mx_num';

%... Plot the results
    mesh(tout, z, x_num)
