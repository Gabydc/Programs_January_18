%... The MatMol Group (2016)
%... Solution of brusselator problem using the POD
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
    nd = 501;                 %... Discretization points
    z  = linspace(0,1,nd)';   

%... Parameters
    a     = 1;
    b     = 3;
    alpha = 0.02;

%... FEM matrices
    [MM, D2_int] = matfem (z, 'dir', 'dir', 'MM', 'DM', 'BM');
    iMM      = eye(size(MM))/MM;

%... POD data (previously computed using file podcomputation.m)
    load pod_data_brus
    neig_u = 6;
    neig_v = 4;
    phiu   = pods_u.phi(:,1:neig_u);
    phiv   = pods_v.phi(:,1:neig_v);

%... Projections
    A_u  = -alpha*phiu'*D2_int*phiu;
    A_v  = -alpha*phiv'*D2_int*phiv;
    proj_op_u = phiu'*MM;
    proj_op_v = phiv'*MM;

%... Initial conditions
    u0  = 1+sin(2*pi*z);
    mu0 = proj_op_u*u0;
    v0  = 3*ones(nd,1);
    mv0 = proj_op_v*v0;
    mx0 = [mu0; mv0];

%... Boundary conditions
    q     = 1e3;
    u_inf = 1;
    v_inf = 3;

%... Time span
    tlist = 0:0.1:20;

%... ODE integration
    options    = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [t, m_num] = ode15s(@pde_brusselator_pod, tlist, mx0, options,...
                        neig_u, neig_v, a, b, q, phiu, phiv,...
                        proj_op_u, proj_op_v, A_u, A_v, nd, u_inf,...
                        v_inf);   

%... Field separation
    mu   = m_num(:,1:neig_u)';
    mv   = m_num(:,1+neig_u:neig_u+neig_v)';
    upod = phiu*mu;
    vpod = phiv*mv;

%... Plot the results
    mesh(z, t, upod')
    figure
    mesh(z, t, vpod')
 