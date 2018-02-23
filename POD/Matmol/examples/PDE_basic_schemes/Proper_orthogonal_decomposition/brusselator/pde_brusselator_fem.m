%... The MatMol Group (2016)
    function dx = pde_brusselator_fem(t, x, lapop, nd, A, B, q, iMM)

%... Field separation
    u = x(1:nd,1);
    v = x(1+nd:2*nd,1);

%... Boundary conditions
    Gu(1,1)  = q*(1-u(1));
    Gu(nd,1) = q*(1-u(nd));
    Gv(1,1)  = q*(3-v(1));
    Gv(nd,1) = q*(3-v(nd));
    iGu      = iMM*Gu;
    iGv      = iMM*Gv;

%... PDE
    du = lapop*u + A - (B + 1)*u + u.^2.*v + iGu;
    dv = lapop*v + B*u - u.^2.*v + iGv;
    dx = [du; dv];