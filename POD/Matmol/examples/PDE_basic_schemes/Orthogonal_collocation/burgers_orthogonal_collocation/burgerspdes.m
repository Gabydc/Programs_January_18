%... The MatMol Group (2016)
    function xt = burgerspdes(t,u)

%... Global variables
    global mu
    global n z0 zL D2 u1

    u1 = u;

%... Linear part
    ut = mu*D2*u;

%... Nonlinear part
    y1 = feval('integrand',-1/sqrt(3));
    y2 = feval('integrand',1/sqrt(3));

%... Time derivative
    ut(2:2:2*n-2,1) = ut(2:2:2*n-2,1)-y1;
    ut(3:2:2*n-1,1) = ut(3:2:2*n-1,1)-y2;

%... Boundary conditions and construction of xt
    xt(1,1)       = u(1)- burgers_exact(z0,t);
    xt(2:2*n-1,1) = ut(2:2*n-1,1);
    xt(2*n,1)     = u(2*n-1) - burgers_exact(zL,t);


