%... The MatMol Group (2016)
%...
%... Burgers equation
%...
%... The well known Burgers equation:
%...
%... x  = -x*x  + mu*x   
%...  t       z       zz
%...
%...
%... where
%...
%...  x      dependent variable
%...
%...  t       time
%...
%...  z       space
%...
%...  mu      diffusion coefficient
%...
%...
%... an exact solution is given by
%...
%... x(z,t) = burgers_exact(z,t)
%...
%... the boundary conditions (BCs) for Burgers' equation are taken as
%...
%... x(z=z0,t) = burgers_exact(z0,t)     (Dirichlet BCs)
%...
%... x(z=zL,t) = burgers_exact(zL,t)               
%...
%... where
%...
%...  zL      spatial domain
%...
%... the initial condition (IC) for Burgers' equations are taken as
%...
%... x(z,t=0) = burgers_exact(z,0)
%...
%... the following code computes a solution to Burgers' equation
%...

    close all
    clear all

%... Global variables (shared by other files)
    global mu
    global wquad xquad
    global n dz z0 zL D2

%... Finite element spatial grid
    z0  = 0;
    zL  = 1;
    n   = 201;
    nel = n-1;
    dz  = (zL - z0)/(n-1);
    z   = (z0:dz:zL)';

%... Second order differentiation matrix
    D2 = lagrlinD2(dz,n);

%... Problem parameters
    mu = 0.01;

%... Computation of the weigths and the abcissa used in the numerical 
%... integration via the Gauss-Legendre quadrature
    nquad     = 2;
    beta      = 0.5./sqrt(1-(2*(1:nquad)).^(-2));
    T         = diag(beta,1) +diag(beta,-1);
    [V,D]     = eig(T);
    xquad     = diag(D);
    [xquad,i] = sort(xquad);
    wquad     = 2*V(1,i).^2;

%... Initial conditions
    x = zeros(1,n);
    for ii = 1:n
        x(ii) = burgers_exact(z(ii),0);
    end

%... Time instants at which the IVP solver will save the solution
    dt   = 0.1;
    time = (0:dt:1);
    nt   = length(time);

%... Time integration with the IVP solver
    ne      = 1;
    M       = lagrlinmass(dz,n,ne); %... FEM mass matrix
    options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',1e-6);
    [timeout,yout]= ode15s(@burgers_pdes,time,x,options);


%... Plot the results
    figure
    hold
    plot(z,yout,'.-k')
    yexact=zeros(n,length(time));
    for k=1:length(timeout)
        for i=1:n
            yexact(i,k)=burgers_exact(z(i),time(k));
        end
    end
    plot(z,yexact,'r')
