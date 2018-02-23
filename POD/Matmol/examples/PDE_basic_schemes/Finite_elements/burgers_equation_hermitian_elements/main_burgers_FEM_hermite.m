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

%... Global variables
    global mu
    global wquad xquad
    global n dz z0 zL D2

%... Finite element spatial grid
    z0  = 0;
    zL  = 1;
    n   = 101;
    nel = n-1;
    dz  = (zL - z0)/(n-1);
    z   = linspace(z0,zL,n);

%... Second order differentiation matrix
    D2 = hermiteD2(dz,n);

%... Problem parameters
    mu = 0.001;

%... Computation of the weigths and the abcissa to be employed in the
%... numerical integration via the Gauss-Legendre quadrature
    nquad     = 2;
    beta      = 0.5./sqrt(1-(2*(1:nquad)).^(-2));
    T         = diag(beta,1) +diag(beta,-1);
    [V,D]     = eig(T);
    xquad     = diag(D);
    [xquad,i] = sort(xquad);
    wquad     = 2*V(1,i).^2;

%... Initial conditions
    x  = zeros(1,n);
    xz = zeros(1,n);
    for i = 1:n
        x(i) = burgers_exact(z(i),0);
        xz(i) = derburgers_exact(z(i),0);
    end
    xx(1:2:2*n-1)=x;
    xx(2:2:2*n)=dz*xz/2;

%... Time instants at which the IVP solver will save the solution
    dt   = 0.1;
    time = (0:dt:1);
    nt   = length(time);

%... Time integration with the IVP solver
    M       = hermitemass(dz,n); %... Finite element mass matrix
    options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',1e-6);
    [timeout, xout] = ode15s(@burgerspdes_hermite,time,xx,options);

%... Recovering the states from the ode15s output (xout)
    x_direct = xout(:,1:2:2*n-1);

%... Recovering the solution using the hermitian polynomials
    dzvis = 0.5;
    zvis = -1:dzvis:1;
    nvis = length(zvis);
%... Basis functions
    [phi1, phi2, phi3, phi4] = trialfunctionsher2n(zvis);
    z_absc   = linspace(z0,zL,(n-1)*(nvis-1)+1);

    for k=1:nel
        xx = xout(:,2*k-1)*phi1+xout(:,2*k)*phi2+...
             xout(:,2*k+1)*phi3+xout(:,2*k+2)*phi4; 
        if k == 1
            n1 = 1;
            n2 = 5;
            ordo(:,n1:n2) = xx;
        else
            n1 = n2+1;
            n2 = n2+4;
            ordo(:,n1:n2) = xx(:,2:end);
        end    
    end

%... Plot the results
    figure
    hold
    ztheor=(0:.001:1);
    ntheor=length(ztheor);
    for k=1:length(timeout)
        for i=1:ntheor
            yexact(k,i)= burgers_exact(ztheor(i),timeout(k));
        end
    end
    plot(z,x_direct,'-k')
    hold on
    plot(ztheor,yexact,'r')
    axis([0 1 0 1.4])
    hold off

    plot(z_absc,ordo,'-k')
    hold on
    plot(ztheor,yexact,'r')
    axis([0 1 0 1.4])
    hold off

