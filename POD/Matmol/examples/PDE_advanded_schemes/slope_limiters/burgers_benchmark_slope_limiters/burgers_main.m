%...  The MatMol Group (2016)
%... The well known Burgers equation:
%...
%... x  = - (0.5*x^2)  + mu*x
%...  t              z       zz
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
%... the boundary conditions (BCs) for Burgers' equation are
%... taken as 
%... 
%... x(z=z0,t) = burgers_exact(z0,t)     (Dirichlet BCs)
%...
%... x(z=zL,t) = burgers_exact(zL,t)
%...
%... where
%...
%...  zL      spatial domain
%...
%... the initial condition (IC) for Burgers' equations are
%... taken as 
%...
%... x(z,t=0) = burgers_exact(z,0)
%...
%... the following code computes a solution to Burgers'
%... equation 

    close all
    clear all

%... start a stopwatch timer
    tic

%... set global variables
    global mu;
    global ne z0 zL z n method D2;

%... spatial grid
    z0 = 0.0;
    zL = 1.0;
    ne = 1;
    n  = 201;
    dz = (zL-z0)/(n-1);
    z  = [z0:dz:zL]';

%... model parameter
     mu=0.0005;
%    mu = 0.000001;

%... initial conditions
    for i = 1:n,
        x(i) = burgers_exact_2(z(i),0);
    end;

%... select slope limiter approximation
    %  method = 'Van_Leer'       %... Van Leer slope limiter
    %  method = 'superbee'       % Superbee slope limiter
    %  method = 'smart'          % Smart slope limiter
    %  method = 'mc'             % monotonized center limiter
    %  method = 'minmod'         % Minmod slope limiter
    % method = 'Koren'           % Koren slope limiter
      method = 'KT'             % Kurganov Tadmor centered scheme 

%... differentiation matrix (diffusive term - direct
%... differentiation) 
    D2 = three_point_centered_D2(z);

%... call to ODE solver
    t            = [0:0.1:1];
    options      = odeset('RelTol',1e-3,'AbsTol',1e-3,'Stats','on'); 
    [tout, xout] = ode15s('burgers_pde',t,x,options);

%... Plot the solution
    plot(z,xout,'k');
    xlabel('z');
    ylabel('x(z,t)');
%...     title('Burgers equation')
    hold on
%...
    for k=1:length(tout),
        for i=1:n,
            xexact(k,i) = burgers_exact_2(z(i),tout(k));
        end;
        plot(z,xexact,'r')
    end;
    axis([0 1 0 1.2]);
%...
    figure(3)
    surf(z,tout,xout);
    view(10,70);
    grid off;
    colormap winter
    axis([z0 zL 0 1 min(min(xout)) max(max(xout))]);
    xlabel('z');
    ylabel('t');
    zlabel('x(z,t)');

%... read the stopwatch timer
    tcpu=toc;


