%... The MatMol Group (2016)
%... Algae Bloom (solution using spectral methods with 4 points)
    close all
    clear all

%... start a stopwatch timer
    tic

%... set global variables
    global mu r K;
    global z0 zL z n D1;
    global tresh fac vectorized

%... model parameters
    mu = 0.01;
    r = 5;
    K = 1;

%... spatial grid
    z0 = 0.0;

%... select initial diameter of the algae patch (when zL >
%... pi*sqrt(mu/r), the patch will grow, otherwise it will decay)
    diameter = 'large';
    % diameter = 'small';
    switch diameter
        %... the patch is large at the initial time
        case('large')
            zL = 1.5*pi*sqrt(mu/r);
            %... the patch is too small
        case('small')
            zL = 0.9*pi*sqrt(mu/r);
    end
    n = 4;   %... Number of spatial grid points

%... spectral differentiation matrix on a nonuniform grid (stagewise
%... differentiation)
    [D1,z] = chebyshev_spectral_D1(z0,zL,n);

%... initial conditions
    x0 = 0.1*sin(pi*z/zL);

%... call to ODE solver (comment/decomment one of the methods to select a solver)
    t0         = 0;
    tf         = 3;
    Dtplot     = 0.3;
    fac        = [];
    tresh      = 1e-12;
    vectorized = 1;

    options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    options = odeset(options,'Vectorized','on');
    options = odeset(options,'Jacobian',@jacobian_algae);
    figure(1);
    colordef white;
    options      = odeset(options,'OutputFcn','odeplot');
    [tout, xout] = ode15s(@algae_bloom_pde,t0:Dtplot:tf,x0,options);

%... plot results
    set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
           'defaultlinelinewidth',1.5,'defaultpatchlinewidth',.7);
    colordef white;
    figure(2)
    zint = linspace(z0,zL,100);
    for i = 1:length(tout),
        xint(i,:) = polyval(polyfit(z,xout(i,:)',n),zint);
        plot(zint,xint(i,:),'-g');
        hold on
    end
    xlabel('z');
    ylabel('u');
    title('Algae bloom')
%...
    figure(4)
%...     mesh(z,tout,yout);
    surf(z,tout,xout);
    grid off;
%...     colormap([0 0 0]);
    axis([z0 zL t0 tf min(min(xout)) max(max(xout))]);
    xlabel('z');
    ylabel('t');
    zlabel('\rho(z,t)');
%...
%... read the stopwatch timer
    tcpu=toc;
