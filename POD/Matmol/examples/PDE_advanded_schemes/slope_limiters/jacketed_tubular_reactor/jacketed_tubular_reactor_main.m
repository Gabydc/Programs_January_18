%... The MatMol Group (2016)
%...
%... Jacketed tubular reactor
%...
%...Transient simulation of a jacketed tubular reactor with dispersion 
%...and Danckwerts boundary conditions
%...
%...   x1t = Diff1*x1zz - v*x1z + alpha*delta*(1-x2)*exp(gamma*x1/(1+x1))
%...
%...         + beta*(u-x1)
%...
%...
%...   x2t = Diff2*x2zz - v*x2z + alpha*(1-x2)*exp(gamma*x1/(1+x1))
%...
%...
%...   C.I. : x1(z,0) = 0
%...
%...          x2(z,0) = 1
%...
%...
%...   C.L. : en z = 0 : D1*x1z = v*x1
%...
%...                     D2*x2z = v*x2
%...
%...          en z = L : x1z = 0
%...
%...                     x2z = 0
%...
%...   x1 and x2 are normalized variables : x1 = (T-Tin)/Tin      
%...                                        x2 = (cin-c)/cin
%...
%...   Numerical values :
%...
%...   Diff1 = Diff2 = v*L/Pe     with   Pe = 0.01, 100, 1e6
%...   v = 0.1
%...   L = 1
%...   alpha = k0*exp(-E/(R*Tin))
%...   k0 = 1e6
%...   E = 11250
%...   R = 1.986
%...   Tin = 340
%...   cin = 0.02
%...   delta = 0.25
%...   gamma = E/(R*Tin)
%...   beta = 0.2
%...   u = (Tw-Tin)/Tin    with   Tin = 340    
%...                       and     Tw = 400   pour  z <= zl
%...                               Tw = 280   pour  zl < z
%...
%...               Pe   |   zl 
%...              -------------
%...               0.01 |  0.58
%...               100  |  0.54
%...               1e6  |  0.54
%...
%...         "F. Logist, I. Smets, A. Vande Wouwer, J. Van Impe 2005.
%...          Optimal control of dispersive tubular chemical reactors:
%...          Part II" Procs of 16th IFAC World Congress, DVD-ROM, 6p.
%...
%... the following code computes a solution to Burgers' equation
%...
    close all
    clear all
%...
%... start a stopwatch timer
    tic
%...
%... set global variables
    global Diff1 Diff2 v alpha delta gamma beta u uL uR
    global D1 D2 n dz nl z ne method
%...
%... model parameters
    v = 0.1;
    L = 1;
    k0 = 1e6;
    E = 11250;
    R = 1.986;
    Tin = 340;
    delta = 0.25;
    beta = 0.2;
    Pe = 1e6;%100;%0.01;%
    zl = 0.54;%0.58;%
    Diff1 = v*L/Pe;     
    Diff2 = v*L/Pe;     
    alpha = k0*exp(-E/(R*Tin));
    gamma = E/(R*Tin);
    uL = (400-Tin)/Tin;
    uR = (280-Tin)/Tin;
    cin = 0.02;
    ne = 2;
%...
%... spatial grid
    n = 601;
    nl = fix(n*zl/L);
    dz = L/(n-1);
    z = [0:dz:L]';
    u = zeros(n,1);
%...
%... select slope limiter approximation
%      method = 'Van_Leer'         %... Van Leer slope limiter
%      method = 'superbee'         %... Superbee slope limiter
%      method = 'smart'            %... Smart slope limiter
%      method = 'mc'               %... monotonized center limiter
%      method = 'minmod'           %... Minmod slope limiter
%      method = 'Koren'            %... Koren slope limiter
       method = 'KT'               %... Kurganov Tadmor centered scheme
%...
%...
%... differentiation matrix
    D1 = five_point_biased_upwind_D1(z,1);%three_point_centered_D1(z);%two_point_upwind_D1(z,1),%
    D2 = three_point_centered_D2(z);
%...
%... initial conditions
    x1 = zeros(n,1);
    x2 = ones(n,1);
    x = [x1;x2];
%...  
%... call to ODE solver 
%...
    time = 0:1:50;
%...
%... flindex = 1 if upwind flux limiter and zero otherwise
%...
    flindex = 0;
    option = odeset('JPattern',JPattern(n,ne,flindex));
    [timeout,xout] = ode45(@jacketed_tubular_reactor_pdes,time,x,option);
%...
%...  impose BCs
%...
    xout(:,1) = xout(:,2)/(1+dz*v/Diff1);
    xout(:,n+1) = xout(:,n+2)/(1+dz*v/Diff2);
    xout(:,n) = 2*xout(:,n-1)-xout(:,n-2);
    xout(:,2*n) = 2*xout(:,2*n-1)-xout(:,2*n-2);
%...
%...  plot results
%...
    figure
    hold
    plot(z,Tin*(1+xout(:,1:n)))
    xlabel('z');
    ylabel('T(z,t)');
%...
    figure
    hold
    plot(z,cin*(1-xout(:,n+1:2*n)))
    xlabel('z');
    ylabel('c(z,t)');
%...
%... read the stopwatch timer
    toc