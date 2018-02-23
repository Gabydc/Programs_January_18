%... The MatMol Group (2016)
%... Distributed parameter settler Model
%... R. David, P. Saucez, J.-L. Vasel, A. Vande Wouwer, 
%... Modeling and Numerical Simulation of Secondary Settlers: a Method of Lines Strategy
%... Water Research 43 (2009), 319-330.

    close all
    clear all

%... Global variables
    global nu0 rp rh nuprim0 Cmin
    global A
    global Cf Qf qun qov
    global n n1 n2
    global vec_D
    global D1_1  D1_2  D2_1  D2_2 D11 D12 D1_1_

%... Eefine the model parameters
    L   = 4;         %... m
    A   = 500;       %... m^2
    In  = 1.8;       %... m
    Qe  = 250;       %... m^3/h
    Qun = 200;       %... m^3/h
    Qf  = Qe + Qun;  %... m^3/h
    qun = Qun/A;     %... m/h (also denoted q_S)
    qov = Qe/A;      %... m/h (also denoted q_W)

%... JEPPSSON, U. et DIEHL, S. (1996). An evaluation of a dynamic model of the secondary
%... clarifier. Water Science and Technology, 34(5-6):19-26.
%... nu0     = 145/24;    %... m/h
%... rp      = 5e-3;      %... m^3/g
%... rh      = 0.42e-3;   %... m^3/g
%... nuprim0 = 100/24;    %... m/h
%... Cf      = 5.2*1000;  %... g/m^3
%... Cmin    = 10;        %... g/m^3
%
%... COST Benchmark
    nu0     = 474/24;
    nuprim0 = 250/24;
    rh      = 5.76e-4;
    rp      = 2.86e-3;
    fns     = 0.00228;
    Cf      = 5.2*1000;  %... g/m^3
    Cmin    = fns*Cf;    %... g/m^3

    D       = 13/24;     %... m^2/h
    D1      = D;
    D2      = D;

%... Spatial grid
    n     = 210;
    n1    = round((n)*In/L)+1;
    n2    = round((n)*(L-In)/L)+1;
    vec_D = [D1*ones(n1,1); D2*ones(n2,1)]; %... vector of diffusion on the spatial grid
    z01   = 0;
    zL1   = In;
    z02   = zL1;
    zL2   = L;
    dz1   = (zL1-z01)/(n1-1);
    z1    = [z01 : dz1 : zL1]';
    dz2   = (zL2-z02)/(n2-1);
    z2    = [z02 : dz2 : zL2]';
    z     = [z1 ; z2];

%... Differentiation matrices in zone 1 (for convective and diffusive terms)
    D1_1_ = five_point_biased_upwind_D1(z1,-1);
    D1_1  = five_point_biased_upwind_D1(z1,1);
    D11   = five_point_centered_D1(z1);
    D2_1  = five_point_centered_D2(z1);

%... Differentiation matrices in zone 2 (for convective and diffusive terms)
    D1_2 = five_point_biased_upwind_D1(z2,1);
    D12  = five_point_centered_D1(z2);
    D2_2 = five_point_centered_D2(z2);

%... Initial conditions
    X0 = 0.0;
    x  = X0*ones( n1+n2, 1);

%... Plot times
    t = [0, 0.1, 0.2, 0.4,1, 2, 3, 10, 1000];

%... ODEs solver parameters
    epsilon = 0;
    options = odeset('Mass',diag([ epsilon ones(1,n1-2) epsilon epsilon ones(1,n2-2) epsilon]),...
                     'RelTol', 1e-3, 'AbsTol', 1e-3);

%... Call to ODE solver             
    [tout, yout] = ode15s(@settler_pdes, t, x, options);

%... Plot the results
    y_max = 11000;
    y_min = 0;

    figure(1)
    label01 = {[''] [''] [''] [''] [''] [''] ['\bf z  [m]'] ['\bf z  [m] '] ['\bf z  [m]']};
    label02 = {['\bf C  [g/m^3]'] [''] [''] ['\bf C  [g/m^3]'] [''] [''] ['\bf C  [g/m^3]'] [' '] ['']};
    for i = 1:9
        subplot(3,3,i)
        plot(z1, yout(i , 1:n1), '-b', 'MarkerSize', 2, 'LineWidth',1);
        hold on
        plot(z2, yout(i, n1+1:n1+n2), '-b', 'MarkerSize', 2, 'LineWidth',1);
        text(0.05, 9000, ['\bf t = ',num2str(t(i)),'\bf h'], 'FontSize', 16);
        set(gca, 'XTick',[0 In L],'XLim',[0 L], 'YLim',[y_min y_max], 'FontSize', 16);
        grid
        xlabel(label01{i}, 'FontSize', 16);
        ylabel(label02{i}, 'FontSize', 16);
    end
