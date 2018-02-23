%... The MatMol Group (2016)
    function TFBBR_MatMOL_Elimination(Pe,n0,orderD1,orderD2,Reltol,Abstol)

%... Function....... : Transient simulation of a tubular biochemical reactor
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "M. Laabissi, J.J. Winkin, D. Dochain and M.E. Achhab 2005. 
%...                   Dynamical analysis of a tubular biochemical reactor 
%...                   infinite-dimensional nonlinear model" Procs of the 44th 
%...                   IEEE Conference on Decision and Control & the European 
%...                   Control Conference, 5965-5970.
%... Technique...... : Solution with the MatMOL toolbox using explicit ODEs 
%...                   solvers (elimination)
%... Input...........: Pe:       Peclet number
%...                   n0:       number of vertex centered discritisation points 
%...                   orderD1:  approximation order of the first order
%...                             spatial derivatives
%...                   orderD2:  approximation order of the second order
%...                             spatial derivatives
%...                   RelTol:   Relative integration tolerance
%...                   AbsTol:   Absolute integration tolerance
%**************************************************************************
    close all

    global D1 v Sin mu0 k0 kd KS KI
    global n Dz Dzz

    if nargin < 1, Pe       = 100.0;    end  
    if nargin < 2, n0       = 101;      end  
    if nargin < 3, orderD1  = 1;        end  
    if nargin < 4, orderD2  = 2;        end  
    if nargin < 5, Reltol   = 1.0e-4;   end  
    if nargin < 6, Abstol   = 1.0e-6;   end  

    tic %... initialise clock

%... 1 Parameters
%... ************

%... Process
    v		=  1.0;                 %... [m/h]     fluid velocity
    L		=  1.0;                 %... [m]       reactor length
    Pe1	    =  Pe;					%... [-]       Peclet number
    D1		=  v*L/Pe1;				%... [m^2/h]   dispersion coefficient

    Sin     =  20.0;                %... [g/L]     feed concentration

    mu0     =  0.4;                 %... [1/h]     maximum growth rate
    k0      =  2.0;                 %... [g/L]     yield factor
    kd      =  0.01;                %... [1/h]     death rate
    KS      =  1.0;                 %... [g/L]     substrate inhibition constant
    KI      =  1.0;                 %... [g/L]     biomass inhibition constant

%... Spatial grid
    n   	=   n0;                 %... [m]       # of discretisation points
    zL      =   0.0;                %... [m]       left border
    zR      =   L;                  %... [m]       right border
    dz      =   (zR-zL)/(n-1);      %... [m]       spatial discretisation length
    Z       =   zL:dz:zR;           %... [m]       position vector

%... Temporal grid
    t_start =   0.0;                %... [h]       initial time
    t_end   =   2.0;                %... [h]       final time
    dt_plot =   0.1;                %... [h]       time interval for plots


%... 2 Numerical computation

%... 2a Differentiation matrices
    switch orderD1
        case 1
            Dz 	= two_point_upwind_D1(Z,v);
        case 4
            Dz  = five_point_biased_upwind_D1(Z,v);
        otherwise
            disp('Invalid selection for order of first spatial derivative.')
            disp('Choose 1 for low-order or 4 for high-order scheme.')
    end

    switch orderD2
        case 2
            Dzz = three_point_centered_D2(Z);
        case 4
            Dzz = five_point_centered_D2(Z);
        otherwise
            disp('Invalid selection for order of second spatial derivative.')
            disp('Choose 2 for low-order or 4 for high-order scheme.')
    end

%... 2b Initial conditions
    S0 =   0.0*ones(length(Z)-2,1);
    X0 = 300.0*ones(length(Z)-0,1);
    IC = [S0; X0];
    clear S0 X0
    
%... 2c Solution of the DAEs
    options = odeset('RelTol',Reltol,'AbsTol',Abstol);
    sol = ode45(@ODE_MatMOL_Explicit2,[t_start t_end],IC',options);


%... 2d Reactor properties
    t    = t_start:dt_plot:t_end;   %... time points to be plotted
    sol2 = deval(sol,t);            %... evaluate solution at time points
    S1   = sol2(1:n-2,:);
    S    = [(D1*Dz(1,2:n-1)*S1+v*Sin)./(v-D1*Dz(1,1)); S1; ...
            -(Dz(end,2:n-1)*S1       ./Dz(end,end))];
    X    = sol2(n-1:2*n-2,:);
    clear S1


%... 2e Integrator statisctics
    Stats   = sol.stats; %... solution statistics
    CPUtime = toc;       %... read clock
    clear sol sol2

    disp('Integration statistics:')
    fprintf('  %g successful steps\n', Stats.nsteps);
    fprintf('  %g failed attempts\n',  Stats.nfailed);
    fprintf('  %g function evaluations\n', Stats.nfevals);
    fprintf('  %g s (total) computation time\n',  CPUtime);

%... 3 Plot results
%... **************

    figure(1)
    mesh(t,Z,S);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [h]','Fontsize',18)
    zlabel('Substrate S [g/L]','Fontsize',18)
    title(['Substrate evolution (MatMOL: Pe=' ...
            num2str(Pe1) ', n=' num2str(n) ')'],'Fontsize',20)
    axis([t_start t_end zL zR -1.0 21.0])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')

    figure(2)
    mesh(t,Z,X);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [h]','Fontsize',18)
    zlabel('Biomass X [g/L]','Fontsize',18)
    title(['Biomass evolution (MatMOL: Pe=' ...
            num2str(Pe1) ', n=' num2str(n) ')'],'Fontsize',20)
    axis([t_start t_end zL zR 295.0 305.0])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')
    
    figure(3)
    plot(Z,S);
    xlabel('Space z [m]','fontsize',18);
    ylabel('Substrate S [g/L]','Fontsize',18)
    title(['MatMOL: Pe=' num2str(Pe1) ', n=' ...
            num2str(n) ', RT=' num2str(Reltol) ', AT=' num2str(Abstol)],'Fontsize',18)
    axis([zL zR -1.0 21.0])
    grid on
    set(gca,'Fontsize',18)
   
    figure(4)
    plot(Z,X);
    xlabel('Space z [m]','fontsize',18);
    ylabel('Biomass X [g/L]','Fontsize',18)
    title(['MatMOL: Pe=' num2str(Pe1) ', n=' ...
            num2str(n) ', RT=' num2str(Reltol) ', AT=' num2str(Abstol)],'Fontsize',18)
    axis([zL zR 295 305.0])
    grid on
    set(gca,'Fontsize',18)

    end %... End of function TFBBR_MatMOL_Elimination




    function dX = ODE_MatMOL_Explicit2(t,Y)
%... Provide the right hand side of the ODEs with eliminated boundary 
%... conditions

    global D1 v Sin mu0 k0 kd KS KI
    global n Dz Dzz

%... Independent variables
    S1 = Y(1:n-2);
    X  = Y(n-1:2*n-2);

%... Algebraic introduction of Danckwerts BC
    S = [(D1*Dz(1,2:n-1)*S1+v*Sin)./(v-D1*Dz(1,1)); S1; ...
          -Dz(end,2:n-1)*S1       ./Dz(end,end)];

%... Reaction rates
    mu = mu0*S./(KS*X + S + 1/KI*S.^2);

%... Spatial derivatives
    Sz = Dz*S;
    Szz = Dzz*S;

%... Temporal derivatives
    St = D1*Szz - v*Sz - k0*mu.*X;
    Xt =               - kd*X + mu.*X;

%... RHS of ODEs
    dX = [St(2:n-1); Xt(1:n)];

    end %... End of function ODE_MatMOL_Explicit2
