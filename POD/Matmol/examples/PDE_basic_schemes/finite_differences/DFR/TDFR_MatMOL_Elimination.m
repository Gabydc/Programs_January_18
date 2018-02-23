%... The MatMol Group (2016)
    function TDFR_MatMOL_Elimination(Pe,n0,zswitch,orderD1,orderD2,Reltol,Abstol)
%... Function....... : Transient simulation of a jacketed tubular reactor with
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "F. Logist, I. Smets, A. Vande Wouwer, J. Van Impe 2005.
%...                   Optimal control of dispersive tubular chemical reactors:
%...                   Part II" Procs of 16th IFAC World Congress, DVD-ROM, 6p.
%... Technique...... : Solution with the MatMOL toolbox using explicit ODEs 
%...                   solvers (elimination)
%... Input...........: Pe:       Peclet number
%...                   n0:       number of vertex centered discritisation points 
%...                   zswitch:  position where the jacket temperature switches 
%...                             from its maximum to its minimum value 
%...                   orderD1:  approximation order of the first order 
%...                             spatial derivatives
%...                   orderD2:  approximation order of the second order 
%...                             spatial derivatives
%...                   RelTol:   Relatitve integration tolerance
%...                   AbsTol:   Absolute integration tolerance

    close all

    global D1 D2 Cin Tin beta delta v k0 E R Tw
    global n Dz Dzz

    if nargin < 1, Pe       = 100.0;    end  
    if nargin < 2, n0       = 101;      end  
    if nargin < 3, zswitch  = 0.54;     end  
    if nargin < 4, orderD1  = 1;        end  
    if nargin < 5, orderD2  = 2;        end  
    if nargin < 6, Reltol   = 1.0e-4;   end  
    if nargin < 7, Abstol   = 1.0e-6;   end  

    tic %... initialise clock


%... 1 Parameters
%... ************

%... Process
    v		=  0.1;                 %... [m/s]     fluid velocity
    L		=  1.0;                 %... [m]       reactor length
    Pe1	    =  Pe;					%... [-]       Peclet number
    Pe2	    =  Pe;					%... [-]       Peclet number
    D1		=  v*L/Pe1;				%... [m^2/s]   dispersion coefficient
    D2      =  v*L/Pe2;				%... [m^2/s]   dispersion coefficient

    Tin     =  340.0;               %... [K]       feed temperature
    Twmin   =  280.0;               %... [K]       minimum jacket temperature
    Twmax   =  400.0;               %... [K]       maximum jacket temperature
    Cin     =  0.02;                %... [mol/L]   feed concentration

    z1      =  zswitch;             %... [m]       switching position
    E		=  11250.0;             %... [cal/mol] activation energy
    R		=  1.986;               %... [cal/(mol*K)] gass constant
    k0		=  1.0e+6;              %... [1/s]     rate constant

    beta	=	0.2;                %... [1/s]     dimensionless heat transfer
    delta	=   0.25;               %... [-]       dimensionless heat of reaction

%... Spatial grid
    n   	=   n0;                 %... [-]       # of discretisation points
    zL      =   0.0;                %... [m]       left border
    zR      =   L;                  %... [m]       right border
    dz      =   (zR-zL)/(n-1);      %... [m]       spatial discretisation length
    Z       =   zL:dz:zR;           %... [m]       position vector

%... Temporal grid
    t_start =   0.0;                %... [s]       initial time
    t_end   =   20.0;               %... [s]       final time
    dt_plot =   1.0;                %... [s]       time interval for plots


%... 2 Numerical computation
%... ***********************

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
    T0  = Tin*ones(length(Z)-2,1);
    C0	=   0*ones(length(Z)-2,1);
    IC  = [T0;C0];
    Tw  = [Twmax*ones(1+round(z1*(n-1)),1); ...
           Twmin*ones(n-(1+round(z1*(n-1))),1)];
    clear T0 C0
    
%... 2c Solution of the DAEs
    options = odeset('RelTol',Reltol,'AbsTol',Abstol);
    sol = ode45(@ODE_MatMOL_Explicit2,[t_start t_end],IC',options);

%... 2d Reactor properties
    t    = t_start:dt_plot:t_end;   %... time points to be plotted
    sol2 = deval(sol,t);            %... evaluate solution at time points
    T1   = sol2(1:n-2,:);
    C1   = sol2(n-1:2*n-4,:);
    T    = [(D1*Dz(1,2:n-1)*T1+v*Tin)./(v-D1*Dz(1,1)); T1; ...
            -(Dz(end,2:n-1)*T1       ./Dz(end,end))];
    C    = [(D2*Dz(1,2:n-1)*C1+v*Cin)./(v-D2*Dz(1,1)); C1; ...
            -(Dz(end,2:n-1)*C1       ./Dz(end,end))];
    clear T1 C1

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
    mesh(t,Z,T);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [s]','Fontsize',18)
    zlabel('Temperature T [K]','Fontsize',18)
    title(['Temperature evolution (MatMOL Explicit2: z1=' num2str(z1) ...
        ' m, Pe1=' num2str(Pe1) ', Pe2=' num2str(Pe2) ', n=' num2str(n) ')'],...
        'Fontsize',20)
    axis([t_start t_end zL zR 270.0 450.0])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')

    figure(2)
    mesh(t,Z,C);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [s]','Fontsize',18)
    zlabel('Concentration C [mole/L]','Fontsize',18)
    title(['Concentration evolution (MatMOL Explicit2: z1=' num2str(z1) ...
            ' m, Pe1=' num2str(Pe1) ', Pe2=' num2str(Pe2) ', n=' num2str(n) ')'],...
            'Fontsize',20)
    axis([t_start t_end zL zR 0.0 0.02])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')
    
    figure(3)
    plot(Z,T);
    xlabel('Space z [-]','fontsize',18);
    ylabel('Temperature T [K]','Fontsize',18)
    title(['MatMOL Explicit2: z1=' num2str(z1) ' m, Pe1=Pe2=' num2str(Pe2) ...
        ', n=' num2str(n) ', RT=' num2str(Reltol) ', AT=' num2str(Abstol)],...
        'Fontsize',18)
    axis([zL zR 270.0 450.0])
    grid on
    set(gca,'Fontsize',18)
   
    figure(4)
    plot(Z,C);
    xlabel('Space z [-]','fontsize',18);
    ylabel('Concentration C [mole/K]','Fontsize',18)
    title(['MatMOL Explicit2: z1=' num2str(z1) ' m, Pe1=Pe2=' num2str(Pe2) ...
        ', n=' num2str(n) ', RT=' num2str(Reltol) ', AT=' num2str(Abstol)],...
        'Fontsize',18)
    axis([zL zR -0.001 0.02])
    grid on
    set(gca,'Fontsize',18)

    end %... End of function TDFR_MatMOL_Elimination
%... *************************************************************************



    function dX = ODE_MatMOL_Explicit2(t,X)
%... Provide the right hand side of the ODEs with eliminated boundary 
%... conditions

    global D1 D2 Cin Tin beta delta v k0 E R Tw
    global n Dz Dzz

%... Independent variables
    T1 = X(1:n-2);
    C1 = X(n-1:2*n-4);

%... Algebraic introduction of Danckwerts BC
    T = [(D1*Dz(1,2:n-1)*T1+v*Tin)./(v-D1*Dz(1,1)); T1; ...
          -Dz(end,2:n-1)*T1       ./Dz(end,end)];
    C = [(D2*Dz(1,2:n-1)*C1+v*Cin)./(v-D2*Dz(1,1)); C1; ...
          -Dz(end,2:n-1)*C1       ./Dz(end,end)];

%... Spatial derivatives
    Tz = Dz*T;
    Cz = Dz*C;
    Tzz = Dzz*T;
    Czz = Dzz*C;

%... Temporal derivatives
    Tt = D1*Tzz - v*Tz + Tin/Cin*delta*k0*C.*exp(-E./(R*T)) + beta.*(Tw-T);
    Ct = D2*Czz - v*Cz - k0*C.*exp(-E./(R*T));

%... RHS of ODE
    dX	= [Tt(2:n-1);Ct(2:n-1)];

    end %... End of function ODE_MatMOL_Explicit2
%... *************************************************************************