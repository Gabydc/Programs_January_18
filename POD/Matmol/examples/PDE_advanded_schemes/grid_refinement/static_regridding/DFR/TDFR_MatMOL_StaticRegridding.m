%... The MatMol Group (2016)
    function TDFR_MatMOL_StaticRegridding(Pe,n0,zswitch,Reltol,Abstol)

%... Function....... : Transient simulation of a jacketed tubular reactor with
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "F. Logist, I. Smets, A. Vande Wouwer, J. Van Impe 2005.
%...                   Optimal control of dispersive tubular chemical reactors:
%...                   Part II" Procs of 16th IFAC World Congress, DVD-ROM, 6p.
%... Technique...... : Solution with the MatMOL toolbox and static regridding
%...                   AGEREG
%... Input...........: Pe:       Peclet number
%...                   n0:       number of vertex centered discritisation points 
%...                   zswitch:  position where the jacket temperature switches 
%...                             from its maximum to its minimum value 
%...                   RelTol:   Relatitve integration tolerance
%...                   AbsTol:   Absolute integration tolerance
%**************************************************************************
    close all

    global D1 D2 Cin Tin beta1 delta1 v k0 E R Tw
    global n Dz Dzz dz1 nsteps maxsteps

    if nargin < 1, Pe       = 1.0e+6;   end
    if nargin < 2, n0       = 51;       end
    if nargin < 3, zswitch  = 0.54;     end
    if nargin < 4, Reltol   = 1.0e-4;   end
    if nargin < 5, Abstol   = 1.0e-6;   end

    tic %... initialise clock


%... 1 Parameters

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

%... Note: 
%... To avoid confusion with parameters from the regridding code, the reactor 
%... related parameters beta and delta are given an additional index 1:
    beta1	=	0.2;                %... [1/s]     dimensionless heat transfer
    delta1	=   0.25;               %... [-]       dimensionless heat of reaction

%... Spatial grid
    n   	=   n0;                 %... [-]       # of discretisation points
    nz      =   n;                  %... [-]       # of discretisation points
    zL      =   0.0;                %... [m]       left border
    zR      =   L;                  %... [m]       right border
    dz      =   (zR-zL)/(nz-1);     %... [m]       spatial discretisation length
    Z       =   [zL:dz:zR]';        %... [m]       position vector

%... Temporal grid
    t_start =   0.0;                %... [s]       initial time
    t_end   =   20.0;               %... [s]       final time
    dt_plot =   1.0;                %... [s]       time interval for plots
    t_plot  =   t_start:dt_plot:t_end+dt_plot;  %... [s] time points to be plotted
     
%... Parameters of the adaptive grid (agereg)
    npdes   =   2;                  %... [-]       nr of PDEs
    nzmax   =   301;                %... [-]       maximum nr of grid points
    maxsteps    = 10;               %... [-]       regridding after ... steps

%... Values for monitor function 1: f = sqrt(alpha+sum(x_z)^2)
    imesh   =   1;                  %... [-]       select monitorfunction
                                    %...           =0: sqrt(alpha+sum(x_z)^2)
                                    %...           =1: sqrt(alpha+max(dabs(x_zz)))
    ilim    =   1;                  %... [-]       select monitor function limiter
                                    %...           =0: no limit
                                    %...           =1: limit related to beta 
    alpha   =   0.01;               %... [-]       agereg parameter to limit
                                    %...           the maximum interval length
    beta    =   300.0;              %... [-]       agereg parameter to limit 
                                    %...           clustering of grid points
    tolz    =   0.04;               %... [-]       equidistribution tolerance
    bound   =   1.2;                %... [-], >1   agereg parameter defining a 
                                    %...           locally bounded grid

%... %... Values for monitor function 2: f = sqrt(alpha+max(dabs(x_zz)))
%... imesh   =   0;
%... ilim    =   1;
%... alpha   =   2;
%... beta    =   1800.0;
%... tolz    =   0.02;
%... bound   =   1.2;


%... 2 Numerical computation

%... 2a Differentiation matrices
    Dz  = five_point_biased_upwind_D1(Z,1);
    Dzz = three_point_centered_D2(Z);

%... 2b Initial conditions
    T0             = Tin*ones(size(Z));
    C0	           = 0.0*ones(size(Z));
    dz1            = Z(2)-Z(1);
    IC(1:2:2*nz-1) = T0;
    IC(2:2:2*nz)   = C0; 
    Tw(Z-z1<=0,1)  = Twmax;
    Tw(Z-z1> 0,1)  = Twmin;
    clear T0 C0

%... 2c Initialise integration
    tk                 = t_start;
    i_plot             = 1;      
    Profiles(i_plot).t = t_start;
    Profiles(i_plot).n = nz;
    Profiles(i_plot).z = Z;
    Profiles(i_plot).T = IC(1:2:2*nz-1);
    Profiles(i_plot).C = IC(2:2:2*nz);
    i_plot             = i_plot+1;
    i_loop             = 0;
    Stats.nsteps       = 0;
    Stats.nfailed      = 0;
    Stats.nfevals      = 0;
    Stats.npds         = 0;
    Stats.ndecomps     = 0;
    Stats.nsolves      = 0;

%... 2d Run the integration loop
    while tk <= 0.9999*t_end
    %... Do the integration for maxsteps steps in a loop until t becomes 
    %... larger than t_end. However, to allow t_end to be part of the solution
    %... the integration span is extended to t_end+dt_plot.
        tspan   = [tk t_plot(i_plot):dt_plot:t_end+dt_plot];

    %... Initialize step counter
        nsteps = 0;

    %... Integration
        options = odeset('RelTol',Reltol,'AbsTol',Abstol);
        options = odeset(options,'JPattern',jpattern(nz));
        options = odeset(options,'Events',@events_StatRegrid);
        sol = ode23s(@ODE_MatMOL_StaticRegridding,tspan,IC,options);

    %... Evaluate solution at elements of tspan lower then integration end
        Y  = (deval(sol,[tspan(tspan<sol.x(end)),sol.x(end)]))';
        X  = Y(end,:);
        X1 = X(1:2:2*nz-1);
        X2 = X(2:2:2*nz);
        tk = sol.x(end);
        clear IC
    
    %... Integrator statisctics
        Stats = AddIntStats(Stats,sol.stats,sol.solver); %... solution statistics
        clear sol sol2
    
    %... Save intermediate results
        i1 = 2;
        while tk >= t_plot(i_plot) 
            Profiles(i_plot).t = tspan(i1);
            Profiles(i_plot).n = nz;
            Profiles(i_plot).z = Z;
            Profiles(i_plot).T = Y(i1,1:2:2*nz-1);
            Profiles(i_plot).C = Y(i1,2:2:2*nz);
            i_plot  = i_plot+1;
            i1      = i1 + 1;           
        end
    
    %... Compute dimensionless independent variables
        Xdimless = zeros(1,2*nz);
        Xdimless(1:2:end-1) = (X(end,1:2:end-1)-Tin)/Tin;
        Xdimless(2:2:end  ) = (Cin-X(end,2:2:end  ))/Cin;
        
    %... Refine the grid
        [Z_new,nz_new,ier,tolz_new]=agereg(Z,Xdimless,npdes,nzmax,alpha,...
                                           beta,tolz,bound,imesh,ilim);

    %... Interpolate the dependent variables
        X1_new = spline(Z,X1,Z_new);
        X2_new = spline(Z,X2,Z_new);
        clear Z X X1 X2

    %... Update the grid, its variables and the dependent variables
        Z    = Z_new';
        nz   = nz_new;
        n    = nz;
        tolz = tolz_new;
    
        Tw            = [];
        Tw(Z-z1<=0,1) = Twmax;
        Tw(Z-z1> 0,1) = Twmin;
    
        IC(1:2:2*nz-1) = X1_new;
        IC(2:2:2*nz  ) = X2_new;
        clear Z_new X1_new X2_new;
     
    %... Compute new differentiation matrices
        Dz  = []; Dzz = [];
        Dz  = five_point_biased_upwind_D1(Z,1);      
        Dzz = three_point_centered_D2(Z);
    
    %... Update loop counter
        i_loop = i_loop+1;
    end
    CPUtime = toc; %... read clock

    disp('Integration statistics:')
    fprintf('  %g solver restarts\n',  i_loop);
    fprintf('  %g successful steps\n', Stats.nsteps);
    fprintf('  %g failed attempts\n',  Stats.nfailed);
    fprintf('  %g function evaluations\n', Stats.nfevals);
    fprintf('  %g partial derivatives\n',  Stats.npds);
    fprintf('  %g LU decompositions\n', Stats.ndecomps);
    fprintf('  %g solutions of linear systems\n', Stats.nsolves);
    fprintf('  %g s (total) computation time\n',  CPUtime);

%... 3 Plot results

    for i1=1:i_plot-1
        figure(1)
        subplot('position',[0.1 0.3 0.8 0.6])
        plot(Profiles(i1).z,Profiles(i1).T,'g',Profiles(i1).z,Profiles(i1).T,'.g')
        ylabel('Temperature T [K]','Fontsize',18)
        title(['MatMOL Static Regridding: Pe1=Pe2=' num2str(Pe2)],'Fontsize',18)
        axis([zL zR 270.0 450.0])
        set(gca,'Fontsize',16)
        hold on
    
        subplot('position',[0.1 0.1 0.8 0.15])
        plot(Profiles(i1).z,Profiles(i1).t*ones(length(Profiles(i1).z),1),'.g')
        xlabel('Space z [m]','fontsize',18);
        ylabel('Time [s]','Fontsize',18)
        axis([zL zR t_start t_end])       
        set(gca,'Fontsize',16)
        hold on
    
        figure(2)
        subplot('position',[0.1 0.3 0.8 0.6])
        plot(Profiles(i1).z,Profiles(i1).C,'g',Profiles(i1).z,Profiles(i1).C,'.g')
        ylabel('Concentration C [mole/K]','Fontsize',18)
        title(['MatMOL Static Regridding: Pe1=Pe2=' num2str(Pe2)],'Fontsize',18)
        axis([zL zR -0.001 0.02])
        set(gca,'Fontsize',16)
        hold on
    
        subplot('position',[0.1 0.1 0.8 0.15])
        plot(Profiles(i1).z,Profiles(i1).t*ones(length(Profiles(i1).z),1),'.g')
        xlabel('Space z [m]','fontsize',18);
        ylabel('Time [s]','Fontsize',18)
        axis([zL zR t_start t_end])    
        set(gca,'Fontsize',16)
        hold on
    
        figure(3)
        plot(Profiles(i1).t,Profiles(i1).n,'g*');
        xlabel('Time [s]','Fontsize',18)
        ylabel('Nr of grid points [-]','Fontsize',18)    
        set(gca,'Fontsize',18)
        axis([t_start t_end 0 nzmax])
        hold on 
    end

    end %... End of function TDFR_MatMOL_StaticRegridding


    function dX = ODE_MatMOL_StaticRegridding(t,X)
%... Provides the right hand side of the ODEs

    global D1 D2 Cin Tin beta1 delta1 v k0 E R Tw
    global n Dz Dzz dz1

%... Independant variables
    T(1:n,1) = X(1:2:2*n-1);
    C(1:n,1) = X(2:2:2*n);

%... Boundary conditions at z = 0
    T(1,1) = (T(2,1)+v/D1*dz1*Tin)/(1+v/D1*dz1);
    C(1,1) = (C(2,1)+v/D2*dz1*Cin)/(1+v/D2*dz1);

%... Spatial derivatives
    Tz  = Dz*T;
    Cz  = Dz*C;
    Tzz = Dzz*T;
    Czz = Dzz*C;

%... Temporal derivatives
    Tt = D1*Tzz - v*Tz + Tin*delta1*k0*C.*exp(-E./(R*T))/Cin + beta1*(Tw-T);
    Ct = D2*Czz - v*Cz - k0*C.*exp(-E./(R*T));

    dX(1:2:2*n-1,1) = Tt;
    dX(2:2:2*n,1)   = Ct;   

    end %... End of function ODE_MatMOL_StaticRegridding

    
    function S = jpattern(n)
%... Sparsity pattern of the Jacobian matrix

    bloc0 = zeros(2,2);
    bloc1 = ones(2,2);
    bloc2 = eye(2);
    l1 = [bloc1 repmat(bloc2,1,4) repmat(bloc0,1,n-5)];
    l2 = [bloc2 bloc1 repmat(bloc2,1,3) repmat(bloc0,1,n-5)];
    l3 = [repmat(bloc2,1,2) bloc1 repmat(bloc2,1,2) repmat(bloc0,1,n-5)];
    for i=1:n-5
        l(1+(i-1)*2:i*2,1:n*2)=[repmat(bloc0,1,i-1) repmat(bloc2,1,3) ...
                                bloc1 repmat(bloc2,1,2) repmat(bloc0,1,n-i-5)];
    end
    lm2 = [repmat(bloc0,1,n-5) repmat(bloc2,1,3) bloc1 bloc2];
    lm1 = [repmat(bloc0,1,n-5) repmat(bloc2,1,4) bloc1];
    X = cat(1,l1,l2,l3,l,lm2,lm1);
    S = sparse(X);

    end %... End of function jpattern
