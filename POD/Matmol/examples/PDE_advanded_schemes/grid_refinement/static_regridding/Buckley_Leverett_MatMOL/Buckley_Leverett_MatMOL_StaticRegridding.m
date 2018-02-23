%... The MatMol Group (2016)
     function Buckley_Leverett_MatMOL_StaticRegridding(eps0,n0,Reltol,Abstol)
%...
%... Function....... : Transient simulation of the classic 1D Buckley-Leverett
%...                   equation for oil reservoirs
%... Technique...... : Solution with the MatMOL toolbox and static regridding
%...                   AGEREG
%... Input...........: eps0:     cappilary pressure term
%...                   n0:       number of vertex centered discritisation points
%...                   RelTol:   Relatitve integration tolerance
%...                   AbsTol:   Absolute integration tolerance

    close all

    global eps
    global n Dz Dzc dz1 nsteps maxsteps

    if nargin < 1, eps0     = 0.001;    end  
    if nargin < 2, n0       = 101;      end  
    if nargin < 3, Reltol   = 1.0e-6;   end  
    if nargin < 4, Abstol   = 1.0e-6;   end  

    tic    % initialise clock


%... 1 Parameters

%... Process
    eps = eps0;               %... [-]       capillary pressure term
    clear eps0

%... Spatial grid
    n  =   n0;                 %... [-]       # of discretisation points
    nz =   n;                  %... [-]       # of discretisation points
    zL =   0.0;                %... [-]       left border position
    zR =   1.0;                %... [-]       right border position
    dz =   (zR-zL)/(n-1);      %... [-]       spatial discretisation length 
    Z  =   [zL:dz:zR]';        %... [-]       dimensionless position vector

%... Temporal grid
    t_start =   0.0;                %... [-]       initial time
    t_end   =   0.4;                %... [-]       final time
    dt_plot =   0.01;               %... [-]       time interval for plots
    t_plot  =   t_start:dt_plot:t_end+dt_plot;  %... [-] time points to be plotted

%... Parameters of the adaptive grid (agereg)
    npdes    = 1;                %... [-]       nr of PDEs
    nzmax    = 301;              %... [-]       maximum nr of grid points
    maxsteps = 40;               %... [-]       regridding after ... steps

%... Values for monitor function 1: f = sqrt(alpha+sum(x_z)^2)
    imesh = 1;                  %... [-]       select monitorfunction
                                %...           =0: sqrt(alpha+sum(x_z)^2)
                                %...           =1: sqrt(alpha+max(dabs(x_zz)))
    ilim  = 1;                  %... [-]       select monitor function limiter
                                %...           =0: no limit
                                %...           =1: limit related to beta 
    alpha = 0.04;               %... [-]       agereg parameter to limit
                                %...           the maximum interval length
    beta  = 400;                %... [-]       agereg parameter to limit 
                                %...           clustering of grid points
    tolz  = 0.04;               %... [-]       equidistribution tolerance
    bound = 2.0;                %... [-], >1   agereg parameter defining a 
                                %...           locally bounded grid

%... Values for monitor function 2: f = sqrt(alpha+max(dabs(x_zz)))
%... imesh   =   0;
%... ilim    =   1;
%... alpha   =   2;
%... beta    =   1800;
%... tolz    =   0.02;
%... bound   =   1.2;


%... 2 Numerical computation

%... 2a Differentiation matrices
    Dzc = four_point_biased_upwind_D1(Z,1);
    Dz  = three_point_centered_D1(Z);

%... 2b Initial conditions
    S0  = max(zeros(n,1),interp1([0 1/3 1],[1 0 0],Z));
    dz1 = Z(2)-Z(1);
    IC  = S0;
    clear S0

%... 2c Initialise integration
    tk      = t_start;
    i_plot  = 1;      
    Profiles(i_plot).t = t_start;
    Profiles(i_plot).n = nz;
    Profiles(i_plot).z = Z;
    Profiles(i_plot).S = IC(1:npdes:npdes*nz);
    i_plot  = i_plot+1;

    i_loop  = 0;
    Stats.nsteps   = 0;
    Stats.nfailed  = 0;
    Stats.nfevals  = 0;
    Stats.npds     = 0;
    Stats.ndecomps = 0;
    Stats.nsolves  = 0;

%... 2d Run the integration loop
    while tk <= 0.9999*t_end
    %... Do the integration for maxsteps steps in a loop until t becomes 
    %...  larger than t_end. However, to allow t_end to be part of the solution
    %... the integration span is extended to t_end+dt_plot.
        tspan   = [tk t_plot(i_plot):dt_plot:t_end+dt_plot];
    
    %... Initialize step counter
        nsteps = 0;

    %... Integration
        options = odeset('RelTol',Reltol,'AbsTol',Abstol,'JPattern',...
                         jpattern(nz),'Events',@events_StatRegrid);
        sol = ode23s(@ODE_BL_MatMOL_StaticRegridding,tspan,IC,options);

    %... Evaluate solution at elements of tspan lower then integration end
        Y  = (deval(sol,[tspan(tspan<sol.x(end)),sol.x(end)]))';
        X  = Y(end,:);
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
            Profiles(i_plot).S = Y(i1,1:npdes:npdes*nz);
            i_plot  = i_plot+1;
            i1      = i1 + 1;           
        end
            
    %... Refine the grid
        [Z_new,nz_new,ier,tolz_new] = agereg(Z,X,npdes,nzmax,alpha,...
                                             beta,tolz,bound,imesh,ilim);

    %... Interpolate the dependent variables
        X_new = spline(Z,X,Z_new);
        clear Z X

    %... Update the grid, its variables and the dependent variables
        Z    = Z_new';
        nz   = nz_new;
        n    = nz;
        tolz = tolz_new;
        
        IC(1:npdes:npdes*nz) = X_new;
        clear Z_new X_new;
     
    %... Compute new differentiation matrices
        Dzc  = []; Dz = [];
        Dzc = four_point_biased_upwind_D1(Z,1);
        Dz  = three_point_centered_D1(Z);
    
    %... Update loop counter
        i_loop = i_loop+1;
    end
    CPUtime = toc; %... read clock

%... 2e Integrator statisctics
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
        plot(Profiles(i1).z,Profiles(i1).S,'g',Profiles(i1).z,Profiles(i1).S,'.g')
        ylabel('Saturation [-]','Fontsize',18)
        title(['MatMOL Static Regridding: \epsilon=' num2str(eps)],'Fontsize',18)
        axis([zL zR -0.1 1.1])
        set(gca,'Fontsize',16)
        hold on
    
        subplot('position',[0.1 0.1 0.8 0.15])
        plot(Profiles(i1).z,Profiles(i1).t*ones(length(Profiles(i1).z),1),'.g')
        xlabel('Space [-]','fontsize',18);
        ylabel('Time [-]','Fontsize',18)
        axis([zL zR t_start t_end])
        set(gca,'Fontsize',16)
        hold on
       
        figure(2)
        plot(Profiles(i1).t,Profiles(i1).n,'g*');
        xlabel('Time [s]','Fontsize',18)
        ylabel('Nr of grid points [-]','Fontsize',18)    
        set(gca,'Fontsize',18)
        axis([t_start t_end 0 nzmax])
        hold on 
    end

    end %... End of function Buckley_Leverett_MatMOL_StaticRegridding


    function dX = ODE_BL_MatMOL_StaticRegridding(t,X)
%... Provides the right hand side of the ODEs

    global eps
    global Dz Dzc

%... Independent variables
    S=X;

%... Spatial derivatives
    Sz = Dz*S;
    f1  = 4*S.*(1-S).*Sz;
    f1z = Dz*f1;

    f2  = (S.^2)./(S.^2+((1-S).^2));
    Szc = Dzc*f2;

%... Temporal derivative
    St = -Szc + eps*f1z;

%... RHS of ODEs
    dX = [St];
     
    end %... End of function ODE_BL_MatMOL_StaticRegridding


    function S = jpattern(nz)
%... Provides the sparsity pattern of the Jacobian.

    global Dz Dzc

    M   = diag(abs(Dz)*ones(nz,1)) + diag(ones(1,nz))*abs(Dz);
    Jp1 = abs(Dz)*M;
    Jp2 = abs(Dzc)*diag(ones(nz,1));

    S = sparse(sign(Jp1+Jp2)); 

    end %... End of function jpattern

