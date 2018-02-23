%... The MatMol Group (2016)
    function TFBBR_MatMOL(Pe,n0,orderD1,orderD2,Reltol,Abstol)

%... Function....... : Transient simulation of a tubular biochemical reactor
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "M. Laabissi, J.J. Winkin, D. Dochain and M.E. Achhab 2005. 
%...                   Dynamical analysis of a tubular biochemical reactor 
%...                   infinite-dimensional nonlinear model" Procs of the 44th 
%...                   IEEE Conference on Decision and Control & the European 
%...                   Control Conference, 5965-5970.
%... Technique...... : Solution with the MatMOL toolbox
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
    S0 =   0.0*ones(length(Z),1);
    X0 = 300.0*ones(length(Z),1);
    IC = [S0; X0];
    clear S0 X0
    
%... 2c Solution of the DAEs
%... options = odeset('Mass',@mass,'MassSingular','yes','Reltol',Reltol,...
%...           'Abstol',Abstol,'JPattern',Jpattern_MatMOL(orderD1,orderD2),...
%...           'Jacobian',@Jacobian_MatMOL);      
    options = odeset('Mass',@mass,'MassSingular','yes','Reltol',Reltol,...
                      'Abstol',Abstol,'JPattern',Jpattern_MatMOL(orderD1,orderD2)); 
    sol = ode15s(@ODE_MatMOL,[t_start:dt_plot:t_end],IC',options);

%... 2d Reactor properties
    t    = t_start:dt_plot:t_end;   %... time points to be plotted
    sol2 = deval(sol,t);            %... evaluate solution at time points
    S    = sol2(1:n,:);
    X    = sol2(n+1:2*n,:);

%... 2e Integrator statisctics
    Stats   = sol.stats; %... solution statistics
    CPUtime = toc;       %... read clock
    clear sol sol2

    disp('Integration statistics:')
    fprintf('  %g successful steps\n', Stats.nsteps);
    fprintf('  %g failed attempts\n',  Stats.nfailed);
    fprintf('  %g function evaluations\n', Stats.nfevals);
    fprintf('  %g partial derivatives\n',  Stats.npds);
    fprintf('  %g LU decompositions\n', Stats.ndecomps);
    fprintf('  %g solutions of linear systems\n', Stats.nsolves);
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

    end %... End of function TFBBR_MatMOL
%**************************************************************************

    function dX = ODE_MatMOL(t,Y)
%... Provides the right hand side of the ODEs

    global D1 v Sin mu0 k0 kd KS KI
    global n Dz Dzz

%... Independent variables
    S=Y(1:n);
    X=Y(n+1:2*n);

%... Reaction rates
    mu = mu0*S./(KS*X + S + 1/KI*S.^2);

%... Spatial derivatives
    Sz = Dz*S;
    Szz = Dzz*S;

%... Temporal derivatives
    St = D1*Szz - v*Sz - k0*mu.*X;
    Xt =               - kd*X + mu.*X;

%... Danckwerts BC at reactor inlet
    St(1) = -D1*Dz(1,:)*S + v*(S(1)-Sin);

%... Danckwerts BC at reactor outlet
    St(n) = Dz(n,:)*S;

%... RHS of ODEs
    dX = [St; Xt];

    end %... End of function ODE_MatMOL
%**************************************************************************

    function M = mass(t,X)
%... Provides the mass matrix M

    global n

    M = diag([0 ones(1,n-2) 0 ones(1,n) ],0);

    end %... End of function mass
%**************************************************************************

    function Jp = Jpattern_MatMOL(orderD1,orderD2)
%... Provides the sparsity pattern of the Jacobian.

    global n

    switch orderD1
        case 1
            Sparsity_D1      = diag(ones(1,n-1),-1)+diag(ones(1,n),0);
            Sparsity_D1(1,2) = 1;
        case 4
            Sparsity_D1      = diag(ones(1,n-3),-3)+diag(ones(1,n-2),-2)+ ...
                               diag(ones(1,n-1),-1)+diag(ones(1,n  ), 0)+ ...
                           diag(ones(1,n-1),+1);
           Sparsity_D1(1:3,1:5)   = ones(3,5);
            Sparsity_D1(3,3)       = 0;
            Sparsity_D1(end,end-4) = 1;
    end

    switch orderD2
        case 2
            Sparsity_D2      = diag(ones(1,n-1),-1)+diag(ones(1,n),0)+ ...
                               diag(ones(1,n-1),+1);
            Sparsity_D2(1,1:4)         = ones(1,4);
            Sparsity_D2(end,end-3:end) = ones(1,4);
        case 4
            Sparsity_D2      = diag(ones(1,n-2),-2)+diag(ones(1,n-1),-1)+ ...
                               diag(ones(1,n  ),-0)+diag(ones(1,n-1),+1)+ ...
                               diag(ones(1,n-2),+2);
            Sparsity_D2(1:2,1:6)             = ones(2,6);
            Sparsity_D2(end-1:end,end-5:end) = ones(2,6); 
    end

    Sparsity = sparse((abs(Sparsity_D1+Sparsity_D2)) > 1d-9);
    Jp       = sparse([Sparsity eye(n); eye(n) eye(n)]);

    end %... End of function Jpattern_MatMOL
%**************************************************************************
    function dFdX = Jacobian_MatMOL(t,Y)
%... Provides the Jacobian of the right hand side of the ODEs
     
    global D1 v Sin mu0 k0 kd KS KI
    global n Dz Dzz

%... Independent variables
    S=Y(1:n);
    X=Y(n+1:2*n);

%... Jacobian of RHS of ODE
    dF1dX1 = D1*Dzz - v*Dz - diag(k0.*mu0.*X./(KS.*X+S+1./KI.*S.^2) ...
             -k0.*mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*(1+2./KI.*S),0);    
    dF1dX2 = -diag(-k0.*mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*KS ... 
             +k0.*mu0.*S./(KS.*X+S+1./KI.*S.^2),0);
    dF2dX1 = +diag(mu0.*X./(KS.*X+S+1./KI.*S.^2) ...
             -mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*(1+2./KI.*S),0);
    dF2dX2 = +diag(-kd -mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*KS ...
             +mu0.*S./(KS.*X+S+1./KI.*S.^2),0);

%... Danckwerts BC at reactor inlet
    dF1dX1(1,1:n) = -D1*Dz(1,:) + v*[1,zeros(1,n-1)];

%... Danckwerts BC at reactor outlet
    dF1dX1(n,1:n) = Dz(n,:);

    dFdX = sparse([dF1dX1 dF1dX2; dF2dX1 dF2dX2]);

    end %... End of function Jacobian_MatMOL
%**************************************************************************