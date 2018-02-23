%... The MatMol Group (2016)
     function TFBBR_SeqMeth_Implicit(Pe,n0,orderD2,sequence,Reltol,Abstol)

%... Function....... : Transient simulation of a tubular biochemical reactor
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "M. Laabissi, J.J. Winkin, D. Dochain and M.E. Achhab 2005. 
%...                   Dynamical analysis of a tubular biochemical reactor 
%...                   infinite-dimensional nonlinear model" Procs of the 44th 
%...                   IEEE Conference on Decision and Control & the European 
%...                   Control Conference, 5965-5970.
%... Technique...... : Solution with the Sequencing Method using matlab's 
%...                   implicit integrator ode15s for the reaction part
%
%... Input...........: Pe:       Peclet number
%...                   n0:       number of cell centered discritisation points 
%...                   orderD2:  approximation order of the second order 
%...                             spatial derivatives
%...                   sequence: sequence in which the different phenomena are
%...                             applied within one time step
%...                             'CRD' convection-reaction-diffusion
%...                             'CDR' convection-diffusion-reaction
%...                   RelTol:   Relative integration tolerance
%...                   AbsTol:   Absolute integration tolerance

    close all

    global D1 v Sin mu0 k0 kd KS KI
    global n

    if nargin < 1, Pe       = 100.0;    end  
    if nargin < 2, n0       = 100;      end  
    if nargin < 3, orderD2  = 2;        end
    if nargin < 4, sequence = 'CRD';    end
    if nargin < 5, Reltol   = 1.0e-3;   end  
    if nargin < 6, Abstol   = 1.0e-3;   end
    
    tic %... initialise clock


%... 1 Parameters

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
    n   	=   n0;                 %... [-]       # of discretisation points
    dz  	=   L/n;                %... [m]       spatial discretisation length
    Z       =   L/(2*n):L/n:L;      %... [m]       position vector

%... Temporal grid
    t_start =   0.0;                %... [s]       initial time
    t_end   =   2.0;                %... [s]       final time
    dt  	=   dz/v;               %... [s]       temporal discretisation length
    t   	=   t_start:dt:t_end;   %... [s]       time vector
    dt_plot =   0.1;                %... [s]       time interval for plots


%... 2 Numerical computation (Sequencing Method)
    
%... Results initialisation
    resultX1 =   0.0*ones(n,length(t));
    resultX2 = 300.0*ones(n,length(t));
    
%... Input
    X1in    =  Sin.*ones(1,length(t));
    
%... Diffusion matrix calculation
    switch orderD2
        case 2
            z   = linspace(0+dz/2,L-dz/2,n);
            Dzz = three_point_centered_D2(z);
            Dzz(1,:) = 1/dz^2*[-1 1 zeros(1,n-2)];
            Dzz(n,:) = 1/dz^2*[zeros(1,n-2) 1 -1];        
        case 4
            z   = linspace(0+dz/2,L-dz/2,n);
            Dzz = five_point_centered_D2(z);
            Dzz(1,:) = 1/dz^2*[-1 1 zeros(1,n-2)];
            Dzz(n,:) = 1/dz^2*[zeros(1,n-2) 1 -1];        
        otherwise
            disp('Invalid selection for order of second spatial derivative.')
            disp('Choose 2 for low-order or 4 for high-order scheme.')
    end

%... Exponential matrix calculation
    matexpX1=expm(Dzz*D1*dt);

%... Initialisation statistics
    Stats.nsteps   = 0;
    Stats.nfailed  = 0;
    Stats.nfevals  = 0;
    Stats.npds     = 0;
    Stats.ndecomps = 0;
    Stats.nsolves  = 0;

%... Simulation loop
    disp('Sequence within each time step:')
    switch sequence
        case 'CRD'
            disp('1) Convection  2) Reaction 3) Diffusion')
            for i=2:length(t)
            %... Convection
                resultX1(:,i) = [X1in(i-1);resultX1(1:n-1,i-1)];
                resultX2(:,i) = [          resultX2(:    ,i-1)];
            %... Reaction
                options	      = odeset('RelTol',Reltol,'AbsTol',Abstol,...
                                        'Jpattern',@Jacobian_SM_pattern,...
                                        'Jacobian',@Jacobian_SM,'maxstep',dt); 
                IC     	      = [(resultX1(:,i))',(resultX2(:,i))'];
                sol           = ode15s(@ODE_SM,[0:dt:dt],IC,options);
                resultX1(:,i) = (sol.y(1:n,end))';
                resultX2(:,i) = (sol.y(n+1:2*n,end))';
                Stats         = AddIntStats(Stats,sol.stats,sol.solver);
                clear  sol IC
            %... Diffusion
                resultX1(:,i) = matexpX1*resultX1(:,i);
                resultX2(:,i) =          resultX2(:,i);
            end;
                      
        case 'CDR'
            disp('1) Convection  2) Diffusion 3) Reaction')
            for i=2:length(t)
            %... Convection
                resultX1(:,i) = [X1in(i-1);resultX1(1:n-1,i-1)];
                resultX2(:,i) = [          resultX2(:    ,i-1)];
            %... Diffusion
                resultX1(:,i) = matexpX1*resultX1(:,i);
                resultX2(:,i) =          resultX2(:,i);
            %... Reaction
                options	      = odeset('RelTol',Reltol,'AbsTol',Abstol,...
                                        'Jpattern',@Jacobian_SM_pattern,...
                                        'Jacobian',@Jacobian_SM,'maxstep',dt); 
                IC     	      = [(resultX1(:,i))',(resultX2(:,i))'];
                sol           = ode15s(@ODE_SM,[0:dt:dt],IC,options);
                resultX1(:,i) = (sol.y(1:n,end))';
                resultX2(:,i) = (sol.y(n+1:2*n,end))';
                Stats         = AddIntStats(Stats,sol.stats,sol.solver);
                clear  sol IC
            end;
        otherwise
            disp('Choose CRD for 1) Convection  2) Reaction 3) Diffusion, or ')
            disp('CDR for 1) Convection  2) Diffusion 3) Reaction')
    end
%... Reactor properties
    S = resultX1;
    X = resultX2;

%... Statistics
    CPUtime         = toc;       %... read clock

    disp('Statistics:')
    fprintf('  %g successful steps\n', Stats.nsteps);
    fprintf('  %g failed attempts\n',  Stats.nfailed);
    fprintf('  %g function evaluations\n', Stats.nfevals);
    fprintf('  %g partial derivatives\n',  Stats.npds);
    fprintf('  %g LU decompositions\n', Stats.ndecomps);
    fprintf('  %g solutions of linear systems\n', Stats.nsolves);
    fprintf('  %g s (total) computation time\n',  CPUtime);

%... 3 Plot results

%... Interpolate
    tinterp = t_start:dt_plot:t_end;
    for m=1:n
        Sinterp(m,:) = interp1(t,S(m,:),tinterp);
        Xinterp(m,:) = interp1(t,X(m,:),tinterp);
    end

%... Plot figures
    figure(1)
    mesh(tinterp,Z,Sinterp);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [h]','Fontsize',18)
    zlabel('Substrate S [g/L]','Fontsize',18)
    title(['Substrate evolution (Seq Meth: Pe=' ...
          num2str(Pe1) ', n=' num2str(n) ')'],'Fontsize',20)
    axis([t_start t_end 0 L -1.0 21.0])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')

    figure(2)
    mesh(tinterp,Z,Xinterp);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [h]','Fontsize',18)
    zlabel('Biomass X [g/L]','Fontsize',18)
    title(['Biomass evolution (Seq Meth: Pe=' ...
           num2str(Pe1) ', n=' num2str(n) ')'],'Fontsize',20)
    axis([t_start t_end 0 L 295.0 305.0])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')
    
    figure(3)
    plot(Z,Sinterp(:,1:length(tinterp)));
    xlabel('Space z [m]','fontsize',18);
    ylabel('Substrate S [g/L]','Fontsize',18)
    title(['Seq Meth: Pe=' num2str(Pe1) ', n=' num2str(n)],'Fontsize',18)
    axis([0 L -1.0 21.0])
    grid on
    set(gca,'Fontsize',18)
   
    figure(4)
    plot(Z,Xinterp(:,1:length(tinterp)));
    xlabel('Space z [m]','fontsize',18);
    ylabel('Biomass X [g/L]','Fontsize',18)
    title(['Seq Meth: Pe=' num2str(Pe1) ', n=' num2str(n)],'Fontsize',18)
    axis([0 L 295.0 305.0])
    grid on
    set(gca,'Fontsize',18)

    end %... End of function TFBBR_SeqMeth_Implicit

    
    
    
    function dX = ODE_SM(t,Y)
%... Provide right hand side for reaction ODEs

    global D1 v Sin mu0 k0 kd KS KI
    global n

%... Independent variables
    S=Y(1:n);
    X=Y(n+1:2*n);

%... Reaction rates
    mu = mu0* S./(KS*X + S + 1/KI*S.^2);

%... Temporal derivatives
    St = - k0*mu.*X;
    Xt = - kd*X + mu.*X;

    dX	= [St;Xt];

    end %... End of function ODE_SM

    
    
    
    function Jpat = Jacobian_SM_pattern
%... Provide sparsity pattern of Jacobian of right hand side for reaction ODEs

    global n

    Jpat = sparse(repmat(diag(ones(1,n)),2,2));

    end %... End of function Jacobian_SM_pattern

    
    
    
    
    function dFdX = Jacobian_SM(t,Y)
%... Provide Jacobian of right hand side for reaction ODEs

    global D1 v Sin mu0 k0 kd KS KI
    global n

%... Independent variables
    S=Y(1:n);
    X=Y(n+1:2*n);

%... Jacobian of RHS of ODE
    dF1dX1 = -(k0.*mu0.*X./(KS.*X+S+1./KI.*S.^2) ...
             -k0.*mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*(1+2./KI.*S));    
    dF1dX2 = -(-k0.*mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*KS ... 
             +k0.*mu0.*S./(KS.*X+S+1./KI.*S.^2));
    dF2dX1 = +(mu0.*X./(KS.*X+S+1./KI.*S.^2) ...
             -mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*(1+2./KI.*S));
    dF2dX2 = +(-kd -mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*KS ...
             +mu0.*S./(KS.*X+S+1./KI.*S.^2));

    dFdX   = sparse([diag(dF1dX1) diag(dF1dX2); diag(dF2dX1) diag(dF2dX2)]);

    end %... End of function Jacobian_SM
