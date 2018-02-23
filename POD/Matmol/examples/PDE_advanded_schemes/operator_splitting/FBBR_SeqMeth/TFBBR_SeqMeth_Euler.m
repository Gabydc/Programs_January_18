%... The MatMol Group (2016)
    function TFBBR_SeqMeth_Euler(Pe,n0,orderD2,sequence)

%... Function....... : Transient simulation of a tubular biochemical reactor
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "M. Laabissi, J.J. Winkin, D. Dochain and M.E. Achhab 2005. 
%...                   Dynamical analysis of a tubular biochemical reactor 
%...                   infinite-dimensional nonlinear model" Procs of the 44th 
%...                   IEEE Conference on Decision and Control & the European 
%...                   Control Conference, 5965-5970.
%... Technique...... : Solution with the Sequencing Method without using
%...                   Matlab's ODE suite for the reaction part but by 
%...                   approximating it by an explicit Euler scheme.
%
%... Input...........: Pe:       Peclet number
%...                   n0:       number of cell centered discritisation points 
%...                   orderD2:  approximation order of the second order 
%...                             spatial derivatives
%...                   sequence: sequence in which the different phenomena are
%...                             applied within one time step
%...                             'CRD' convection-reaction-diffusion
%...                             'CDR' convection-diffusion-reaction

    close all

    global D1 v Sin mu0 k0 kd KS KI
    global n

    if nargin < 1, Pe       = 100.0;    end  
    if nargin < 2, n0       = 100;      end  
    if nargin < 3, orderD2  = 2;        end
    if nargin < 4, sequence = 'CRD';    end
 
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
                X     	      = [(resultX1(:,i));(resultX2(:,i))];
                dX            = ODE_SM(X);
                resultX1(:,i) = resultX1(:,i)+dX(1:n,1)*dt;
                resultX2(:,i) = resultX2(:,i)+dX(n+1:2*n,1)*dt;
                clear  X dX  
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
                X     	      = [(resultX1(:,i));(resultX2(:,i))];
                dX            = ODE_SM(X);
                resultX1(:,i) = resultX1(:,i)+dX(1:n,1)*dt;
                resultX2(:,i) = resultX2(:,i)+dX(n+1:2*n,1)*dt;
                clear  X dX         
            end;
        otherwise
            disp('Choose CRD for 1) Convection  2) Reaction 3) Diffusion, or ')
            disp('CDR for 1) Convection  2) Diffusion 3) Reaction')
    end
%... Reactor properties
    S = resultX1;
    X = resultX2;

%... Statistics
    Stats.nsteps    = t_end/dt;
    Stats.nfailed   = 0;
    Stats.nfevals   = t_end/dt;
    Stats.npds      = 0;
    Stats.ndecomps  = 0;
    Stats.nsolves   = 0;
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

    end %... End of function TDFR_SeqMeth_Euler


    
    
    function dX = ODE_SM(Y)
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
