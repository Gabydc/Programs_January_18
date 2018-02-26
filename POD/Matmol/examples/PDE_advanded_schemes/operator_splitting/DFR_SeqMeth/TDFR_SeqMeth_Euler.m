%... The MatMol Group (2016)									
     function TDFR_SeqMeth_Euler(Pe,n0,zswitch,orderD2,sequence)

%... Function....... : Transient simulation of a jacketed tubular reactor with
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "F. Logist, I. Smets, A. Vande Wouwer, J. Van Impe 2005.
%...                   Optimal control of dispersive tubular chemical reactors:
%...                   Part II" Procs of 16th IFAC World Congress, DVD-ROM, 6p.
%... Technique...... : Solution with the Sequencing Method without using
%...                   Matlab's ODE suite for the reaction part but by 
%...                   approximating it by an explicit Euler scheme.
%
%... Input...........: Pe:       Peclet number
%...                   n0:       number of cell centered discritisation points 
%...                   zswitch:  position where the jacket temperature switches 
%...                             from its maximum to its minimum value 
%...                   orderD2:  approximation order of the second order 
%...                             spatial derivatives
%...                   sequence: sequence in which the different phenomena are
%...                             applied within one time step
%...                             'CRD' convection-reaction-diffusion
%...                             'CDR' convection-diffusion-reaction

    close all

    global Cin Tin beta delta k0 E R Tw
    global n

    if nargin < 1, Pe       = 100.0;    end  
    if nargin < 2, n0       = 100;      end  
    if nargin < 3, zswitch  = 0.54;     end
    if nargin < 4, orderD2  = 2;        end
    if nargin < 5, sequence = 'CRD';    end
 
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

    beta	=	0.2;                %... [1/s]     dimensionless heat transfer
    delta	=   0.25;               %... [-]       dimensionless heat of reaction

%... Spatial grid
    n   	=   n0;                 %... [-]       # of discretisation points
    dz  	=   L/n;                %... [m]       spatial discretisation length
    Z       =   L/(2*n):L/n:L;      %... [m]       position vector

%... Temporal grid
    t_start =   0.0;                %... [s]       initial time
    t_end   =   20.0;               %... [s]       final time
    dt  	=   dz/v;               %... [s]       temporal discretisation length
    t   	=   t_start:dt:t_end;   %... [s]       time vector
    dt_plot =   1.0;                  %... [s]       time interval for plots


%... 2 Numerical computation (Sequencing Method)
    
%... Results initialisation
    resultX1 = Tin*ones(n,length(t));
    resultX2 = 0.0*ones(n,length(t));
    
%... Input
    X1in    =  Tin.*ones(1,length(t));
    X2in    =  Cin.*ones(1,length(t));
    
    Tw(1:round(z1*n),1)     = Twmax;
    Tw(round((z1*n)+1):n,1) = Twmin;

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
    matexpX2=expm(Dzz*D2*dt);
    
%... Simulation loop
    disp('Sequence within each time step:')
    switch sequence
        case 'CRD'
            disp('1) Convection  2) Reaction 3) Diffusion')
            for i=2:length(t)
            %... Convection
                resultX1(:,i) = [X1in(i-1);resultX1(1:n-1,i-1)];
                resultX2(:,i) = [X2in(i-1);resultX2(1:n-1,i-1)];
            %... Reaction
                X     	      = [(resultX1(:,i));(resultX2(:,i))];
                dX            = ODE_SM(X);
                resultX1(:,i) = resultX1(:,i)+dX(1:n,1)*dt;
                resultX2(:,i) = resultX2(:,i)+dX(n+1:2*n,1)*dt;
                clear  X dX  
            %... Diffusion
                resultX1(:,i) = matexpX1*resultX1(:,i);
                resultX2(:,i) = matexpX2*resultX2(:,i);
            end;
                      
        case 'CDR'
            disp('1) Convection  2) Diffusion 3) Reaction')
            for i=2:length(t)
            %... Convection
                resultX1(:,i) = [X1in(i-1);resultX1(1:n-1,i-1)];
                resultX2(:,i) = [X2in(i-1);resultX2(1:n-1,i-1)];
            %... Diffusion
                resultX1(:,i) = matexpX1*resultX1(:,i);
                resultX2(:,i) = matexpX2*resultX2(:,i);
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
    T = resultX1;
    C = resultX2;

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
        Tinterp(m,:) = interp1(t,T(m,:),tinterp);
        Cinterp(m,:) = interp1(t,C(m,:),tinterp);
    end

%... Plot figures
    figure(1)
    mesh(tinterp,Z,Tinterp);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [s]','Fontsize',18)
    zlabel('Temperature T [K]','Fontsize',18)
    title(['Temperature evolution (Seq Meth: z1=' num2str(z1) ' m, Pe1=' ...
          num2str(Pe1) ', Pe2=' num2str(Pe2) ', n=' num2str(n) ')'],'Fontsize',20)
    axis([t_start t_end 0.0 L 270.0 450.0])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')

    figure(2)
    mesh(tinterp,Z,Cinterp);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [s]','Fontsize',18)
    zlabel('Concentration C [mole/L]','Fontsize',18)
    title(['Concentration evolution (Seq Meth: z1=' num2str(z1) ' m, Pe1=' ...
          num2str(Pe1) ', Pe2=' num2str(Pe2) ', n=' num2str(n) ')'],'Fontsize',20)
    axis([t_start t_end 0.0 L 0.0 0.02])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')
    
    figure(3)
    plot(Z,Tinterp(:,1:length(tinterp)));
    xlabel('Space z [m]','fontsize',18);
    ylabel('Temperature T [K]','Fontsize',18)
    title(['Seq Meth: z1=' num2str(z1) ' m, Pe1=Pe2=' num2str(Pe2) ', n=' ...
          num2str(n)],'Fontsize',18)
    axis([0.0 L 270.0 450.0])
    grid on
    set(gca,'Fontsize',18)
   
    figure(4)
    plot(Z,Cinterp(:,1:length(tinterp)));
    xlabel('Space z [m]','fontsize',18);
    ylabel('Concentration C [mole/K]','Fontsize',18)
    title(['Seq Meth: z1=' num2str(z1) ' m, Pe1=Pe2=' num2str(Pe2) ', n=' ...
           num2str(n)],'Fontsize',18)
    axis([0.0 L -0.001 0.02])
    grid on
    set(gca,'Fontsize',18)

    end %... End of function TDFR_SeqMeth_Euler

%...
    
    function dX = ODE_SM(X)
%... Provide right hand side for reaction ODEs

    global Cin Tin beta delta k0 E R Tw
    global n

%... Independent variables
    T=X(1:n);
    C=X(n+1:2*n);

%... Temporal derivatives
    Tt = Tin/Cin*delta*k0*C.*exp(-E./(R*T)) + beta.*(Tw-T);
    Ct = - k0*C.*exp(-E./(R*T));

    dX	= [Tt;Ct];

end %... End of function ODE_SM