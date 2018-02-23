%... The MatMol Group (2016)									
     function TDFR_SeqMeth_Transition(Pe,n0,zswitch,orderD2,sequence)
%... Function....... : Transient simulation of a jacketed tubular reactor with
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "F. Logist, I. Smets, A. Vande Wouwer, J. Van Impe 2005.
%...                   Optimal control of dispersive tubular chemical reactors:
%...                   Part II" Procs of 16th IFAC World Congress, DVD-ROM, 6p.
%... Technique...... : Solution with the Sequencing Method without using 
%...                   Matlab's ODE suite for the reaction part but by 
%...                   approximating it by a transition matrix formulation.
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
%**************************************************************************
    close all

    global Cin Tin beta delta k0 E R
    global dt

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
    dt_plot =   1.0;                %... [s]       time interval for plots


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
                u = Twmax; %... max part [0,zswitch]   
                for j=1:round(n*z1)
                    x1h             =   resultX1(j,i);
                    x2h             =   resultX2(j,i);
                    Y               =   Reaction_SM(x1h,x2h,u);
                    resultX1(j,i)   =   x1h+Y(1);
                    resultX2(j,i)   =   x2h+Y(2);
                end
                u = Twmin; %...  min part ]zswitch,L]
                for j=round(n*z1)+1:n
                    x1h             =   resultX1(j,i);
                    x2h             =   resultX2(j,i);
                    Y               =   Reaction_SM(x1h,x2h,u);
                    resultX1(j,i)   =   x1h+Y(1);
                    resultX2(j,i)   =   x2h+Y(2);
                end   
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
                u = Twmax; %... max part [0,zswitch]   
                for j=1:round(n*z1)
                    x1h             =   resultX1(j,i);
                    x2h             =   resultX2(j,i);
                    Y               =   Reaction_SMd(x1h,x2h,u);
                    resultX1(j,i)   =   x1h+Y(1);
                    resultX2(j,i)   =   x2h+Y(2);
                end
                u = Twmin; %...  min part ]zswitch,L]
                for j=round(n*z1)+1:n
                    x1h             =   resultX1(j,i);
                    x2h             =   resultX2(j,i);
                    Y               =   Reaction_SMd(x1h,x2h,u);
                    resultX1(j,i)   =   x1h+Y(1);
                    resultX2(j,i)   =   x2h+Y(2);
                end   
            end;
        otherwise
            disp('Choose CRD for 1) Convection  2) Reaction 3) Diffusion, or ')
            disp('CDR for 1) Convection  2) Diffusion 3) Reaction')
    end

%... Reactor properties
    T   = resultX1;
    C   = resultX2;

%... Statistics
    Stats.nsteps    = t_end/dt*n;
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
    title(['Temperature evolution (Seq Meth without ODE: z1=' num2str(z1),... 
           'm, Pe1=' num2str(Pe1) ', Pe2=' num2str(Pe2) ', n=' num2str(n) ')'],...
           'Fontsize',20)
    axis([t_start t_end 0.0 L 270.0 450.0])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')

    figure(2)
    mesh(tinterp,Z,Cinterp);
    ylabel('Space z [m]','Fontsize',18)
    xlabel('Time t [s]','Fontsize',18)
    zlabel('Concentration C [mole/L]','Fontsize',18)
    title(['Concentration evolution (Seq Meth without ODE: z1=' num2str(z1),...
            ' m, Pe1=' num2str(Pe1) ', Pe2=' num2str(Pe2) ', n=' num2str(n) ')'],...
            'Fontsize',20)
    axis([t_start t_end 0.0 L 0.0 0.02])
    grid on
    set(gca,'Fontsize',18,'Ydir','rev')
    
    figure(3)
    plot(Z,Tinterp(:,1:length(tinterp)));
    xlabel('Space z [m]','fontsize',18);
    ylabel('Temperature T [K]','Fontsize',18)
    title(['Seq Meth without ODE: z1=' num2str(z1) ' m, Pe1=Pe2=' ...
            num2str(Pe2) ', n=' num2str(n)],'Fontsize',18)
    axis([0.0 L 270.0 450.0])
    grid on
    set(gca,'Fontsize',18)
   
    figure(4)
    plot(Z,Cinterp(:,1:length(tinterp)));
    xlabel('Space z [m]','fontsize',18);
    ylabel('Concentration C [mole/K]','Fontsize',18)
    title(['Seq Meth without ODE: z1=' num2str(z1) ' m, Pe1=Pe2=' ...
            num2str(Pe2) ', n=' num2str(n)],'Fontsize',18)
    axis([0.0 L -0.001 0.02])
    grid on
    set(gca,'Fontsize',18)

    end     %... End of function TDFR_SeqMeth_Transition



    function Y=Reaction_SM(x1h,x2h,uh)
%... Solve the reaction ODE by a transition matrix formulation

%... Rationale behind the method:
%... ****************************

%... 1 for linear systems:
%... Y'   = A Y + B    with Y(0)=0                                   (Eqn (1))
%... General solution is of the form
%... Y(t) = expm(A*t)*(Y(0)+ (int_{0}^{t} expm(-A*xi)*B dxi))        (Eqn (2))

%... 2 Now we want to calculate the solution for x1 and x2 over the time 
%... interval t* until t* + delta t*

%... From the balance equations we know:
%... dx1/dt = Fa(x1,x2)
%... dx2/dt = Fa(x1,x2)

%... Linearize the solution x1(t) and x2(t) for every time step around its 
%... 'initial' state (i.e. the value at the begin of each time interval t*) 
%... x1* and x2* 

%... dx1/dt = Fa(x1*,x2*) + dFa/dx1(x1*,x2*)*(x1-x1*) 
%...                      + dFa/dx2(x1*,x2*)*(x2-x2*)
%... dx2/dt = Fb(x1*,x2*) + dFb/dx1(x1*,x2*)*(x1-x1*) 
%...                      + dFb/dx2(x1*,x2*)*(x2-x2*)

%... and by the introduction of deviation variables y1=x1-x1* and y2=x2-x2*
%... dy1/dt = dFa/dx1(x1*,x2*)*y1 + dFa/dx2(x1*,x2*)*y2 + Fa(x1*,x2*)
%... dy2/dt = dFb/dx1(x1*,x2*)*y1 + dFb/dx2(x1*,x2*)*y2 + Fb(x1*,x2*)

%... Rewriting with 
%... Y=[y1; y2],
%... A=[dFa/dx1(x1*,x2*) dFa/dx2(x1*,x2*); 
%...    dFb/dx1(x1*,x2*) dFb/dx2(x1*,x2*) ],
%... B=[Fa(x1*,x2*); Fb(x1*,x2*)],
%... and Y(0)=0
%... makes this fit perfectly in Eqn (1) for which the solution is given by 
%... Eqn (2). 

%... The solution [x1; x2] is simply is given by [x1*; x2*] + Y(delta t*). The
%... integral in Eqn (2) is computed approximately by employing the trapezium 
%... rule over n+1 points within the interval.

    global Cin Tin beta delta k0 E R
    global dt

%n1      = 1;           
    Fa     = Tin/Cin*delta*k0*x2h.*exp(-E./(R*x1h)) + beta.*(uh-x1h);
    dFadx1 = Tin/Cin*delta*k0*x2h.*exp(-E./(R*x1h)).*E./(R*x1h.^2) - beta;
    dFadx2 = Tin/Cin*delta*k0   *exp(-E./(R*x1h));
    Fb     = - k0*x2h.*exp(-E./(R*x1h));
    dFbdx1 = -k0*x2h.*exp(-E./(R*x1h)).*E./(R*x1h.^2);
    dFbdx2 = -k0  .*exp(-E./(R*x1h));
    A      = [dFadx1 dFadx2; dFbdx1 dFbdx2];
    B      = [Fa;Fb];
%Y0      = [0;0];

%... %... When n1 = 1
%... IntM    = (eye(2) + expm(-A*dt))*B/2*dt;
%... Y       = expm(A*dt)*(Y0+IntM);
%... which can be simplified for the current case to:
    Y = (expm(A*dt)+eye(2))*B/2*dt;

%... %... When n1 >= 1 use trapezium rule in Matlab
%... t_vec = [0:dt/n1:dt];
%... for p=1:length(t_vec)
%...     m_vec(1:2,p) = expm(-A*t_vec(p))*B;
%... end
%... IntM(1:2,1) = trapz(t_vec,m_vec,2);
%... Y           = expm(A*dt)*(Y0+IntM);

    end %... End of function Reaction_SM