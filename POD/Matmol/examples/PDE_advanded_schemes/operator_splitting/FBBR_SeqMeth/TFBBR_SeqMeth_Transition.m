%... The MatMol Group (2016)
     function TFBBR_SeqMeth_Transition(Pe,n0,orderD2,sequence)

%... Function....... : Transient simulation of a tubular biochemical reactor
%...                   with dispersion and Danckwerts boundary conditions.
%...                   "M. Laabissi, J.J. Winkin, D. Dochain and M.E. Achhab 2005. 
%...                   Dynamical analysis of a tubular biochemical reactor 
%...                   infinite-dimensional nonlinear model" Procs of the 44th 
%...                   IEEE Conference on Decision and Control & the European 
%...                   Control Conference, 5965-5970.
%... Technique...... : Solution with the Sequencing Method without using 
%...                   Matlab's ODE suite for the reaction part but by 
%...                   approximating it by a transition matrix formulation .
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
    global dt

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
                for j=1:n
                    x1h             =   resultX1(j,i);
                    x2h             =   resultX2(j,i);
                    Y               =   Reaction_SM(x1h,x2h);
                    resultX1(j,i)   =   x1h+Y(1);
                    resultX2(j,i)   =   x2h+Y(2);
                end            
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
                for j=1:n
                    x1h             =   resultX1(j,i);
                    x2h             =   resultX2(j,i);
                    Y               =   Reaction_SM(x1h,x2h);
                    resultX1(j,i)   =   x1h+Y(1);
                    resultX2(j,i)   =   x2h+Y(2);
                end            
            end;
        otherwise
            disp('Choose CRD for 1) Convection  2) Reaction 3) Diffusion, or ')
            disp('CDR for 1) Convection  2) Diffusion 3) Reaction')
    end
%... Reactor properties
    S = resultX1;
    X = resultX2;

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

    end %... End of function TDFR_SeqMeth_Transition


%...

    function Y=Reaction_SM(x1h,x2h)
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

    global D1 v Sin mu0 k0 kd KS KI
    global dt

    S   = x1h;
    X   = x2h;
    mu  = mu0* S./(KS*X + S + 1/KI*S.^2);

%n1      = 1;           
    Fa     = - k0*mu.*X;
    dFadx1 = -(k0.*mu0.*X./(KS.*X+S+1./KI.*S.^2) ...
             -k0.*mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*(1+2./KI.*S));    
    dFadx2 = -(-k0.*mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*KS ... 
             +k0.*mu0.*S./(KS.*X+S+1./KI.*S.^2));
    Fb     = - kd*X + mu.*X;
    dFbdx1 = +(mu0.*X./(KS.*X+S+1./KI.*S.^2) ...
             -mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*(1+2./KI.*S));
    dFbdx2 = +(-kd -mu0.*S.*X./(KS.*X+S+1./KI.*S.^2).^2.*KS ...
             +mu0.*S./(KS.*X+S+1./KI.*S.^2));
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