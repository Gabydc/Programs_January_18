%... The MatMol Group (2016)
     function Buckley_Leverett_MatMOL_DynamicRegridding(eps0,n0,Reltol,Abstol)
%...
%... Function....... : Transient simulation of the classic 1D Buckley-Leverett
%...                   equation for oil reservoirs
%... Technique...... : Solution with the MatMOL toolbox and dynamic regridding
%...                   MOVGRID
%... Input...........: eps0:     capilary pressure term
%...                   n0:       number of vertex centered discritisation points
%...                   RelTol:   Relatitve integration tolerance
%...                   AbsTol:   Absolute integration tolerance

close all

%... Set global variables
     global zL zR n ne X eq choice
     global alpha kappa mu tau
     global u z
     global eps

%... Start a stopwatch timer
     tic

     if nargin < 1, eps0     = 0.001;    end  
     if nargin < 2, n0       = 101;      end  
     if nargin < 3, Reltol   = 1.0e-6;   end  
     if nargin < 4, Abstol   = 1.0e-6;   end  

%... Temporal conditions
     t = [0:0.05:.5];

%... Moving grid parameters
     alpha = 1;
     kappa = 2; 
     mu    = kappa*(kappa+1);
     tau   = .05;

%... PDE parameters
     ne  = 1;
     eps = eps0;
     
%... Initial spatial grid z : if n is the number of moving grid points,
%... dim(z) = n+2
     zL   = 0.0;
     zlim = 1/3;
     zR   = 1.0;
     n    = n0-2;
     dz   = (zR-zL)/(n+1);
     z    = [zL:dz:zR]';

%... Initial dependent variables u(n+2,ne)
     for jj = 1:ne,
        for ii = 1:n+2
            if z(ii) < zlim
                u(ii,jj) = 1-3*z(ii);
            else
               u(ii,jj) = 0 ;
            end
        end
     end

%... Global variables for JPat : function JPat implements the sparsity
%... pattern of the jacobian in X 
     X = zeros(n*(ne+1),n*(ne+1));

%... see right_vect file :
%... diffusion term : stagewise of D1 = 3 pts centered finite difference
%... convection term : kurganov flux-limiter : same pattern as 5 pts
%... centered finite difference
     eq(:,:,1) = [-2 2 1];

%... Monitoring choice :
%... Choose between :
%...                                     1    ne (u(i+1,j)-u(i,j))^2
%...     der1'   : m(i) = sqrt[ alpha + ---  sum ------------------- ]
%...                                     ne  j=1   (z(i+1)-z(i))^2
%...        
%...                                     1    ne     uzz(i+1,j)+uzz(i,j)
%...     'der2'  : m(i) = sqrt[ alpha + ---  sum abs(-------------------)]
%...                                     ne  j=1           2  
%...
%...                                     1    ne (fl[u(i+1,j)]-fl[u(i,j)])^2
%...     'derfl' : m(i) = sqrt[ alpha + ---  sum --------------------------- ]
%...                                     ne  j=1      (z(i+1)-z(i))^2
     choice = 'derfl';

%... Adapt the initial grid 
%... compute the spatial smoothing matrix
     Ainit          = diag((1+2*mu)*ones(n+1,1),0)+diag(-mu*ones(n,1),1)+...
                      diag(-mu*ones(n,1),-1);
     Ainit(1,1)     = Ainit(1,1)-mu;
     Ainit(n+1,n+1) = Ainit(n+1,n+1)-mu;
%... iterative initial grid adaptation
    for iad=1:6
     %... compute the monitoring function mon(i)
        mon = monitor(u,z,choice,t);

     %... compute the new grid
        suminv = 0;
        for ii = 1:n+1,
            suminv = suminv+1/mon(ii);
        end
        for ii = 1:n+1,
            g(ii) = (zR-zL)/(mon(ii)*suminv);
        end
        delz = Ainit\g';
        for ii = 2:n+1,
            z(ii) = z(ii-1)+delz(ii-1);
        end

    %... new initial dependent variables u(n+2,ne)
        for jj = 1:ne,
            for ii = 1:n+2
                if z(ii) < zlim
                    u(ii,jj) = 1-3*z(ii);
                else
                    u(ii,jj) = 0 ;
                end
            end
        end
    end

%... Global initial condition
     for j = 1:ne,
        x(j:ne+1:n*(ne+1)) = u(2:n+1,j);
     end
    x(ne+1:ne+1:n*(ne+1))  = z(2:n+1);

%... call to ODE solver 
    options = odeset('RelTol',Reltol,'AbsTol',Abstol,'Mass',@mass,...
                     'MStateDependence','strong','JPattern',JPat,'MvPattern',...
                     MvPat,'MassSingular','no','stats','on');

    [tout, yout] = ode15s(@right_vect,t,x,options);

%... visualisation of the solution 
    h1 = axes('position',[0.1 0.22 0.8 0.75]);
    h2 = axes('position',[0.1 0.1 0.8 0.10]);
    axes(h1);

%... Separate the dependent variables and the node positions
    absc = zeros(n+2,1);
    ordo = zeros(n+2,ne);

    for ii = 1:length(tout),
        for jj = 1:n,
            absc(jj+1) = yout(ii,(ne+1)*jj);
            for kk = 1:ne,
                ordo(jj+1,kk) = yout(ii,(ne+1)*jj-(ne+1-kk));
            end
        end
        ordo(1,1)   = 1;
        ordo(n+2,1) = 0;
        absc(1)     = zL;
        absc(n+2)   = zR;
    %...  plot the solution
        plot(absc,ordo(:,1),'k')
        hold on
        plot(absc,ordo(:,1),'.k')
    end
    absc = zeros(n+2,1);
    ordo = zeros(n+2,ne);
    xlabel('z');
    ylabel('u(t,z)');
    axis([0 1 -.1 1.1]);

    axes(h2);
    for ii = 1:length(tout),
        for jj = 1:n,
            absc(jj+1) = yout(ii,(ne+1)*jj);
        end
        absc(1)   = zL;
        absc(n+2) = zR;
        plot(absc,((ii-1)*1.)*ones(1,n+2),'.k')
        hold on
    end
    xlabel('z');
    ylabel('t');

%... Read the stopwatch timer
    tcpu=toc