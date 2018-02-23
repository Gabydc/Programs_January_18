%... The MatMol Group (2016)
%...
%... Shock-tube problem

%...   u1  = -f1  + eps*u1
%...     t      z         zz
%...
%...   u2  = -f2  + eps*u2
%...     t      z         zz
%...
%...   u3  = -f3  + eps*u3
%...     t      z         zz
%...
%... with
%...
%...   f1 = u2
%...
%...                     (gam-3)*(u2^2)
%...   f2 = (gam-1)*u3 - --------------
%...                           2*u1
%...
%...                  (gam-1)*(u2^2)  u2
%...   f3 = [gam*u3 - --------------]*--
%...                        2*u1      u1
%...
%...
%... where
%...
%... u1, u2 u3 dependent variable
%...
%... t         time
%...
%... z         space
%...
%... eps       diffusion coefficient
%...
%... gam       isentropic coefficient
%...
%...
%... I.C. : u1(z,0) = 1      0 < z < 0.5
%...                  0.125  0.5 < z < 1
%...
%...        u2(z,0) = 0      O < z < 1
%...
%...        u3(z,0) = 2.5    0 < z < 0.5
%...                  0.25   0.5 < z < 1
%...
%... B.C. : u1(0,t)  = 0        u1(1,t) = 0
%...               z                   z
%...
%...        u2(0,t) = 0         u2(1,t) = 0
%...
%...        u3(O,t) = 0         u3(1,t)  = 0 
%...               z                   z
%...
%... numerical values :
%...
%... gam = 1.4
%...
%... eps = 0.1, 0.01, ...
%...
%... Implementation of MOVGRID - Finite differences version
%...
    close all
    clear all
%...
%... start a stopwatch timer
%...
    tic
%... 
%... set global variables
%...
    global zL zR n ne X eq choice
    global alpha kappa mu tau
    global u z
    global gam eps
%...
%... temporal conditions
    t = [0 .15 .23 .28 .32];
%...
%... moving grid parameters
%...
    alpha = 0.05;
    kappa = 1;
    mu    = kappa*(kappa+1);
    tau   = 1e-4;
%...
%... PDE parameters
    ne   = 3;
    eps  = 10^(-3);
    gam  = 1.4;
%...
%... initial spatial grid z : if n is the number of moving grid points,
%... dim(z) = n+2
    zL  = 0;
    zLL = 0.5-15*eps;
    zmL = 0.5-5*eps;
    zmR = 0.5+5*eps;
    zRR = 0.5+15*eps;
    zR  = 1;
    nLL = 15;
    nL  = 15;
    nm  = 60;
    nR  = 25;
    nRR = 25;
    dz1 = (zLL-zL)/(nLL-1);
    dz2 = (zmL-zLL)/(nL-1);
    dz3 = (zmR-zmL)/(nm-1);
    dz4 = (zRR-zmR)/(nR-1);
    dz5 = (zR-zRR)/(nRR-1);
    z   = [zL:dz1:zLL zLL+dz2:dz2:zmL zmL+dz3:dz3:zmR zmR+dz4:dz4:zRR zRR+dz5:dz5:zR]';
    n   = length(z)-2;
%...
%... initial dependent variables u(n+2,ne)
%...
    for i=1:n+2,
        if z(i)<=zmL
            u(i,1)=1;
            u(i,2)=0;
            u(i,3)=2.5;
        elseif z(i)<zmR
            u(i,1)=1+(.125-1)*(z(i)-zmL)/(zmR-zmL);
            u(i,2)=0;
            u(i,3)=2.5+(.25-2.5)*(z(i)-zmL)/(zmR-zmL);
        else
            u(i,1)=.125;
            u(i,2)=0;
            u(i,3)=.25;
        end
    end
%...
%... global variables for JPat : function JPat implements the sparsity
%... pattern of the jacobian in X ; informations depending on the
%... choosen spatial derivative matrices are stocked in the ne eq(ne,ne,3) matrices.
%... 
%...   example : hypothetic problem : 
%...
%...  u1  = -f1  + eps1*u1   + k1*u2    where  f1 depends on u3
%...    t      z          zz
%...
%...  u2  = -f2  + eps2*u2              where  f2 depends on u2 and u3
%...    t      z          zz
%...
%...  u3  = -f3  + eps3*u3              where  f3 depends on u1, u2 and u3
%...    t      z          zz
%...
%...   where the 1st (2nd) derivative is computed by a 3 (5) pts centered finite difference
%... 
%...   It results :
%...
%...  eq(:,:,1) = [-2 2 1; 0 0 1;-1 1 1];
%...  eq(:,:,2) = [ 0 0 0;-2 2 1;-1 1 1];
%...  eq(:,:,3) = [-1 1 1;-1 1 1;-2 2 1];

    X = zeros(n*(ne+1),n*(ne+1));
%... 
%... see right_vect file :
%... convection terms : kurganov flux limiter : same pattern as 5 pts centered finite difference
%... diffusion terms : 5 pts centered finite difference
    eq(:,:,1) = [-2 2 1;-2 2 1; 0 0 0];
    eq(:,:,2) = [-2 2 1;-2 2 1;-2 2 1];
    eq(:,:,3) = [-2 2 1;-2 2 1;-2 2 1];
%...
%... monitoring choice
%...
%...  1    ne (u(i+1,j)-u(i,j))^2
%...   choose between : 'der1'  : m(i) = sqrt[ alpha + ---  sum ------------------- ]
%...  ne  j=1   (z(i+1)-z(i))^2
%...
%...  
%...
%...  1    neuzz(i+1,j)+uzz(i,j)
%... 'der2'  : m(i) = sqrt[ alpha + ---  sum abs(-------------------)]
%...  ne  j=1 2  
%...
%...
%...
%...  1    ne (fl[u(i+1,j)]-fl[u(i,j)])^2
%... 'derfl' : m(i) = sqrt[ alpha + ---  sum --------------------------- ]
%...  ne  j=1 (z(i+1)-z(i))^2
%...

    choice = 'der1';
%...
%... 
%... Global initial condition
%...
    for j=1:ne,
        x(j:ne+1:n*(ne+1)) = u(2:n+1,j);
    end
    x(ne+1:ne+1:n*(ne+1))  = z(2:n+1);
%...
%... %%%%%%%%%%%%%%%%%%%%%%
%... %... call to ODE solver %
%... %%%%%%%%%%%%%%%%%%%%%%
%...
    options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Mass',@mass,'MStateDependence','strong',...
                     'JPattern',JPat,'MvPattern',MvPat,'MassSingular','no','stats','on');
%...
    [tout, yout] = ode15s(@right_vect,t,x,options);
%...
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%... %... visualisation of the solution %
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
%... separate the dependent variables and the node positions
%...
    h1=axes('position',[0.05 0.55 0.4 0.4]);
    h2=axes('position',[0.55 0.55 0.4 0.4]);
    h3=axes('position',[0.05 0.05 0.4 0.4]);
    h4=axes('position',[0.55 0.05 0.4 0.4]);

    absc=zeros(n+2,1);
    ordo=zeros(n+2,2);
%...
    for i=1:length(tout),
        for j=1:n,
            absc(j+1)=yout(i,(ne+1)*j);
            for k=1:ne;
                ordo(j+1,k)=yout(i,(ne+1)*j-(ne+1-k));
            end
        end
        ordo(1,1)=ordo(2,1);
        ordo(1,2)=ordo(2,2);
        ordo(1,3)=ordo(2,3);
        ordo(n+2,1)=ordo(n+1,1);
        ordo(n+2,3)=ordo(n+1,3);
        absc(1)=zL;
        absc(n+2)=zR;
%...
%... plot the solution
%...
        axes(h1);
%...
        plot(absc,ordo(:,1),'g');
        hold on;
        xlabel('z');
        ylabel('u');
        axis([0 1 0 1.1]);
%...
        axes(h2);
%...
        plot(absc,ordo(:,2),'g');
        hold on;
        xlabel('z');
        ylabel('v');
        axis([0 1 0 .5]);
%...
        axes(h3);
%...
        plot(absc,ordo(:,3),'g');
        hold on;
        xlabel('z');
        ylabel('w');
        axis([0 1 0 2.6]);
    end
%...
    axes(h4);
    for i=1:length(tout),
        for j=1:n,
            absc(j+1)=yout(i,(ne+1)*j);
        end
        absc(1)=zL;
        absc(n+2)=zR;
        plot(absc,((i-1)*.05)*ones(1,n+2),'.k')
        hold on
    end

xlabel('z');
ylabel('t');
%
%... read the stopwatch timer
%
tcpu=toc    
%... 
