%...   Solution of the 2D FitzHugh-Nagumo problem using the
%...   finite element method with linear Lagrange elements
%...   FitzHugh-Nagumo
%...
%...   ut = k*(uxx+uyy) + u*(a-u)*(u-1)-v
%...   vt = k*rho*(vxx+vyy) + epsilon*(beta*u-gamma*v+delta)
%...   
%...   0 < x < 200     0 < y < 200
%...
%...   t > 0 
%...
%...   k = 1;   a = 0.1;   rho = 0;  epsilon = 0.01;  beta = 0.5
%...
%...   gamma = 1;  delta = 0
%...
%...   ICs :  u(x,y,0) = 0 for all x > 0;   u(0,y,0) = 1 
%...          v(x,y,0) = 0
%...
%...   BCs :  grad(u(x,y,t)) = grad(v(x,y,t)) = 0 on the boundary

    clear all
    clc
    
%... Global variables
    global alpha epsilon delta gamma beta nd AA
    warning off
  

%... Load the FEM data (previously computed for this problem to save time)
    load data/matrix_fem_1925

%... System parameters
    alpha   = 0.1;
    epsilon = 0.01;
    beta    = 0.5;
    delta   = 0;
    gamma   = 1;

%... FEM spatial operator
    nd = size(MM,1);
    AA = -iMM*DM;

%... Spatial coordinates
    XX      = ms(: , 1);
    YY      = ms(: , 2);
clear ms MM DM iMM 

%... Initial conditions
%... Front
    % kk      = find(XX==0);
    % v0      = zeros(nd , 1);
    % v0(kk)  = 1;
    % w0      = zeros(nd , 1);
%... Spiral
    load data/front_data_201311
    kk     = 180;  % time at which the front breaks
    v0     = vv(:,kk);
    w0     = ww(:,kk);
    kk     = find(YY>=100);
    v0(kk) = 0;
    clear vv ww kk
    y0     = [v0; w0];
    neq    = length(y0);

%... Problem integration
    tf      = 180;
    deltat  = 1;
    tl      = 0 : deltat : tf;
    pasos   = size(tl,2) - 1;
    vv      = zeros(nd, pasos+1);
    ww      = zeros(nd, pasos+1);
    vv(:,1) = v0;
    ww(:,1) = w0;

%... Integration
    resname    = 'ode_fhn';
    jacname    = '[]';              itol       = 1;        rtol       = 1.0e-6; 
    atol       = 1.0e-6;            itask      = 1;        istate     = 1;             
    mf         = 222;               lrw        = 1e7;      rwork      = zeros(lrw,1);  
    iopt       = 0;                 liw        = 5e5;      iwork      = zeros(liw,1);      

    for kk = 1:pasos
        tl(kk)
        
        %... Integration
        [t , y]    = Lsodes(resname, jacname, neq, y0, tl(kk), tl(kk+1),...
                            itol, rtol, atol, itask, istate, iopt, rwork, lrw,...
                            iwork, liw, mf);
    
         %... New initial conditions
         y0 = y(end,:)';
     
         %... Store solution
         vv(: , kk+1) = y0(1:nd);
         ww(: , kk+1) = y0(1+nd:2*nd);
    end
