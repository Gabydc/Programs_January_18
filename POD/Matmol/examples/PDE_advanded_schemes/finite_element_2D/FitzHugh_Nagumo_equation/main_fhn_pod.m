%...   Solution of the 2D FitzHugh-Nagumo problem using the
%...   proper orthogonal decomposition
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
%...   
    clear all
    clc

%... Global varibles    
    global alpha epsilon delta gamma beta 
    global AA nv nw phi_v phi_w proj_phiv proj_phiw
    warning off

%... Load the FEM and POD basis data
    load data/matrix_fem_1925 % Previously computed for this problem to save time
    load data/pod_basis       % Previously computed for this problem to save time using file podcomputation.m
    clear fem iMM

%... System parameters
    alpha   = 0.1;
    epsilon = 0.01;
    beta    = 0.5;
    delta   = 0;
    gamma   = 1;

%... Number of basis function used
    nv = 80;   % number of POD basis to be used in the field v
    nw = 30;   % number of POD basis to be used in the field w
    phi_v = phi_v(:,1:nv);
    phi_w = phi_w(:,1:nw);

%... FEM spatial operator and projection operators
    nd        = size(MM,1);
    AA        = -phi_v'*DM*phi_v;
    proj_phiv = phi_v'*MM;
    proj_phiw = phi_w'*MM;

%... Spatial coordinates
    XX = ms(: , 1);
    YY = ms(: , 2);
    clear ms MM DM

%... Initial conditions
    % type_sol = 'front';
    type_sol = 'spiral';
    switch type_sol
        case{'front'}
            kk      = find(XX==0);
            v0      = zeros(nd , 1);
            v0(kk)  = 1;
            mv0     = proj_phiv*v0;
            w0      = zeros(nd , 1);
            mw0     = proj_phiw*w0;
        case{'spiral'}
            load data/front_data_201311
            kk      = 180;  % time at which the fron breaks
            v0      = vv(:,kk);
            w0      = ww(:,kk);
            kk      = find(YY>=100);
            v0(kk)  = 0;
            mv0     = proj_phiv*v0;
            mw0     = proj_phiw*w0;
            clear vv ww kk
    end
    my0     = [mv0; mw0];
    neq     = length(my0);

%... Problem integration
    tf      = 320;
    deltat  = 1;    % For comparison with the FEM it must be 1
    tl      = 0 : deltat : tf;
    pasos   = size(tl,2) - 1;
    mv      = zeros(nv, pasos+1);
    mw      = zeros(nw, pasos+1);
    mv(:,1) = mv0;
    mw(:,1) = mw0;

%... Integration
    resname    = 'ode_fhn_pod';
    jacname    = '[]';              itol       = 1;        rtol       = 1.0e-6; 
    atol       = 1.0e-6;            itask      = 1;        istate     = 1;             
    mf         = 222;               lrw        = 1e7;      rwork      = zeros(lrw,1);  
    iopt       = 0;                 liw        = 5e5;      iwork      = zeros(liw,1);      

    tic
    for kk = 1:pasos
        tl(kk)
        
        %... Integration
        [t , my]    = Lsodes(resname, jacname, neq, my0, tl(kk), tl(kk+1),...
                             itol, rtol, atol, itask, istate, iopt, rwork, lrw,...
                             iwork, liw, mf);
    
         %... New initial conditions
         my0 = my(end,:)';
     
         %... Store solution
         mv(: , kk+1) = my0(1:nv);
         mw(: , kk+1) = my0(1+nv:nv+nw);
    end
    toc

%... Recovery the original fields
    vv_pod = phi_v*mv;
    ww_pod = phi_w*mw;

%... Comparison with the FEM solution
%... Load the FEM data
    switch type_sol
        case{'front'}
            load data/front_data_201311
            %... tl is also loaded so we recompute it
            tl      = 0 : deltat : tf;
        case{'spiral'}
            load data/spiral_data_201311
            %... tl is also loaded so we recompute it
            tl      = 0 : deltat : tf;
    end
    vv_fem = vv(:,1:tf+1);  % Same length as vv_pod
    mv_fem = proj_phiv*vv_fem;
    ww_fem = ww(:,1:tf+1);  % Same length as ww_pod
    mw_fem = proj_phiw*ww_fem;
    clear vv ww
%... Error computation
    error_v = abs(vv_pod-vv_fem);
    error_w = abs(ww_pod-ww_fem);
    [Mev_v, Mev_v_mc] = max(error_v);
    [Me_v, Me_v_mc]   = max(Mev_v);
    [Mev_w, Mev_w_mc] = max(error_w);
    [Me_w, Me_w_mc]   = max(Mev_w);
    fprintf('Max error v\t Field value v\n')
    fprintf('%2.4f    \t %2.4f\n\n',Me_v, vv_fem(Mev_v_mc(Me_v_mc),Me_v_mc))
    fprintf('Max error w\t Field value w\n')
    fprintf('%2.4f    \t %2.4f\n\n',Me_w, ww_fem(Mev_w_mc(Me_w_mc),Me_w_mc))


%... Plot the solution
    figure
    plot(tl,mv_fem(1:3,:))
    hold on
    plot(tl(1:5:end),mv(1:3,1:5:end),'*')