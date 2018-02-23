%... The MatMol Group (2016)
    function dmx = pde_brusselator_pod(t, mx, neig_u, neig_v, a, b, q,...
                                       phiu, phiv, proj_op_u, proj_op_v,...
                                       plapop_u, plapop_v, nd, u_inf, v_inf)

%... Mode separation
    mu = mx(1:neig_u,1);
    mv = mx(1+neig_u:neig_u+neig_v,1);

%... Field computation
    u = phiu*mu;
    v = phiv*mv;

%... Boundary conditions
    Gu(1,1)  = q*(u_inf-u(1));
    Gu(nd,1) = q*(u_inf-u(nd));
    Gv(1,1)  = q*(v_inf-v(1));
    Gv(nd,1) = q*(v_inf-v(nd));
    iGu      = phiu'*Gu;
    iGv      = phiv'*Gv;

%... Nonlinear functions
    nlu  = a - (b + 1)*u + u.^2.*v;
    nlv  = b*u - u.^2.*v;
    pnlu = proj_op_u*nlu;
    pnlv = proj_op_v*nlv;

%... PDE
    dmu = plapop_u*mu + pnlu + iGu;
    dmv = plapop_v*mv + pnlv + iGv;
    dmx = [dmu; dmv];