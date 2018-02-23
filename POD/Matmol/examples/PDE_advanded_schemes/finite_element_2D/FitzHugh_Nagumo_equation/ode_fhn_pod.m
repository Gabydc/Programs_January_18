    function dmy = ode_fhn_pod(t,my)

    global alpha epsilon delta gamma beta 
    global AA nv nw phi_v phi_w proj_phiv proj_phiw

%... The time dependent coefficients
    mv  = my(1:nv , 1);
    mw  = my(1+nv:nv+nw , 1);

%... The field
    v = phi_v*mv;
    w = phi_w*mw;

%... Nonlinear terms
    ff  = proj_phiv*((alpha-v).*(v-1).*v - w);
    gg  = proj_phiw*(epsilon*(beta*v-gamma*w-delta));

%... ODEs
    dmv = AA*mv + ff;
    dmw = gg;
    dmy = [dmv;dmw];