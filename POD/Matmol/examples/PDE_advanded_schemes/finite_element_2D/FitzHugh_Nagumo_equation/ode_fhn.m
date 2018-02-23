    function dy = ode_fhn(t,y)

    global alpha epsilon delta gamma beta nd AA

%... The fields
    v  = y(1 : nd , 1);
    w  = y(1 + nd : 2*nd , 1);

%... Nonlinear terms
    ff  = (alpha-v).*(v-1).*v - w;

%... ODEs
    dv = AA*v + ff;
    dw = epsilon*(beta*v-gamma*w-delta);
    dy = [dv;dw];