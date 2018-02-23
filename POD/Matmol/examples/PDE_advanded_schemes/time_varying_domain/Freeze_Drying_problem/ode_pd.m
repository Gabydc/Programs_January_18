%... The MatMol Group (2016)
    function dy = ode_pd(t, y, Adv_op_dr, Adv_op_fr, Lap_op_dr, Lap_op_fr,...
                          inv_MMdr, inv_MMfr, ze, ye, ndz, ndy, L, alpha,...
                          Kfr, Kdr, sigma, ep, fp, Tcl, ...
                          rhodr, cpdr, rhofr, cpfr, kd1, kd2,...
                          hL1, hL2, hL3, Pc, Tcr)
    t
%... The states
    Rdr = y(1:ndz);
    Rfr = y(1+ndz:ndz+ndy);
    x   = y(ndz+ndy+1);

%... States gradients
    Rdrz = Adv_op_dr*Rdr;
    Rfry = Adv_op_fr*Rfr;

%... Heat transfer coefficient
    hL  = hL1 + hL2*Pc/(1+Pc/hL3);

%... Front velocity
    w = alpha*(Kfr/(L-x)*Rfry(1) - Kdr/x*Rdrz(ndz));

%... Permeability
    Kd = 1/(kd1*Pc+kd2*x);

%... Front pressure and temperature
    Pfront = Pc + ((rhofr-rhodr)*w)/(Kd);
    Tfront = funcclapeyron(Pfront,'T');

%... Boundary conditions
    Gdr      = zeros(ndz,1);
    Gdr(1)   = sigma*ep*fp*(Tcl^4-Rdr(1)^4)/(x*rhodr*cpdr);
    Gdr(ndz) = 1e6*(Tfront-Rdr(ndz));
    Gfr      = zeros(ndy,1);
    Gfr(1)   = 1e6*(Tfront-Rfr(1));
    Gfr(ndy) = hL*(Tcr-Rfr(ndy))/((L-x)*rhofr*cpfr);
    BVdr     = inv_MMdr*Gdr;
    BVfr     = inv_MMfr*Gfr;

%... Diffusivities
    difdr = Kdr/(rhodr*cpdr);
    diffr = Kfr/(rhofr*cpfr);

%... ODEs
    dRdr = w/x*ze.*Rdrz + (difdr/x^2*Lap_op_dr)*Rdr + BVdr;
    dRfr = w/(L-x)*(1-ye).*Rfry + (diffr/(L-x)^2*Lap_op_fr)*Rfr + BVfr;
    dx   = w;
    dy   = [dRdr; dRfr; dx];
