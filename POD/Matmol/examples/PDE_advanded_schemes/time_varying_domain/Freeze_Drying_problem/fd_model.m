%... The MatMol Group (2016)
%... Solution of the 1D Freeze-Drying problem using the Landau
%... transform to fix the moving boundary and the finite
%... element method with linear Lagrange elements
%...
%... Freeze-Drying model
%...   dRdr = w/x*z.*Rdrz + (difdr/x^2)*Rdrzz
%...   dRfr = w/(L-x)*(1-y).*Rfry + (diffr/(L-x)^2)*Rfryy
%...   dx   = w
%... Heat transfer coefficient
%...   hL  = hL1 + hL2*Pc/(1+Pc/hL3)
%... 
%... Front velocity
%...   w = alpha*(Kfr/(L-x)*Rfry(1) - Kdr/x*Rdrz(ndz))
%...
%... Permeability
%...   Kd = 1/(kd1*Pc+kd2*x)
%... 
%... Front pressure
%...   Pfront = Pc + ((rhofr-rhodr)*w)/(Kd)
%... 
%... Boundary conditions
%...   Rdrz(1)   = sigma*ep*fp*(Tcl^4-Rdr(1)^4)/(x*rhodr*cpdr)
%...   Rdr(ndz)  = Tfront
%...   Rfr(1)    = Tfront
%...   Rfry(ndy) = hL*(Tcr-Rfr(ndy))/((L-x)*rhofr*cpfr);
%... Diffusivities
%...   difdr = Kdr/(rhodr*cpdr);
%...   diffr = Kfr/(rhofr*cpfr);

    clear all
    clc

%... Fixed parameters
%... Dried region parameters
    ep    = 0.78;            % Emisivity on the product top [-]
    fp    = 0.98;            % Geometrical correction factor on the top [-]
    rhodr = 200.31;          % Density [kg m^-3]
    cpdr  = 1254;            % Specific heat [J kg^-1 K^-1]
    Kdr   = 0.0129;      % Heat conductivity in the dried region
    ndz   = 13;              % Dried region discretization points
%... Frozen region parameters
    Kfr   = 2.4;             % Heat conductivity [W K^-1 m^-1]
    rhofr = 1001.6;          % Density [kg m^-3]
    cpfr  = 1818.8;          % Specific heat [J kg^-1 K^-1]
    ndy   = 13;              % Dried region discretization points
%... Other parameters
    TK      = 273.15;          % Transformation Celsius <--> Kelvin
    L       = 5.75e-3;         % Sample length [m]
    DHs     = 2791.2e3;        % Sublimation heat [J kg^-1]
    sigma   = 5.6704e-8;       % Stefan-Boltzmann constant [W m^-2 K^-4]
    alpha   = 1/((rhofr-rhodr)*DHs);
    Tcl     = 20+TK;                 % Chamber temperature [K]
    kd1   = 4.75e3;   % Parameter in the mass resistance equation
    kd2   = 6.051e7;  % Parameter in the mass resistance equation
    hL1   = 3.85;     % Parameter in the heat transfer coefficient
    hL2   = 0.352;    % Parameter in the heat transfer coefficient
    hL3     = 34;                   % Parameter in the heat transfer coefficient

%... Control variables
    Pc  = 20;    % Chamber pressure [Pa]
    Tcr = TK+20; % Shelf temperature [K]

%... Time data
    tf  = 26000; % final time [s]
    tl  = linspace(0,tf,200);  

%... Initial conditions (read experiments)
    T_ini       = 244*ones(25,1);                % Dried region initial temperature
    x0          = L*0.015;           % Initial front position (x0>0)
    xi          = linspace(0,L,length(T_ini)); 
    ze          = linspace(0,1,ndz)';   % Dryed region (Landau) transformed spatial coordinates [-]
    ye          = linspace(0,1,ndy)';   % Frozen region (Landau) transformed spatial coordinates [-]
%... Initial conditions (adaptation to Landau)
    xi_ndz = linspace(0,x0,ndz);
    xi_ndy = linspace(x0,L,ndy);
    Rdr0   = interp1(xi, T_ini, xi_ndz, 'pchip')';
    Rfr0   = interp1(xi, T_ini, xi_ndy, 'pchip')';
    y0     = [Rdr0; Rfr0; x0];
    clear T0 T_ini

%... FEM matrices
    [MMdr, DMdr, CMdr] = matfem (ze,'neu','neu','MM','DM','CM');
    [MMfr, DMfr, CMfr] = matfem (ye,'neu','neu','MM','DM','CM');
    inv_MMdr = MMdr\eye(size(MMdr));
    inv_MMfr = MMfr\eye(size(MMfr));
%... Spatial operators
    Lap_op_dr = -inv_MMdr*DMdr;
    Adv_op_dr = inv_MMdr*CMdr;
    Lap_op_fr = -inv_MMfr*DMfr;
    Adv_op_fr = inv_MMfr*CMfr;


%... Problem integration with ode15s
%... Call the integrator
    atol       = 1e-6;
    rtol       = 1e-6;
    optionspd  = odeset('AbsTol',atol,'RelTol',rtol);
    [tpd, ypd] = ode15s(@ode_pd, tl, y0, optionspd,...
                        Adv_op_dr, Adv_op_fr, Lap_op_dr, Lap_op_fr,...
                        inv_MMdr, inv_MMfr, ze, ye, ndz, ndy, L, alpha,...
                        Kfr, Kdr, sigma, ep, fp, Tcl, ...
                        rhodr, cpdr, rhofr, cpfr, kd1, kd2,...
                        hL1, hL2, hL3, Pc, Tcr);
                    
%... Store the data
    Rdr      = ypd(:,1:ndz)';
    Rfr      = ypd(:,1+ndz:ndz+ndy)';
    xx       = ypd(:,ndz+ndy+1)';
    
%... Transformed spatial coordinates
    ntpd   = length(tl);
    zdr    = zeros(ndz,ntpd);
    zfr    = zeros(ndy,ntpd);
    for ii = 1:length(xx)
        zdr(:,ii) = xx(ii)*ze;
        zfr(:,ii) = (L-xx(ii))*ye+xx(ii);
    end

%... Computation of the temperature (undo the landau transform)
    ndxi  = 11;
    pts   = linspace(0,L,ndxi);
    zz    = [zdr;zfr(2:end,:)];
    RR    = [Rdr;Rfr(2:end,:)];
    Tpd   = zeros(ndxi,size(zz,2));
    for jj = 1 : size(zz,2)
        Tpd(:,jj) = interp1(zz(:,jj),RR(:,jj),pts,'pchip');
    end

%... Plot the solution
    plot(tl(1:10:end)/3600,Tpd(1:2:end,1:10:end),'LineWidth',2)
    title('Product temperature')
    xlabel('Time(hours)','Fontsize',16)
    ylabel('Temperature (K) ','Fontsize',16)
