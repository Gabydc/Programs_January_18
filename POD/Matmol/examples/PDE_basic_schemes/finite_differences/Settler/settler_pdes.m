%... The MatMol Group (2016)
    function xt = settler_pdes(t, x)

%... The Matmol Group 2009

%... Set global variables
    global A
    global Cf Qf qun qov
    global D1_1  D1_2  D2_1  D2_2 D11 D12 D1_1_
    global vec_D
    global n1 n2

%... Transfer dependent variables
    C1 = x(1:n1);
    C2 = x(n1+1:n1+n2);

%... Compute the settling velocities
    nu1 = nu_takacs(C1);
    nu2 = nu_takacs(C2);

%... Compute the solid fluxes
    J1s = nu1'.*C1;
    J1q = -qov'.*C1;
    v2  = qun + nu2;
    J2  = v2'.*C2;

%... Compute the spatial derivatives
    J1z  = D1_1_*J1q + D1_1*J1s;
    J2z  = D1_2*J2;
    C1z  = D11*C1;
    C1zz = D2_1*C1;
    C2z  = D12*C2;
    C2zz = D2_2*C2;

%... Assemble the PDEs
    C1t = -J1z + vec_D(1:n1).*C1zz;
    C2t = -J2z + vec_D(n1+1:n1+n2).*C2zz;

%... Boundary conditions
    C1t(1)  = J1s(1) - vec_D(1)*C1z(1);                  %... clear water outlet
    C1t(n1) = C2(1) - C1(n1);                            %... feed level
    C2t(1)  = (  (Qf/A*Cf - (qov - nu1(n1))*C1(n1)) / (qun + nu2(1)) - C2(1)  ); %... feed level
    C2t(n2) = nu2(n2)*C2(n2) - vec_D(n1+n2)*C2z(n2);     %... sludge outlet

%... transfer temporal derivatives
    xt(1:n1)       = C1t;
    xt(n1+1:n1+n2) = C2t;
    xt             = xt';
