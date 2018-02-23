%... The MatMol Group (2016)
    function xt = algae_bloom_pde(t,x)

%... set global variables
    global mu r K kd;
    global z0 zL z n D1;

%... boundary conditions at z = 0
    x(1) = 0;

%... boundary conditions at z = L
    x(n) = 0;

%... spatial derivatives (stagewise differentiation)
    xz  = D1*x;
    xzz = D1*xz;

%... temporal derivatives
    xt    = mu*xzz + r*x.*(1-x/K);
    xt(1) = 0;
    xt(n) = 0;
