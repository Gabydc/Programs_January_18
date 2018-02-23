%... The MatMol Group (2016)
    function ct = advection_pde(t,c)

%... set global variables
     global v D;
     global z0 zL n D1;

%... test the positiveness of c
     c(c<0)=0;

%... spatial derivatives
     cz = D1*c;
     czz = D1*cz;

%... temporal derivatives
     ct = -v.*cz + D*czz;
     ct(1) = 0;    
     ct(n) = -v(end).*cz(n);
