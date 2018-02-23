%... The MatMol Group (2016)
    function xtt = wave2_pde(t,x)

%... set global variables
     global z0 zL z n D1;
%...
%... spatial derivatives
%...    
     xz = D1*x;
     xzz = D1*xz;
%...
%... temporal derivatives
%...
     xtt = xzz;
     xtt(1) = 0;         
     xtt(n) = 0;
