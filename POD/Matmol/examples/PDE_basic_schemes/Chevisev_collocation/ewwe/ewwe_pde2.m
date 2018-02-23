%...  The Matmol group (2016)
     function xt = ewwe_pde2(t,x)
%...
%... set global variables
     global a p mu c;
     global z0 zL z n D1;
     global L U
%...    
%... spatial derivatives
     xz = D1*x;
%...
%... temporal derivatives
%...
     f = -a*(x.^p).*xz;
     xt = U\(L\f);
