%...  The Matmol group (2016)
     function M = mass
%...
%... set global variables
     global a p mu c;
     global z0 zL z n D1;
%...
%... Assemble mass matrix
     M = diag(ones(n,1),0) - mu*D1*D1;
%...     
