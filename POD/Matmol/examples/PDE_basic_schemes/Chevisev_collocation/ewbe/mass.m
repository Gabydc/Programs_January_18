%... The MatMol Group (2016)
     function M = mass(t,x)

%... set global variables
     global a p mu delta c kappa phi0 z0;
     global zL zR z n D1;
%...
%... Assemble mass matrix
     M = diag(ones(n,1),0) - mu*D1*D1;
     M(1,1:n) = 0;
     M(n,1:n) = 0;
