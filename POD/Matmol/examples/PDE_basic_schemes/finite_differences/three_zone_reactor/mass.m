%... The MatMol Group (2016)
     function M = mass
%...
%... set global variables
     global n1 n2 n3;
%...
%... Assemble mass matrix
     M = diag([1 1 1 ones(1,3*n1-6) 0 0 0 0 0 0 1 ones(1,4*n2-8) 0 0 0 1 0 0 0 ones(1,3*n3-6) 0 0 0],0);