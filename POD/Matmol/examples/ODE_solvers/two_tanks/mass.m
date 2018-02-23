%... The MatMol Group (2016)
     function M = mass(t,x)
%...
%... set global variables
     global A B eps K
%...
     M = [1 0  0   0; 
          0 1  0   0; 
          0 0 eps  0;
          0 0  0  eps];