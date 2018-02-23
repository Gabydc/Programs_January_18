%... The MatMol Group (2016)
    function N = logistic_exact(t)
%...
%... Set global variables
     global a b K N0
%...
     N = K./(1+(K/N0-1)*exp(-a*t));
