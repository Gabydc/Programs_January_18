%... The MatMol Group (2016)
     function dfdt = ft(t,x)
%...
%... Time derivative of RHS (in case of a nonautonomous problem)
     dfdt = zeros(size(x));