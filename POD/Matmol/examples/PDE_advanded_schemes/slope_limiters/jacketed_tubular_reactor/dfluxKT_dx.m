%... The MatMol Group (2016)
    function [fx] = dfluxKT_dx(ne,t,x)
%...
    global n
%...    
    fx(:,1,1) = ones(n,1);
    fx(:,1,2) = zeros(n,1);
    fx(:,2,1) = zeros(n,1);
    fx(:,2,2) = ones(n,1);