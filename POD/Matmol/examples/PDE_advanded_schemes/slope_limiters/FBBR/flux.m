%...  The MatMol Group (2016)
    function [f] = flux(n,ne,t,x)
%... Computes the flux

    global v
    f(1:n,1)=v*x(1:n,1);
    f(1:n,2)=0*x(1:n,2);

    end 
