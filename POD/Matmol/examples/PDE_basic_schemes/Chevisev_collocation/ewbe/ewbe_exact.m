%... The MatMol Group (2016)
    function [x] = ewbe_exact(z,t)

%... this function computes an exact solution to Burgers' equation
%...
     global a p mu delta c kappa phi0 z0;
%...
%...     
    lambda = 2/(mu*(p+2));
    c = 2*delta/( mu*(p+4)*sqrt(lambda) );
    kappa = p*sqrt(lambda)/4;
    phi0 = 2*sqrt( (p^2+3*p+2)*delta*sqrt(lambda)/(4*a*(p+4)) );
    x = ( phi0 * ( 1 - tanh(kappa*(z-c*t-z0)) - 0.5 * ( sech(kappa*(z-c*t-z0)) ).^2 ) ).^(1/p);    
