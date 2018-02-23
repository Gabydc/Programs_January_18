%... The MatMol Group (2016)
    function fx = dflux_dx(ne,t,x) 

%...    
    global gam

%...    
    n = length(x);
    
%...    
    fx(1:n,1,2) = 1;
    fx(1:n,2,1) = .5*(gam-3)*(x(1:n,2).^2)./(x(1:n,1).^2);
    fx(1:n,2,2) = -(gam-3)*x(1:n,2)./x(1:n,1);
    fx(1:n,2,3) = gam-1;
    fx(1:n,3,1) = (gam-1)*((x(1:n,2)./x(1:n,1)).^3) - gam*x(1:n,3).*x(1:n,2)./(x(1:n,1).^2);
    fx(1:n,3,2) = gam*x(1:n,3)./x(1:n,1)-1.5*(gam-1)*((x(1:n,2)./x(1:n,1)).^2);
    fx(1:n,3,3) = gam*x(1:n,2)./x(1:n,1);
