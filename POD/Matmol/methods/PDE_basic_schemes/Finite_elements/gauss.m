%...  The Matmol group (2016)
     function I = gauss(f,n)
%...    
%... This code computes the integral of a given function f over the
%... interval x = [-1,1] using the Gauss-Legendre quadrature
%... Input parameters:
%... f:   matlab function with the expression of the function to be
%...      integrated
%... n:   Number of points to be employed in the formula
%...
%... Output parameter:
%... I:   Numerical value of the integral

%... Computation of the abcissa
    beta      = 0.5./sqrt(1-(2*(1:n)).^(-2));
    T         = diag(beta,1) + diag(beta,-1);
    [V,D]     = eig(T);
    xquad     = diag(D);
    [xquad,i] = sort(xquad);

%... Computation of the weights
    w = 2*V(1,i).^2;

%... Integral value
    I = w*feval(f,xquad);