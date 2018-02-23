%... The MatMol Group (2016)
     function Jac = jacobian_complex_diff2(t,x)
%...
%... Jacobian through complex step differentiation

     n = length(x);
     Jac = zeros(n,n);
     h = n*eps;
     for k = 1:n
         x1 = x;
         x1(k) = x1(k)+h*i;
         xt = bioreactor_pdes2(t,x1);
         Jac(:,k)=imag(xt)/h;
     end
     Jac = sparse(Jac);