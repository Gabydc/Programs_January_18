%... The MatMol Group (2016)
     function Jac = jacobian_complex_diff(t,x)
%...
%... Jacobian through complex step differentiation
%...
     nx = length(x);
     Jac = zeros(nx,nx);
     h = nx*eps;
     for k = 1:nx
         x1 = x;
         x1(k) = x1(k)+h*i;
         xt = brusselator_pdes(t,x1);
         Jac(:,k)=imag(xt)/h;
     end
     Jac = sparse(Jac);