%...  The Matmol group (2016)
    function out = integrand2(z)
%... Computation of the sencond integrand on the nonlinear term for 
%... Burgers' equation. This term correspond with 
%...      y2 = phi_2*(xfem*xfem_z)

    global xx n

%... Dependent variable values with the FEM local nomenclature
    x1 = xx(1:n-1);
    x2 = xx(2:n);

%... Preallocation of memmory for the output
    out = zeros(n,length(z));

%... Computation of the integrand
    for i = 1:length(z)
        %... Basis function computation
        [phi_1, phi_2] = trialfunctionslagrlin(z(i));
        %... Finite element approximation of the dependent variable
        xfem  = phi_1*x1 + phi_2*x2;
        %... Finite element approximation of the first spatial derivative
        %... of the dependent variable 
        xfem_z = 1/2*(-x1 + x2);
        %... Final output
        out(2:n,i) = phi_2*((xfem).*(xfem_z));
    end

