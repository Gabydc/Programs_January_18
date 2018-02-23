%...  The Matmol group (2016)
    function out = integrand2_hermite(xi)
%...
%... Computation of the second integrand on the nonlinear term for
%... the Burgers' equation. This term correspond with 
%...      y2 = phi_2*(xfem*xfem_z)

    global xx n

%... Dependent variable values and their derivatives with the FEM
%... local nomenclature
    x1    = xx(1:2:2*n-3);
    x1_xi = xx(2:2:2*n-2);
    x2    = xx(3:2:2*n-1);
    x2_xi = xx(4:2:2*n);

%... Preallocation of memmory for the output
    out = zeros(2*n,length(xi));

    for i = 1:length(xi)
        %... Basis function computation
        [phi1, phi2, phi3, phi4]     = trialfunctionsher2n(xi(i));
        %... First derivative of the basis function computation
        [phi1p, phi2p, phi3p, phi4p] = trialfunctionsher2np(xi(i));
        %... Finite element approximation of the dependent variable
        xfem = phi1*x1 + phi2*x1_xi + phi3*x2 + phi4*x2_xi;
        %... Finite element approximation of the first spatial derivative
        %... of the dependent variable 
        xfem_z = phi1p*x1 + phi2p*x1_xi + phi3p*x2 + phi4p*x2_xi;
        out(2:2:2*n-2,i) = phi2*(xfem.*xfem_z);
    end
