%... The MatMol Group (2016)
%... Numerical solution of the heat equation using series
%... expansion and legendre polynomials

    clear all
    clc
    
%... Spatial coordinates
    nd = 101;                 %... Discretization points
    zz = linspace(0,1,nd)';   

%... Parameters
    k = 0.1;

%... Time span
    tlist = [0:0.05:0.5 0.6:0.1:2 2.5:0.5:4];

%... Analytical solution
    syms z kk nn

%... Initial conditions
    x0   = 5*(z^2/2-z^4/4)+1;
    x0_a = 5*(zz.^2/2-zz.^4/4)+1;

%... Load the basis functions employed in the approximation. They have
%... been previoulsy computed to save time
    nf = 7;
    load shifted_legendre_polynomials
%... Integrate the function (numerically since I did not find the
%... analytical solution) to find the coefficients
    ifnum     = zeros(nf,nd);
    legpolnum = zeros(nd,nf);
    for ii = 1 : nf
        intfunc(ii)     = x0*phi(ii);
        legpolnum(:,ii) = subs(phi(ii),zz);
        m0(ii,1)        = eval(int(intfunc(ii),0,1));
    end

%... Computation of the projection of the spatial operator
    for ii = 1:nf
        for jj = 1:nf
            dphi1  = diff(phi(jj),1);
            dphi2  = diff(phi(ii),1);
            AA(ii,jj) = eval(int(dphi1*dphi2,0,1));
        end
    end

%... Call the IVP solver
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t,m] = ode15s(@equ_ode, tlist, m0, options, AA, k);

%... Recovery the field
    x_num = legpolnum*m';

%... Plot the figure
    mesh(tlist,zz, x_num)
    xlabel('Time')
    ylabel('z')
    zlabel('x(t,z)')
