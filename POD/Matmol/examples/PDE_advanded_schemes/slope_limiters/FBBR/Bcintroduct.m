%...  The MatMol Group (2016)
    function [u,z] = Bcintroduct(t,x)

%... Boundary conditions (they must be introduced in Dirichlet form)

%... Global variables
    global n ne zL zR
    global D1 v Sin

%... Separate dependent variables and node positions
    for jj = 1:ne
        u(2:n+1 , jj) = x(jj:ne+1:n*(ne+1));
    end
    z = [zL x(ne+1:ne+1:n*(ne+1))' zR];

%... Implement the BCs
%... Computation of the derivative matrix
    Dz       = matfd (z, 'non_uni', '1st', 5, 'centered');

%... Boundary conditions for S
    u(1,1)   = BC_dir('rob',u(:,1),Sin,Dz,'begin',v,D1);   %... First point
    u(n+2,1) = BC_dir('neu',u(:,1),0,Dz,'end',0,0);        %... Last point

%... Boundary conditions for X. 
%... Note: The dynamic equation for X is a PDE and not an ODE thus its
%...       resolution does not need the definition of boundary conditions.
%...       Nevertheless, the Dynamic Grid procedure requires the specification
%...       of Boundary conditions. A type of boundary conditions which work
%...       well in the problems we have tested is Neumann non-homogenous:
%
%...             Xz(zL,t) = Xz(zL+dzL,t)     in the first point
%...             Xz(zR,t) = Xz(zR-dzR,t)     in the last point
%
%...       where dzL = z(2)-z(1); and dzR = z(n+2)-z(n+1);
%
%... Computation of the derivative matrix
    Dzr      = matfd(z(2:end-1), 'non_uni', '1st', 5, 'centered');
    Xi1      = Dzr(1,:)*u(2:n+1,2);
    Xi2      = Dzr(end,:)*u(2:n+1,2);
    u(1,2)   = BC_dir('neu',u(:,2),Xi1,Dz,'begin',0,0);  %... First point
    u(n+2,2) = BC_dir('neu',u(:,2),Xi2,Dz,'end',0,0);    %... Last point
    end