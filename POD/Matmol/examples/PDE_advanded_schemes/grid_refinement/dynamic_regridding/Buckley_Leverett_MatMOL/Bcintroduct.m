%...  The MatMol Group (2016)
     function [u z] = Bcintroduct(t,x)
%...
%... Implement the boundary conditions in Dirichlet form and add the boundary
%... nodes
     global n ne zL zR

%... Separate dependent variables and node positions
    for jj = 1:ne,
        u(2:n+1,jj) = x(jj:ne+1:n*(ne+1));
    end
    z = [zL x(ne+1:ne+1:n*(ne+1))' zR];

%... Implement the BCs
    u(1,1)   = 1;
    u(n+2,1) = 0;
