      function [D]=seven_point_centered_D3(z)
%...
%...  The MatMol Group (2009)
%...
%...  function seven_point_centered_D3 returns the differentiation matrix 
%...  for computing the third derivative, xzzz, of a variable x over a nonuniform 
%...  grid z from seven-point, sixth-order finite difference approximations
%...
%...  the following parameters are used in the code:
%...
%...       z               spatial grid
%...
%...       n               number of grid points
%...
%...       zs(n)           stencil of the finite difference scheme
%...
%...       ns              number of points in the stencil
%...
%...       zd              location where the derivative is to be computed
%...
%...       m               highest derivative for which weights are sought
%...
      m=3;      
      ns=7;
%...
%...  sparse discretization matrix      
      n = length(z);
      D = sparse(n,n);
%...
%...  boundary points      
      zs=z(1:ns);
      for i=1:3,
        zd=z(i);
        [w]=weights(zd,zs,ns,m);
        D(i,1:ns)=w(1:ns,m+1)';
      end;
%...
%...  interior points      
      for i=4:n-3,
        zs=z(i-3:i+3);
        zd=z(i);
        [w]=weights(zd,zs,ns,m);
        D(i,i-3:i+3)=w(1:ns,m+1)';
      end;
%...
%...  boundary points      
      zs=z(n-ns+1:n);
      for i=3:-1:1,
        zd=z(n-i+1);
        [w]=weights(zd,zs,ns,m);
        D(n-i+1,n-6:n)=w(1:ns,m+1)';
      end;
