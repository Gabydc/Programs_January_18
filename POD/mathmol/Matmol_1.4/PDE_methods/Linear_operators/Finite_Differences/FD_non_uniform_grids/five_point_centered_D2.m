      function [D]=five_point_centered_D2(z)
%...
%...  The MatMol Group (2009)
%...
%...  function five_point_centered_D2 returns the differentiation matrix for  
%...  computing the second derivative, xzz, of a variable x over a nonuniform grid z
%...  from five-point, fourth-order finite difference approximations
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
%...
      m=2;      
      ns=5;
%...
%...  sparse discretization matrix      
      n = length(z);
      D = sparse(n,n);
%...
%...  boundary points      
      zs=z(1:ns);
      for i=1:2,
        zd=z(i);
        [w]=weights(zd,zs,ns,m);
        D(i,1:ns)=w(1:ns,m+1)';
      end;
%...
%...  interior points      
      for i=3:n-2,
        zs=z(i-2:i+2);
        zd=z(i);
        [w]=weights(zd,zs,ns,m);
        D(i,i-2:i+2)=w(1:ns,m+1)';
      end;
%...
%...  boundary points      
      zs=z(n-ns+1:n);
      for i=2:-1:1,
        zd=z(n-i+1);
        [w]=weights(zd,zs,ns,m);
        D(n-i+1,n-4:n)=w(1:ns,m+1)';
      end;
      
      
