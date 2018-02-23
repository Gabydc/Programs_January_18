      function [D]=three_point_centered_D2(z)
%...
%...  The MatMol Group (2009)
%...
%...  function three_point_centered_D2 returns the differentiation matrix for  
%...  computing the second derivative, xzz, of a variable x over a nonuniform grid z
%...  from three-point, second-order finite difference approximations
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
      m=2;      
      ns=3;
%...
%...  sparse discretization matrix      
      n = length(z);
      D = sparse(n,n);
%...
%...  boundary point     
      zs=z(1:ns);
      zd=z(1);
      [w]=weights(zd,zs,ns,m);
      D(1,1:ns)=w(1:ns,m+1)';
%...
%...  interior points      
      for i=2:n-1,
        zs=z(i-1:i+1);
        zd=z(i);
        [w]=weights(zd,zs,ns,m);
        D(i,i-1:i+1)=w(1:ns,m+1)';
      end;
%...
%...  boundary point      
      zs=z(n-ns+1:n);
      zd=z(n);
      [w]=weights(zd,zs,ns,m);
      D(n,n-2:n)=w(1:ns,m+1)';
      
      
