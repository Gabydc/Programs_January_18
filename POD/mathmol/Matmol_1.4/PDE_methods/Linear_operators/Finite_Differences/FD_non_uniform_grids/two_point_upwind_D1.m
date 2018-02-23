      function [D] = two_point_upwind_D1(z,v)
%...
%...  The MatMol Group (2009)
%...
%...  function two_point_upwind_D1 returns the differentiation matrix  
%...  for computing the first derivative, xz, of a variable x over a nonuniform 
%...  grid z from upwind two-point, first-order finite difference approximations
%...
%...  the following parameters are used in the code:
%...
%...       z               spatial grid
%...
%...       v               fluid velocity (positive from left to right - only the sign is used)
%...
%...       n               number of grid points
%...
%...       zs(ns)          stencil of the finite difference scheme
%...
%...       ns              number of points in the stencil
%...
%...       zd              location where the derivative is to be computed
%...
%...       m               highest derivative for which weights are sought
%...
        m=1;      
        ns=2;
%...
%...  sparse discretization matrix      
        n = length(z);
        D = sparse(n,n);
%...
%...     (1)  finite difference approximation for positive v     
              if v > 0
%...
%...             boundary point     
                   zs=z(1:ns);
                   zd=z(1);
                   [w]=weights(zd,zs,ns,m);
                   D(1,1:2)=w(1:ns,m+1)';
%...
%...             interior points      
                   for i=2:n,
                     zs=z(i-1:i);
                     zd=z(i);
                     [w]=weights(zd,zs,ns,m);
                     D(i,i-1:i)=w(1:ns,m+1)';
                   end;
              end;
%...
%...     (2)  finite difference approximation for negative v
              if v < 0
%...
%...             interior points      
                   for i=1:n-1,
                     zs=z(i:i+1);
                     zd=z(i);
                     [w]=weights(zd,zs,ns,m);
                     D(i,i:i+1)=w(1:ns,m+1)';
                   end;
%...
%...            boundary point      
                  zs=z(n-1:n);
                  zd=z(n);
                  [w]=weights(zd,zs,ns,m);
                  D(n,n-1:n)=w(1:ns,m+1)';
%...		  
              end;                
