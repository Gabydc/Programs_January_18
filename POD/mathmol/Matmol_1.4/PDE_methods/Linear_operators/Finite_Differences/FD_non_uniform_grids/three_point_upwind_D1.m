      function [D]=three_point_upwind_D1(z,v)
%...
%...  The MatMol Group (2009)
%...
%...  function three_point_upwind_D1 returns the differentiation matrix  
%...  for computing the first derivative, xz, of a variable x over a nonuniform 
%...  grid z from upwind three-point, second-order finite difference approximations
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
        ns=3;
%...
%...  sparse discretization matrix      
        n = length(z);
        D = sparse(n,n);
%...
%...     (1)  finite difference approximation for positive v     
              if v > 0
%...
%...             boundary points     
                   zs=z(1:ns);
                   for i=1:ns-1,
                     zd=z(i);
                     [w]=weights(zd,zs,ns,m);
                     D(i,1:ns)=w(1:ns,m+1)';
                   end;  
%...
%...             interior points      
                  for i=ns:n,
                     zs=z(i-(ns-1):i);
                     zd=z(i);
                     [w]=weights(zd,zs,ns,m);
                     D(i,i-(ns-1):i)=w(1:ns,m+1)';
                  end;
%...		  
              end;
%...
%...     (2)  finite difference approximation for negative v
              if v < 0
%...
%...             interior points      
                  for i=1:n-(ns-1),
                     zs=z(i:i+(ns-1));
                     zd=z(i);
                     [w]=weights(zd,zs,ns,m);
                     D(i,i:i+(ns-1))=w(1:ns,m+1)';
                  end;
%...
%...            boundary points      
                  zs=z(n-ns+1:n);
                  for i=(ns-1):-1:1,
                    zd=z(n-i+1);
                    [w]=weights(zd,zs,ns,m);
                    D(n-i+1,n-(ns-1):n)=w(1:ns,m+1)';
                  end;
%...		  
              end;                
