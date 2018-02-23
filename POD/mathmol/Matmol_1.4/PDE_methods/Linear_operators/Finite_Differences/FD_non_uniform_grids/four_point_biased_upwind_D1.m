      function [D]=four_point_biased_upwind_D1(z,v)
%...
%...  The MatMol Group (2009)
%...
%...  function four_point_biased_upwind_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable x over a nonuniform grid z
%...  from biased-upwind four-point, third-order finite difference approximations
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
      ns=4;
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
                 for i=1:2,
                   zd=z(i);
                   [w]=weights(zd,zs,ns,m);
                   D(i,1:ns)=w(1:ns,m+1)';
                 end;
%...
%...             interior points      
                 for i=3:n-1,
                 zs=z(i-2:i+1);
                 zd=z(i);
                 [w]=weights(zd,zs,ns,m);
                 D(i,i-2:i+1)=w(1:ns,m+1)';
                 end;
%...
%...             boundary point     
                 zs=z(n-ns+1:n);
                 zd=z(n);
                 [w]=weights(zd,zs,ns,m);
                 D(n,n-ns+1:n)=w(1:ns,m+1)';
              end;
%...
%...     (2)  finite difference approximation for negative v
              if v < 0
%...
%...             boundary point     
                 zs=z(1:ns);
                 zd=z(1);
                 [w]=weights(zd,zs,ns,m);
                 D(1,1:ns)=w(1:ns,m+1)';
%...
%...             interior points      
                 for i=2:n-2,
                   zs=z(i-1:i+2);
                   zd=z(i);
                   [w]=weights(zd,zs,ns,m);
                   D(i,i-1:i+2)=w(1:ns,m+1)';
                 end;
%...
%...            boundary points      
                zs=z(n-ns+1:n);
                for i=2:-1:1,
                  zd=z(n-i+1);
                  [w]=weights(zd,zs,ns,m);
                  D(n-i+1,n-3:n)=w(1:ns,m+1)';
                end;
              end;                
