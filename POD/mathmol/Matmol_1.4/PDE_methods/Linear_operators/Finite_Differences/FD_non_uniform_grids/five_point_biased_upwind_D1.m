      function [D]=five_point_biased_upwind_D1(z,v)
%...
%...  The MatMol Group (2009)
%...
%...  function five_point_biased_upwind_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable x over a nonuniform grid z
%...  from biased-upwind five-point, fourth-order finite difference approximations
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
      m=1;      
      ns=5;
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
                 for i=1:3,
                   zd=z(i);
                   [w]=weights(zd,zs,ns,m);
                   D(i,1:ns)=w(1:ns,m+1)';
                 end;
%...
%...             interior points      
                 for i=4:n-1,
                 zs=z(i-3:i+1);
                 zd=z(i);
                 [w]=weights(zd,zs,ns,m);
                 D(i,i-3:i+1)=w(1:ns,m+1)';
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
                 for i=2:n-3,
                   zs=z(i-1:i+3);
                   zd=z(i);
                   [w]=weights(zd,zs,ns,m);
                   D(i,i-1:i+3)=w(1:ns,m+1)';
                 end;
%...
%...            boundary points      
                zs=z(n-ns+1:n);
                for i=3:-1:1,
                  zd=z(n-i+1);
                  [w]=weights(zd,zs,ns,m);
                  D(n-i+1,n-4:n)=w(1:ns,m+1)';
                end;
              end;                
