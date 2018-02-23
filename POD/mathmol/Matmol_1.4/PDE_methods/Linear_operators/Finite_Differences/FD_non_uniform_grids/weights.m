      function [w] = weights(zd,zs,ns,m)
%...
%...  The MatMol Group (2009)
%...
%...  weighting coefficients for finite difference approximations, computed
%...  by an algorithm of B. Fornberg (1,2), are used in the following 
%...  approximations.
%...
%...  (1)  Fornberg, B., fast generation of weights in finite difference
%...       formulas, in recent developments in numerical methods and
%...       software for odes/daes/pdes, G. Byrne et al (eds), World
%...       Scientific, River Eedge, NJ, 1992
%...
%...  (2)  Fornberg, B., calculation of weights in finite difference
%...       formulas, Siam Review, vol. 40, no. 3, pp 685-691, September, 
%...       1999
%...
%...
%...  slight adaptations by: W.E. Schiesser, P. Saucez and A. Vande Wouwer
%...  
%...  this function computes the weights of a finite difference scheme
%...  on a nonuniform grid
%...
%...  input Parameters
%...
%...       zd              location where the derivative is to be computed
%...
%...       ns              number of points in the stencil
%...
%...       zs(ns)          stencil of the finite difference scheme
%...
%...       m               highest derivative for which weights are sought
%...
%...  output Parameter
%...
%...       w(1:ns,1:m+1)   weights at grid locations z(1:ns) for derivatives
%...                       of order 0:m, found in w(1:ns,1:m+1)
%...                    
%...
      c1 = 1.0;
      c4 = zs(1)-zd;
      for k=0:m
        for j=0:ns-1
          w(j+1,k+1) = 0.0;
      end
      end
      w(1,1) = 1.0;
      for i=1:ns-1
         mn = min(i,m);
         c2 = 1.0;
         c5 = c4;
         c4 = zs(i+1)-zd;
         for j=0:i-1
           c3 = zs(i+1)-zs(j+1);
           c2 = c2*c3;
           if (j==i-1)
             for k=mn:-1:1
               w(i+1,k+1) = c1*(k*w(i,k)-c5*w(i,k+1))/c2;
             end 
             w(i+1,1) = -c1*c5*w(i,1)/c2;
           end
           for k=mn:-1:1
             w(j+1,k+1) = (c4*w(j+1,k+1)-k*w(j+1,k))/c3;
           end
             w(j+1,1) = c4*w(j+1,1)/c3;
         end
         c1 = c2;
       end

