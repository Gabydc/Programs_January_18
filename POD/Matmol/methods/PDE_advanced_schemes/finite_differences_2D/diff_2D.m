%... The Matmol group (2016)
    function [xz2D] = diff_2D(n1,n2,nD,D,x2D)
%...  function diff_2D computes a partial derivative over a two-
%...  dimensional domain using either centered or biased upwind 
%...  approximations. Also, the coding is very simple because it 
%...  uses the corresponding one-dimensional functions.
%...
%...  argument list
%...
%...     n1        number of grid points for the first independent
%...               variable (input)
%...
%...     n2        number of grid points for the second independent
%...               variable (input)
%...
%...     nD        number of the independent variable for which the
%...               partial derivative is to be computed (input)
%...
%...     D         differentiation matrix corresponding to the independent 
%...               variable for which the partial derivative is to be computed 
%...               (input)
%...
%...     x2D       two-dimensional (n1,n2) array containing the dependent 
%...               variable which is to be differentiated with respect to
%...               independent variable nD (input)
%...
%...  the following one-dimensional arrays contain the dependent
%...  variable (x1D) and its partial derivative (xz1D).  in each
%...  case, one of the independent variables is constant and the
%...  other independent variable varies over its total interval.
%...
      if nD == 1,
%...
%...  ******************************************************************
%...
%...  the partial derivative is to be computed with respect to the
%...  first independent variable defined over an interval consisting
%...  of n1 grid points.  
      for j = 1:n2,
%...
%...     transfer the dependent variable in the two-dimensional array x2D
%...     to the one-dimensional array x1D so that a differentiation matrix
%...     can be used to calculate the partial derivative
         x1D = x2D(:,j);
%...
%...     compute the partial derivative using the differentiation matrix 
         xz1D = D*x1D;
%...
%...     return the partial derivative in the one-dimensional array xz1D
%...     to the two-dimensional array xz2D
         xz2D(:,j) = xz1D;
%...
%...  the partial derivative at a particular value of the second inde-
%...  pendent variable has been calculated.  Repeat the calculation for
%...  the next value of the second independent variable
      end
%...
      elseif nD == 2
%...
%...  ******************************************************************
%...
%...  the partial derivative is to be computed with respect to the
%...  second independent variable defined over an interval consisting
%...  of n2 grid points.
      for i = 1:n1,
%...
%...     transfer the dependent variable in the two-dimensional array x2D
%...     to the one-dimensional array x1D so that a differentiation matrix
%...     can be used to calculate the partial derivative
         x1D = x2D(i,:)';      
%...
%...     compute the partial derivative using the differentiation matrix
         xz1D = D*x1D;
%...
%...     return the partial derivative in the one-dimensional array xz1D
%...     to the two-dimensional array xz2D
         xz2D(i,:) = xz1D';      
%...
%...  the partial derivative at a particular value of the first inde-
%...  pendent variable has been calculated. Repeat the calculation for
%...  the next value of the first independent variable
      end
%...
%...  the partial derivative has been calculated over the entire n1 x
%...  n2 grid.  therefore return to the calling program with the partial
%...  derivative in the two-dimensional array xz2D
      end
