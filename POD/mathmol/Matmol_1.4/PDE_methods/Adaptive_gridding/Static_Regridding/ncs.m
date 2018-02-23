	function s=ncs(x,y)
%   
%  The MatMol Group (2009)
%
%  Natural cubic spline
%
%  This code is based on algorithms detailed in the book Numerical Methods Using Matlab
%  by J.H. Mathews and K.D. Fink (Prentice Hall, 1999)
%  Modifications by W.E. Schiesser, P. Saucez and A. Vande Wouwer (2001)  
%
%  input parameters:
%
%     x             : abscissa vector
%     y             : ordinate vector
%
%  output parameters:
%
%     s(k,:)        : coefficients, in descending order, of the kth cubic interpolant
%
   n=length(x)-1;
   s=zeros(n+1,4);
   m=zeros(n+1,1);   
%
%  Set up the tridiagonal system of equations
%  b = diagonal coefficients, c = offdiagonal coefficients,
%  v = right hand side vector
%     
   h=diff(x);
   d=diff(y)./h;
   b=2*(h(1:n-1)+h(2:n));
   c=h(2:n);
   v=6*diff(d);
%
%  Natural spline endpoint constraints (zero second-order derivatives)
%
%  No changes in b(1), u(1), b(n-1), u(n-1)
%   
%  Step through the forward elimination
   for k=2:n-1;
      temp=c(k-1)/b(k-1);
      b(k)=b(k)-temp*c(k-1);
      v(k)=v(k)-temp*v(k-1);
   end;   
%
%  Step through the back substitution
%
   m(n)=v(n-1)/b(n-1);
   for k=n-2:-1:1;
      m(k+1)=(v(k)-c(k)*m(k+2))/b(k);
   end;
%
%  Endpoint constraints
%
   m(1)=0;
   m(n+1)=0;
%
%  Compute the coefficients s(k,1:4) of the kth cubic spline interpolant
%
   s(1:n,1)=(m(2:n+1)-m(1:n))./(6*h(1:n)');
   s(1:n,2)=m(1:n)/2;
   s(1:n,3)=d(1:n)'-h(1:n)'.*(2*m(1:n)+m(2:n+1))/6;
   s(1:n,4)=y(1:n)';
%
%  Complete the matrix s with the values at the right end-point
%
   s(n+1,1)=s(n,1);
   s(n+1,2)=0;
   s(n+1,3)=s(n,3)+2*s(n,2)*h(n)+3*s(n,1)*h(n)^2;
   s(n+1,4)=y(n+1);




