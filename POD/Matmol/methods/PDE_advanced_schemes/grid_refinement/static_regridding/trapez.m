%... The Matmol group (2016)
    function fint=trapez(z,f)
%... fint=trapez(z,f)
%...
%... trapez computes the integral of f between z(1) and z(i) using 
%... the trapezoidal rule
%...
      n=length(z);
      fint=zeros(1,n);
      for i=2:n;
          fint(i)=fint(i-1)+0.5*(f(i)+f(i-1))*(z(i)-z(i-1));
      end;
