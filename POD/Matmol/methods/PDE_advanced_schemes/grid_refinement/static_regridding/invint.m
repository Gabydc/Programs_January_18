%... The Matmol group (2016)
    function [z_new,nz_new]=invint(z,f,delta)
%...  this subroutine performs the inverse interpolation of (z_i,f_i) 
%...  to (znew_j,j*delta) using linear interpolation
%...
%...  delta should be selected, so that (f(n)-f(1))/delta is an integer <= nzmax
%...
      nz=length(z);
%...
      nz_new=1;
      z_new(nz_new)=z(1);
      fp=delta+f(1);
%...      
      for i=2:nz;
          while and((fp>f(i-1)),(fp<=f(i)))
                nz_new=nz_new+1;
	            z_new(nz_new)=z(i-1)+(fp-f(i-1))*(z(i)-z(i-1))/(f(i)-f(i-1));
                fp=fp+delta;
          end
      end;       
%...
%...  if there exists a small gap between the last node and the right boundary,
%...  then add one more node
%...
      if z_new(nz_new)<(z(nz)-1e-8*(z(nz)-z(nz-1)))
         nz_new=nz_new+1;
      end
%...
%...  the last node is the right boundary
%...
      z_new(nz_new)=z(nz);

