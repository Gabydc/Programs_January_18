     function g = rightmon(mon,dltz,mu)
% The MatMOL group (2009)
%...
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%... % compute the right vector g(n) of the moving grid equation %
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
     global n
%...
%... compute the inverse of mon and dltz
%...
     mo_inv=1./mon;
     dz_inv=1./dltz;
%...
%... compute g(n)
%...
     col=ones(n,1);
     p1=spdiags([-mu*col (1+2*mu)*col -mu*col],1:3,n,n+3);
     p2=spdiags([-mu*col (1+2*mu)*col -mu*col],0:2,n,n+3);
     res1=p1*dz_inv';
     res2=p2*dz_inv';
     dmon1=spdiags(mo_inv(2:n+1)',0,n,n);
     dmon2=spdiags(mo_inv(1:n)',0,n,n);
     g=dmon1*res1-dmon2*res2;
   