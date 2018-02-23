     function mon = monitor(u,z,choice,t)
% The MatMOL group (2009)
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%... % compute the monitor function %
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
%...
     global ne n dltz
     global alpha
%...
%... compute dltz(i) = z(i)-z(i-1)
%...
     dltz(2:n+2)=diff(z);
     dltz(1)=dltz(2);
     dltz(n+3)=dltz(n+2);
     
     switch choice
         
         case('der1')
%...
%...     compute dltu(i,j)=u(i,j)-u(i-1,j)           
%...
             dltu(2:n+2,:)=diff(u);
%...
%...     compute monitor function :
%...
%...                              1    ne (u(i+1,j)-u(i,j))^2
%...        m(i) = sqrt[ alpha + ---  sum ------------------- ]
%...                              ne  j=1   (z(i+1)-z(i))^2
%...
%...
             quot(1:n+1,:)=dltu(2:n+2,:)./repmat(dltz(2:n+2)',1,ne);
             if ne==1
                 mon=sqrt(alpha+(quot.*quot)'/ne);
             else
                 mon=sqrt(alpha+sum((quot.*quot)')/ne);
             end

         case('der2')
%...
%...     compute d2u/dz2 by 3_pt_centered scheme           
%...
             D2 = three_point_centered_D2(z);
             uzz = D2*u;
             uzze = (uzz(1:n+1,:)+uzz(2:n+2,:))/2;
%...
%...     compute monitor function
%...
%...                              1    ne     uzz(i+1,j)+uzz(i,j)
%...        m(i) = sqrt[ alpha + ---  sum abs(-------------------)]
%...                              ne  j=1           2  
%...
%...
             if ne==1
                 mon=sqrt(alpha+(abs(uzze))'/ne);
             else
                 mon=sqrt(alpha+sum((abs(uzze))')/ne);
             end

           
         case('derfl')
%...
%...     compute dltfl[u(i,j)]=fl[u(i+1,j)]-fl[u(i,j)]           
%...
             dltu(2:n+2,:)=diff(flux(n,ne,t,u));
%...
%...     compute monitor function :
%...
%...                              1    ne (fl[u(i+1,j)]-fl[u(i,j)])^2
%...        m(i) = sqrt[ alpha + ---  sum --------------------------- ]
%...                              ne  j=1      (z(i+1)-z(i))^2
%...
%...
             quot(1:n+1,:)=dltu(2:n+2,:)./repmat(dltz(2:n+2)',1,ne);
             if ne==1
                 mon=sqrt(alpha+(quot.*quot)'/ne);
             else
                 mon=sqrt(alpha+sum((quot.*quot)')/ne);
             end
     end