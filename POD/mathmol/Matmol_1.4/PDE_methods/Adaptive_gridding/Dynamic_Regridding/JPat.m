     function jacob_sparse = JPat

% The MatMOL group (2009)
%... % compute the sparsity of the Jacobian %
%...
%... The sparsity of the jacobian is PDEs-depending ; anyhow the sparsity
%...
%... matrix Js is a nxn square-blocks matrix ; the dimension of each block is (ne+1)x(ne+1).
%...
%... So the global dimension of Js is (ne+1)xn.
%...
%...

     
     global n ne X eq choice
%...
     nt=ne+1;
%...
%... finite differences choice depending part
%...
     for k=1:ne
         for i=1:n
             for j=1:ne
                 ibeg=max(0,i-1+eq(j,1,k));
                 if ibeg == 0
                     iend=eq(j,2,k)-eq(j,1,k);
                 else
                     iend=min(n-1,i-1+eq(j,2,k));
                 end
                 if iend == n-1
                     ibeg=iend-(eq(j,2,k)-eq(j,1,k));
                 end
                 X((i-1)*nt+k,ibeg*nt+j:nt:iend*nt+j)=eq(j,3,k);
                 X((i-1)*nt+k,(ibeg+1)*nt:nt:(iend+1)*nt)=eq(k,3,k);
             end
         end
     end
%
%... monitor choice depending part
%...
     switch choice
         
       case('der1')
         for i=1:2
             X(i*nt,1:3*nt)=1;
             X(i*nt,4*nt)=1;
             X(i*nt,5*nt)=1;
         end
         for i=3:n-2
             X(i*nt,(i-2)*nt:nt:(i+2)*nt)=1;
             for k=1:ne
                 X(i*nt,(i-2)*nt+k:nt:i*nt+k)=1;
             end
         end
         for i=n-1:n
             X(i*nt,(n-3)*nt+1:n*nt)=1;
             X(i*nt,(n-3)*nt)=1;
             X(i*nt,(n-4)*nt)=1;
         end
         
       case('der2')
         for i=1:2
             X(i*nt,1:5*nt)=1;
         end
         for i=3:n-2
             X(i*nt,(i-3)*nt+1:(i+2)*nt)=1;
         end
         for i=n-1:n
             X(i*nt,(n-5)*nt+1:n*nt)=1;
         end

       case('derfl')
         for i=1:2
             X(i*nt,1:3*nt)=1;
             X(i*nt,4*nt)=1;
             X(i*nt,5*nt)=1;
         end
         for i=3:n-2
             X(i*nt,(i-2)*nt:nt:(i+2)*nt)=1;
             for k=1:ne
                 X(i*nt,(i-2)*nt+k:nt:i*nt+k)=1;
             end
         end
         for i=n-1:n
             X(i*nt,(n-3)*nt+1:n*nt)=1;
             X(i*nt,(n-3)*nt)=1;
             X(i*nt,(n-4)*nt)=1;
         end

     end
     jacob_sparse=sparse(X);
