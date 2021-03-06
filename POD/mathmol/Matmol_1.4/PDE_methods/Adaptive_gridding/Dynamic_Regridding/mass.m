     function M = mass(t,x)
% The MatMOL group (2009)
%... % compute the left matrix of the global ODEs system %
%...
%... The following code implements the left matrix A(x,t) of the system
%...
%...            A(x,t)xdot = b(x,t)
%...
%... with
%...          1   2         ne      1   2         ne            1         ne
%...    x = [u   u   ...   u   z   u   u   ...   u   z   ...   u   ...   u   z ]'
%...          1   1         1   1   2   2         2   2         n         n   n
%...
%... and where the meshpoints 1, 2, ... n are the interiors nodes. Nevertheless, A(x,t) uses also
%...
%... the boundary points and values z , z   , u  , u   .
%...                                 0   n+1   0    n+1
%...
%... The code   a) calls the boundary conditions of the problem
%...
%...            b) computes the monitor function
%...
%...
%... Following elements of A(x,t) are then implemented :
%...
%...            a) elements from the moving grid equation :
%...
%...                         mu
%...    B(i,i-2) = - --------------------
%...                 M(i-1)(dltz(i-2)**2)
%...
%...
%...                        mu                1 + 2mu                    mu
%...    B(i,i-1) = + ------------------ + -------------------- + --------------------
%...                 M(i)(dltz(i-1)**2)   M(i-1)(dltz(i-1)**2)   M(i-1)(dltz(i-2)**2)
%...
%...    
%...                        mu                1 + 2mu            1 + 2mu                    mu
%...    B(i,i)   = - ------------------ - ---------------- - -------------------- - -----------------
%...                 M(i)(dltz(i-1)**2)   M(i)(dltz(i)**2)   M(i-1)(dltz(i-1)**2)   M(i-1)(dltz(i)**2)
%...
%...
%...                        mu                1 + 2mu                mu
%...    B(i,i+1) = + ------------------ + ---------------- + ------------------
%...                 M(i)(dltz(i+1)**2)   M(i)(dltz(i)**2)   M(i-1)(dltz(i)**2)
%...
%...
%...                        mu
%...    B(i,i+2) = - ------------------
%...                 M(i)(dltz(i+1)**2)
%...
%...
%...            b) elements from the PDEs :
%...
%...              j      j
%...             u    - u
%...              i+1    i-1
%...    D(i,j) = -----------    i = 1, ..., n           j = 1, ..., ne
%...             z    - z
%...              i+1    i-1
%...
%...
%... Structure of A(x,t) : A(x,t) is (ne+1)x(ne+1) block - pentadiagonal with
%...
%...
%...              1  0  0 ... 0  -D(1,1)                        0 ... 0  0
%...              0  1  0 ... 0  -D(1,2)                        0 ... 0  0
%...    A(i,i) =  0  0  1 ... 0  -D(1,3)              A(i,j) =  0 ... 0  0
%...                    .....                                    ......
%...              0  0  0 ... 1  -D(1,ne)                       0 ... 0  0  
%...              0  0  0 ... 0  tau*B(i,i)                     0 ... 0  tau*B(i,j)
%...
%...
%... so that
%...
%...                A(1,1)  A(1,2)  A(1,3)    0       0       0       0       0      . . . . . . . .  0
%...
%...                A(2,1)  A(2,2)  A(2,3)  A(2,4)    0       0       0       0      . . . . . . . .  0
%...
%...                A(3,1)  A(3,2)  A(3,3)  A(3,4)  A(3,5)    0       0       0      . . . . . . . .  0
%...
%...                  0     A(4,2)  A(4,3)  A(4,4)  A(4,5)  A(4,6)    0       0      . . . . . . . .  0
%...
%...                  0       0     A(5,3)  A(5,4)  A(5,5)  A(5,6)  A(5,7)    0      . . . . . . . .  0
%...
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
%...    A(x,t) = 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .                                                                 
%... 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .                                                                 
%... 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .                                                                 
%... 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .                                                                 
%... 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
%...
%...                  0       . . . . . . . . . . . . . . . . . . . . .       0   A(n,n-2) A(n,n-1) A(n,n)
%...
%...
%...
%... N.B. alpha, mu and tau are parameters of the moving grid procedure.
%...
%...
%...
     global zL zR n ne
     global alpha kappa mu tau
     global u z dltz mon dltu
     global mo_inv dz_inv
     global choice
%...     
%     t
%...
%... separate dependent variables and node positions and implement the BCs
%...
     [u z] = Bcintroduct(t,x);
%...
%... Compute the monitoring function mon(i)
%...
     mon = monitor(u,z,choice,t);
%...
%... compute the inverse of mon and dltz
%...
     mo_inv=1./mon;
     dz_inv=1./dltz;
%...
%... compute the left members of the pdes d(i,j)
%...
     d=(u(3:n+2,:)-u(1:n,:))./repmat((z(3:n+2)-z(1:n))',1,ne);
%...         
%... assemble the mass matrix
%...
     v1=mo_inv.*(dz_inv(1:n+1).^2);
     v2=mo_inv.*(dz_inv(2:n+2).^2);
     v3=mo_inv.*(dz_inv(3:n+3).^2);
     diagb1=-mu*v1(3:n);
     diagb2= mu*(v1(3:n+1)+v1(2:n))+(1+2*mu)*v2(2:n);
     diagb3=-mu*(v1(2:n+1)+v3(1:n))-(1+2*mu)*(v2(2:n+1)+v2(1:n));
     diagb4= mu*(v3(1:n-1)+v3(2:n))+(1+2*mu)*v2(2:n);
     diagb5=-mu*v3(2:n-1);

     diagm1=reshape([zeros(ne,n-2);tau*diagb1],(n-2)*(ne+1),1);
     diagm2=reshape([zeros(ne,n-1);tau*diagb2],(n-1)*(ne+1),1);
     diagm3=reshape([ones(ne,n);tau*diagb3],n*(ne+1),1);
     diagm4=reshape([zeros(ne,n-1);tau*diagb4],(n-1)*(ne+1),1);
     diagm5=reshape([zeros(ne,n-2);tau*diagb5],(n-2)*(ne+1),1);
     
     d1=[diagm1' zeros(1,2*(ne+1))]';
     d2=[diagm2' zeros(1,1*(ne+1))]';
     d3=[diagm3']';
     d4=[zeros(1,1*(ne+1)) diagm4']';
     d5=[zeros(1,2*(ne+1)) diagm5']';
     for i=1:ne,
         diaginter(1:n*(ne+1),i)=...
             reshape([zeros(ne-i,n);-d(:,ne+1-i)';zeros(i,n)],n*(ne+1),1);
         dinter(1:n*(ne+1),i)=[zeros(1,i) diaginter(1:n*(ne+1)-i,i)']';
     end

     M=spdiags([d1 d2 d3 dinter d4 d5],...
         [-2*(ne+1) -(ne+1) 0 (1:ne) ne+1 2*(ne+1)],n*(ne+1),n*(ne+1));

%...
