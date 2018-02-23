%... The Matmol group (2016)
    function mass_sparse = MvPat
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%... % compute the sparsity of the mass-matrix %
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
%... The sparsity of the mass-matrix A(x,t) depends only on the number of 
%...
%... PDEs (ne),the number of interior nodes (n) and the monitor function. 
%...
%... The sparsity is described by a As (ne+1)*n square matrix which is 
%...
%... (ne+1)*(ne+1) block-pentadiagonal. Three types of blocks are present :
%...
%...              0 . . . . 0                                1 0 . . 0 1                                 0 . . . . 0
%...              0 . . . . 0                                0 1 0 . 0 1                                 0 . . . . 0
%...   As(i,i) =  . . . . . .       As(i,i+1) = As(i,i-1) =  . . . . . .        As(i,i+2) = As(i,i-2) =  . . . . . .
%...              . . . . . .                                . . . . . .                                 . . . . . .
%...              0 . . . . 0                                0 . . 0 1 1                                 0 . . . . 0
%...              1 1 1 1 1 1                                1 1 1 1 1 1                                 0 . . . 0 1
 
%... so that
%...
%...                As(1,1)  As(1,2)  As(1,3)     0        0        0        0        0      . . . . . . . .  0
%...
%...                As(2,1)  As(2,2)  As(2,3)  As(2,4)     0        0        0        0      . . . . . . . .  0
%...
%...                As(3,1)  As(3,2)  As(3,3)  As(3,4)  As(3,5)     0        0        0      . . . . . . . .  0
%...
%...                   0     As(4,2)  As(4,3)  As(4,4)  As(4,5)  As(4,6)     0        0      . . . . . . . .  0
%...
%...                   0        0     As(5,3)  As(5,4)  As(5,5)  As(5,6)  As(5,7)     0      . . . . . . . .  0
%...
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
%...    As(x,t) = 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
%... 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
%... 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
%... 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
%... 
%...                  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
%...
%...                  0   . . . . . . . . . . . . . . . . . . . . . . . . . . 0      As(n,n-2)  As(n,n-1)  As(n,n)
%...
%...
     
     global n ne
     global choice 

     % General blocs of MvPattern
     bloc0            = zeros(ne+1,ne+1);
     bloc1            = bloc0;
     bloc1(ne+1,:)    = 1;
     bloc2            = eye(ne+1);
     bloc2(:,ne+1)    = 1;
     bloc2(ne+1,:)    = 1;
     bloc3            = bloc0;
     bloc3(ne+1,ne+1) = 1;

     switch choice
         
         case('der1')
             l1=[bloc1 bloc2 bloc3 repmat(bloc0,1,n-3)];
             l2=[bloc2 bloc1 bloc2 bloc3 repmat(bloc0,1,n-4)];
             for i=1:n-4
                 l(1+(i-1)*(ne+1):i*(ne+1),1:n*(ne+1))=[repmat(bloc0,1,i-1) bloc3 bloc2 bloc1...
                     bloc2 bloc3 repmat(bloc0,1,n-4-i)];
             end
             lnm1=[repmat(bloc0,1,n-4) bloc3 bloc2 bloc1 bloc2];
             ln=[repmat(bloc0,1,n-3) bloc3 bloc2 bloc1];
             X=cat(1,l1,l2,l,lnm1,ln);
           
         case('der2')
             l1=[bloc1 bloc2 bloc1 repmat(bloc0,1,n-3)];
             l2=[bloc2 bloc1 bloc2 bloc1 repmat(bloc0,1,n-4)];
             for i=1:n-4
                 l(1+(i-1)*(ne+1):i*(ne+1),1:n*(ne+1))=[repmat(bloc0,1,i-1) bloc1 bloc2 bloc1...
                     bloc2 bloc1 repmat(bloc0,1,n-4-i)];
             end
             lnm1=[repmat(bloc0,1,n-4) bloc1 bloc2 bloc1 bloc2];
             ln=[repmat(bloc0,1,n-3) bloc1 bloc2 bloc1];
             X=cat(1,l1,l2,l,lnm1,ln);

         case('derfl')
             l1=[bloc1 bloc2 bloc3 repmat(bloc0,1,n-3)];
             l2=[bloc2 bloc1 bloc2 bloc3 repmat(bloc0,1,n-4)];
             for i=1:n-4
                 l(1+(i-1)*(ne+1):i*(ne+1),1:n*(ne+1))=[repmat(bloc0,1,i-1) bloc3 bloc2 bloc1...
                     bloc2 bloc3 repmat(bloc0,1,n-4-i)];
             end
             lnm1=[repmat(bloc0,1,n-4) bloc3 bloc2 bloc1 bloc2];
             ln=[repmat(bloc0,1,n-3) bloc3 bloc2 bloc1];
             X=cat(1,l1,l2,l,lnm1,ln);

     end
     mass_sparse=sparse(X);
     