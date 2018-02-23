%...  The Matmol group (2016)
     function [pods] = matpod(V , MM , mth , ener)
%...
%... MATPOD  Computes the POD basis functions.
%...
%... [pods] = matpod(V , MM , mth , ener)
%... computes the POD basis functions and their associated
%... eigenvalues. 
%
%...  Input variables:
%... V:    Set of simulation or experimental data. Each column
%...       corresponds with a snapshot at a given time
%... MM:   Mass matrix obtained from the FEM. The FEM
%...       discretization must coincide with the discretization for
%...       V 
%... mth:  Method used for computing the basis functions. 'd' is
%...       the direct method, 'i' is the indirect method. The
%...       indirect method results more efficient if the number of
%...       discretization points is much larger than the number of
%...       measurements. 
%... ener: Energy captured by the PODs.
%
%... Output variables:
%... pods: This is a structure where pods.phi are the basis
%...       functions and pods.lambda the eigenvalues


%... Default method and energy
     switch (nargin)
         case {2}
             mth     = 'd';
             ener    = 99.99;
         case {3}
             ener    = 99.99;
     end
%... Number of states
     ndisc = size(MM , 1);
     spods = size(V , 1);
     nstat = round(spods/ndisc);

%... Correlation matrices
     L = size(V , 2);
     for i = 1 : nstat
         U(i).V  = V((i-1)*ndisc + 1 : i*ndisc , :);
         switch (mth)
             case {'i'}
                 RD(i).V = 1/L*U(i).V'*MM*U(i).V;
             case {'d'}
                 R(i).V  = 1/L*U(i).V*U(i).V';
                 RD(i).V = R(i).V*MM;
             otherwise
                 'Non existing method'
                 break
         end
     end
     clear V i R

%... Computation of Eigenfunctions and eigenvalues
     for i = 1 : nstat
         [P(i).P , P(i).lambda_nn] = eig(full(RD(i).V));
         switch (mth)
             case {'i'}
                 P(i).phi_nn  = U(i).V*P(i).P;
             case {'d'}
                 P(i).phi_nn  = P(i).P;
             otherwise
                 'Non existing method'
                 break
         end
    
    %... Selection of the more representative POD
         tot_ener(i) = sum(diag(P(i).lambda_nn));
         for j = 1 : L
             if (ener > 100)
                 P(i).phinn      = P(i).phi_nn(: , 1:end);
                 pods(i).lambda  = P(i).lambda_nn(1:end , 1:end);
             else
                 par_ener(j) = sum( diag( P(i).lambda_nn(1:j , 1:j)));
                 cap_ener    = par_ener(j)/tot_ener(i)*100;
                 if (cap_ener >= ener)
                     P(i).phinn      = P(i).phi_nn(: , 1:j);
                     pods(i).lambda  = P(i).lambda_nn(1:j , 1:j);
                     break
                 end
             end
         end
         clear j par_ener cap_ener
     end
     clear i tot_ener par_ener ener RD U

%... Eigenfunctions normalization
     for i = 1 : nstat
         for ii = 1 : size(P(i).phinn,2)
             cc(i).par(ii)    = P(i).phinn(:,ii)'*MM*P(i).phinn(:,ii); 
             pods(i).phi(:,ii)= 1/sqrt(cc(i).par(ii))*P(i).phinn(:,ii);
         end
     end
     clear i ii cc P 