
function [U,S,V,SV]=POD_basis(X,varargin)
% Compute the SVD of a matrix X, the 'Covariance' option computes the
% eigenvalues of the covariance matrix, i.e. X = X-Xm;
% The option 'X_norm', normalizes the vectors of X to the largest, i.e. 
% X_i = X_i/norm(max(X_i))
opt = struct('Covariance', false, 'X_norm', false, 'U_norm', false, 'U_rank', false);
opt   = merge_options(opt, varargin{:});
Covariance = opt.Covariance;
X_norm = opt.X_norm;
U_norm = opt.U_norm;
U_rank = opt.U_rank;
%%
% X = [1  0  3
%      0  2  3
%      0  0  3];
%  X are the snapshots
[lx,ly] = size(X);
if Covariance
    % Compute mean of the snapshots
    Xm = mean(X,2);
    for i=1:ly
        X(:,i)=X(:,i)-Xm;
    end
end

if X_norm
Mx = max(X');
nmx =norm(Mx);
% Normalize the snapshots
for i=1:ly
    X(:,i)=X(:,i)/nmx; 
end
end
%%

% Compute the snapshots correlation matrix (R = 1/m X * X^T), to do this, we
% compute first the eigenvalues of D = 1/m X^T * X
%D = (1/ly)*(X'*X);
D = X' * X;
%Compute the eigenvalues (L) and eigenvectrors (V) of D
%D = sparse(D);
[V,L] = eig(D);
% Allocate S
S = zeros(ly,ly);
 
% The vector U can be computed from V and S as U = X * V *(L^T) ^ (1/2) 
% Compute *(L^T) ^ (1/2) 
g = diag(L);
g = sqrt(abs(g));
S(1:ly,1:ly) = diag(g);
for i = 1 : ly
    if abs(S(i,i))<5e-15
        S(i,i) = 0;
    end
end
% S = sparse(S);
% V = sparse(V);
% X = sparse(X);
% Compute U
U = X * V * S';
SV = S*V';
%U = U (:,1:ly);
%S = diag(S);

if U_norm
    for i=1:ly
        U(:,i)=U(:,i)/norm(U(:,i));
    end
end
if U_rank
    ru = rank(full(U));
    disp(['The rank of U is: ' num2str(ru)])
end


