function [phi , lambda] = matlsd(MM, DM, BM, varargin)

%==========================================================================
% The MatMol Group (2009)
%
% MATLSD  Uses the FEM matrices to compute the eigenfunctions and 
% eigenvalues of the laplacian operator
%
%    [phi , lambda] = matlsd( MM, DM, BM, k, v)
%
%    Input variables:
%    MM:   Mass matrix of the FEM
%    DM:   Diffusion matrix of the FEM
%    BM:   Boundary (homogeneous part) matrix of the FEM
%    k:    Diffusivity
%    v:    Velocity
%
%    Output variables:
%    phi:    Eigenvectors of the laplacian operator
%    lambda: Eigenvalues of the laplacian operator
%==========================================================================

% Number of input parameters
nip = length(varargin);
switch nip
    case {1}
        k   = varargin{1};
        v   = 1;
    case {2}
        k   = varargin{1};
        v   = varargin{2};
    otherwise
        'warning: two much input arguments'
end

% Computation of the inverse of the mass matrix
inv_MM	= inv(MM);

% Discretized (FEM) version of the laplacian operator
P   = inv_MM*(k*DM + v*BM);
clear DM BM

% Eigenvalue Problem
N_eig               = size(P,1);
Nm_iter             = 20000;
opts.maxit          = Nm_iter;
[phi_nnd,lambdad]   = eig(full(P));
diag_lambdad        = diag(lambdad);
clear P lambdad N_max N_eig

% The eigenvalues (an eigenfunctions) are sorted
[diag_lambda , i]   = sort(diag_lambdad);
lambda              = diag(diag_lambda);
phi_nn              = phi_nnd(: , i);
clear phi_nnd diag_lambdad diag_lambda

% Eigenfunctions normalization
for k = 1 : size(phi_nn,2)
    cc(k)       = phi_nn(:,k)'*MM*phi_nn(:,k);
    phi(:,k)    = 1/sqrt(cc(k))*phi_nn(:,k);
end
clear phi_nn

if (lambda(1,1) > lambda(end,end))
    'Warning: The last eigenvalue is larger than the first one'
end
