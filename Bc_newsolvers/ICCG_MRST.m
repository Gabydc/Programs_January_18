% This program solves the linear system Ax = b with deflation and split
% preconditioner. The error, true residual and approximation residual are
% plotted.
%
% programmer: Gabriela Diaz
% e-mail    : diazcortesgb@gmail.com
% date      : 09-11-2017
function [result,flag,res,its,resvec,resulte] = ICCG_MRST(A,b,tol,maxit,M1,M2,x0,varargin)
%warning on backtrace
%warning off verbose
% size(x0)
% size(b)
% size(A)
resulte = [];
opt = struct('Error', 1e-6, ...
    'opts', {{false, false, false, false, false, false,false,false}}, ...
    'wells', []);
% 'Residual', false, ...
% 'x_true', false, ...
% 'Convergence' , false, ...
% 'Amatrix_eigs', false, ...
% 'MAmatrix_eigs', false, ...
% 'Iter_m' , false, ...
% 'A_cn' , false, ...
% 'dir', []

opt   = merge_options(opt, varargin{:});
Error = opt.Error;
opts  = opt.opts;
W     = opt.wells;

% Display message of options
dopts      = opts{1}(1);
% Compute condition number of the matrix A
A_cn       = opts{1}(2);
% Compute true solution
x_true     = opts{1}(3);
% Checks the convergence of the method, in the case of DICCG the residual
% can increas, in that case, the solution will be the solution with
% minimal residual
Convergence = opts{1}(4);
% Save the variables to plot residual
Residual    = opts{1}(5);
% Display the number of iterations
Iter_m      = opts{1}(6);
% Compute eigenvalues of matrix A
Amatrix_eigs = opts{1}(7);
% Compute eigenvalues of matrix M^{-1}A
MAmatrix_eigs = opts{1}(8);

[n,m] = size(A);
nw = numel(W);
if nw > 0
    na = n - nw;
else
    na = n;
end
if (A_cn{1})
    SZ_A =size(A);
    C_A = cond(A);
    if(dopts{1})
        display(dopts)
        disp(['The size of A is: ' num2str(SZ_A)])
        disp(['The condition number of A is: ' num2str(C_A)])
    end
end
A(na+1 : n,:);
b(na+1 : n);
A(:,na+1 : n);
if (m ~= n)
    warning(['NonSquareMatrix'])
    error(message());
end
if ~isequal(size(b),[n,1])
    warning('RSHsizeMatchCoeffMatrix', m);
end
if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning('tooSmallTolerance');
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(['tooBigTolerance']);
    warned = 1;
    tol = 1-eps;
end
if (nargin < 5) || isempty(maxit)
    maxit = min(n,20);
end

if (nw > 0)
    for i = 1 : nw
        i;
        n-nw+i;
        bw = b(n-nw+i,1);
        Aw = A(n-nw+i,n-nw+i);
        pw(n-nw+i,1) = bw / Aw;
    end
    
    A1(1:n-nw,1:n-nw)=A(1:n-nw,1:n-nw);
    M11(1:n-nw,1:n-nw)=M1(1:n-nw,1:n-nw);
    M22(1:n-nw,1:n-nw)=M2(1:n-nw,1:n-nw);
    b1(1:n-nw,1)=b(1:n-nw,1);
    x01(1:n-nw,1)=x0(1:n-nw,1);
    A = A1;
    M1 = M11;
    M2 = M22;
    b = b1;
    x0 = x01;
    clear A1 Z1 M11 M22 b1 x01
end


if(x_true{1})
    xtrue = A\b;
    normxtrue = norm(xtrue);
end



i = 1;
x = x0;
Mb = M1 \ b;
r = b - A * x;
r = M1 \ r;
res(i,1) = norm(r);
%p = r;
p = M2 \ r;
nmb = norm(Mb);
nb = norm(b);
tol = tol*norm(nmb);

if(x_true{1} ) 
    terr(i,1)   = norm(xtrue-x);
    MAx = M1 \ (A * x);
    tresm(i,1) = norm(Mb-MAx);
    tres(i,i)  = norm(b-A*x);
end


if(res(i) < tol)
    disp(['ICCG only needs one iteration, initial residual is, r_1 = ||P*M^{-1}(b-A*xk)||_2' num2str(res(i))])
    if(x_true{1})
        disp(['True residual is, r_1 = ||b-A*xk||_2: ' num2str(tresm(i,1))])
    end
end

while  (i < maxit) && (res(i) > tol)

    i = i+1;
    w = A * p;
    alpha = (r' * r) / (p' * w);
    x = x + alpha * p;
    y = M1 \ w;
    r = r - alpha * y;
    beta = (r' * r)/(res(i-1)^2);
    %z = r;
    z = M2 \ r;
    p = z + beta * p;
    res(i,1) = norm(r);

    if(x_true{1} )
        xk = x;
        terr(i,1) = norm(xtrue-xk);
        MAxk = M1 \ (A * xk);
        tresm(i,1) = norm(Mb-MAxk);
        tres(i,1)  = norm(b-A*xk);
    end

    
    if(Convergence{1})
        E = Error * nmb;

        if ( res(i) < E )
            % If the residual increases, the approximation will be the previous
            % solution
            xacc = x;
            % If the residual increases, the approximation will be the previous
            % solution
            flagr = 1;
            rmin = res(i-1);
            %         plot(i,res(i),'*')
            %             hold on
            if (res(i) < rmin)
                rmin = res(i);
                flagr = 1;
            else if (abs(rmin -res(i)) > E)
                    
                    rmin = res(i-1);
                    flagr = 0;
                end
            end
            if flagr == 0
                warning(['Maximum accuracy is : ' num2str(res(i))])
                x =  xacc;
                break
            end
            
        end
    end
end



if (Iter_m{1} ) 
    disp(['Maximum accuracy is : ' num2str(res(i) / nmb)])
    disp(['Number of iterations is: ' num2str(i)])
end

xk = x;
tr = norm(b-A*xk)/norm(b);
Mr = M1 \ (b-A*xk);
ptr = norm(Mr)/norm(Mb);
result = xk;
flag= 0;
its =i;
resvec = res;

if (nw > 0)
    for i = 1 : nw
        result(n-nw+i,1) = pw(n-nw+i,1);
    end
end


%% Eigenvalues eigenvectors
if(Amatrix_eigs{1} )
    [VA,DA] = eigs(A,n);
    CA = cond(A,2);
end
if(MAmatrix_eigs{1} )
    IM = inv(M1);
    [VMA,DMA] = eigs(IM*A*IM',n);
    CMA = cond(IM*A*IM',2);
end
%% Save values for plots
if(A_cn{1})
    resulte.SZ_A =  SZ_A;
    resulte.C_A  =  C_A;
end
if(x_true{1})
    resulte.xtrue     = xtrue;
    resulte.terr      = terr;
    resulte.tres      = tres;
    resulte.tresm     = tresm;
end
if (Residual{1})
    resulte.res = res;
    resulte.b      = b;
    resulte.Mb     = Mb;
    resulte.tr     = tr;
    resulte.ptr    = ptr;
end
if(Amatrix_eigs{1} )
    resulte.VA = VA;
    resulte.DA = diag(DA);
    resulte.CA = CA;
end
if(MAmatrix_eigs{1} )
    resulte.IM  = IM;
    resulte.VMA = VMA;
    resulte.DMA = diag(DMA);
    resulte.CMA = CMA;
end



end

