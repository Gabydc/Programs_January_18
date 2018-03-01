% This program solves the linear system Ax = b with deflation and split
% preconditioner. The error, true residual and approximation residual are
% plotted.
%
% programmer: Gabriela Diaz
% e-mail    : diazcortesgb@gmail.com
% date      : 09-11-2017
function [result,flag,res,its,resvec,resulte] = DICCG_MRST(A,b,Z,tol,maxit,M1,M2,x0,varargin)
%warning on backtrace
%warning off verbose
resulte = [];
opt = struct('Error', 1e-2, ...
'opts', {{false, false, false, false, false, false,false,false,false,false}}, ...
'wells', []);

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
  % Compute eigenvalues of matrix PM^{-1}A
 PMAmatrix_eigs = opts{1}(9);
   % Compute eigenvalues and condition number of matrix E
 E_cn = opts{1}(10);

nw = numel(W);

[n,m] = size(A);
if (m ~= n)
    warning(['NonSquareMatrix'])
    error(message());
end
if ~isequal(size(b),[n,1])
   warning('RSHsizeMatchCoeffMatrix', m);
end
[nz,mz] = size(Z);
if ~isequal(size(b),[nz,1])
   warning(['WrongDeflationMAtrixSize', nz]);
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
Z1(1:n-nw,:)=Z(1:n-nw,:);
M11(1:n-nw,1:n-nw)=M1(1:n-nw,1:n-nw);
M22(1:n-nw,1:n-nw)=M2(1:n-nw,1:n-nw);
b1(1:n-nw,1)=b(1:n-nw,1);
x01(1:n-nw,1)=x0(1:n-nw,1);
A = A1;
Z = Z1;
M1 = M11;
M2 = M22;
b = b1;
x0 = x01;
clear A1 Z1 M11 M22 b1 x01
n = n-nw ;
end


Z  = sparse(Z);
AZ = sparse(A*Z);
E  = Z'*AZ;
EI = inv(E);

%P = sparse(eye(n)-AZ*EI*Z');
[u]=qvec(Z,EI,b);
%u = Z*EI*Z'*b;

if(x_true{1})
    xtrue = A\b;
    normxtrue = norm(xtrue);
end



i = 1;
x = x0;
Mb = M1 \ b;
r = b - A * x;
nr = norm(r)/norm(b)
r = defvec(Z,EI,A,r);
ndr =norm(r)/norm(Mb)
%r = P * r;
y = M1 \ r;
res(i,1) = norm(y)
%p = r;
p = M2 \ y;
nb = norm(b);
nmb = norm(Mb);
tol = tol * nmb

if(x_true{1} )
    xk = x;
    terr(i,1)   = norm(xtrue-xk);
    MAxk = M1 \ (A * xk);
    tresm(i,1) = norm(Mb-MAxk);
    tres(i,1) = norm(b-A*xk);
end



% if(res(i) < tol) 
%    disp(['DICCG only needs one iteration, initial residual is, r_1 = ||P*M^{-1}(b-A*xk)||_2 =' num2str(res(i))])
%    if(x_true{1})
%    disp(['True residual is, r_1 = ||b-A*xk||_2 = ' num2str(tresm(i))])
%    end
% end


while  (i < maxit) && (res(i) > tol)
   
 %   if(Convergence{1}) & (res(i) < Error)
        % If the residual increases, the approximation will be the previous
        % solution
        xacc = x;
  %  end
    
    i = i+1;
    y0 = y;
    r0 = r;
    w = A * p;
    %PAp = P*(A*p);
    PAp = defvec(Z,EI,A,w);
    alpha = (r'*y)/(p'*PAp);
    x = x+alpha*p;
    r = r - alpha * w;
    
    y = M1 \ r;
    
    beta = (r'*y)/(r0'*y0);
    %z = r;
     z = M2 \ y;
    p = z+beta*p;
    res(i,1) = norm(y);
    
    if(x_true{1} )
        [xk]=tdefvec(Z,EI,A,x);
        [Qb] = qvec(Z,EI,b);
        xk = Qb + xk;
        terr(i,1) = norm(xtrue-xk);
        MAxk = M1 \ (A * xk);
        tresm(i,1) = norm(Mb-MAxk);
        tres(i,1) = norm(b-A*xk);
    end

    if(Convergence{1})
         [rc]=tdefvec(Z,EI,A,r);
          rcn = norm(rc);
          E = Error * nmb;
        if ( rcn < E)
            % If the residual increases, the approximation will be the previous
            % solution
            xacc = x;
            % If the residual increases, the approximation will be the previous
            % solution           
            flagr = 1;
            rmin = res(i-1); %
            %         plot(i,res(i),'*')
            %             hold on
            rdf = norm(rmin -res(i));
            srdf = sign(rmin -res(i));
%             figure(1)
%                          plot(i,rdf,'*')
%                           hold on
%             
            %               figure(2)
            %             plot(i,rmin,'*')
            %                          hold on
            
            if (srdf > 0)
                
                rmin = res(i);
                flagr = 1;
                xacc = x;
            else if (srdf < 0 )
                    
                    flagr = 1;
                    %rmin = res(i-1)
                    if (rdf > E)
                        
                        flagr = 0;
                    end
                end
            end
            if flagr == 0
                warning(['Maximum accuracy is : ' num2str(res(i) / nmb)])
                x =  xacc;
                break
            end
            
        end
    end
   
end
if (Iter_m{1} )
disp(['Accuracy is : ' num2str(res(i) / nmb)])
disp(['Number of iterations is: ' num2str(i)])
end
%xk = (u+P'*x);
[xk] = tdefvec(Z,EI,A,x);
[Qb] = qvec(Z,EI,b);
xk = Qb + xk;
tr = norm(b-A*xk)/norm(b);
Mr = M1 \ (b-A*xk);
ptr = norm(Mr)/norm(Mb);


result = xk;
flag = 0;
its = i;
resvec = res;

if (nw > 0)
    for i = 1 : nw
        result(n-nw+i,1) = pw(n-nw+i,1);
    end
end
[n,m] = size(A);


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
if (E_cn{1})
    C_E = cond(E);
    [VE,DE] = eigs(E);
    disp(['The condition number of E is: ' num2str(C_E)])
end
if(PMAmatrix_eigs{1})
Q = Z * EI * Z';
P = sparse(eye(n)-AZ*EI*Z');
[VPMA,DPMA] = eigs(IM*P*A*IM',n);

DPMA = diag(abs(DPMA));
DPMA1 = real(DPMA(mz+1:n));
DPMA1 = abs(DPMA1);
lmax = max(DPMA1);
lmin = min(DPMA1);
CPMA = lmax/lmin;
end

%% Save values for plots
if(E_cn{1})
    resulte.E   =  E;
    resulte.VE  =  VE;
    resulte.DE  =  DE;
end
if(A_cn{1})
    resulte.SZ_A =  SZ_A;
    resulte.C_A  =  C_A;
end

if(x_true{1})
    resulte.xtrue     = xtrue;
    resulte.terr      = terr;
    resulte.tresm     = tresm;
    resulte.tres      = tres;
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
if(PMAmatrix_eigs{1} )
    resulte.VPMA = VPMA;
    resulte.DPMA = (DPMA);
    resulte.CPMA = CPMA;
end
end





function[Qx]=qvec(Z,EI,x)
Qx=Z'*x;
Qx=EI*Qx;
Qx=Z*Qx;
end
function[Px]=defvec(Z,EI,A,x)
[Qx]=qvec(Z,EI,x);
Px=x-A*Qx;
end
function[Px]=tdefvec(Z,EI,A,x)
Ax=A'*x;
[QAx]=qvec(Z,EI,Ax);
Px=x-QAx;
end