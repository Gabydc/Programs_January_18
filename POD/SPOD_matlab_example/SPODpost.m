function mode = SPODpost(a)
% mode = SPODpost(a)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  mode = SPODpost(a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds mode pairs from a Spectral Proper Orthogonal 
% Decomposition.
% 
% Function input:
% a:    SPOD mode coefficients -> output of SPOD()
%
% Function output:
% mode: Structure where single element (for one mode pair k) consists of:
%       mode(k).ind: Indices [i,j] of the coressponding SPOD modes
%       mode(k).c: Magnitude of the harmonic corellation
%       mode(k).a: Analytic mode coefficient a(:,i) + 1i*a(:,i)
%       mode(k).at: Derivative of analytic mode coefficient
%       mode(k).K: Relative energy content of a single mode pair
%       mode(k).f: Estimated frequency of the mode pair [1/sample]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Moritz Sieber
% Version: 03/2017
% E-Mail: moritz.sieber@fd.tu-berlin.de
%
% Please cite our article if you use this method in your own work:
%
% Sieber, M., Paschereit, C.O. and Oberleithner, K. (2016) ‘Spectral proper 
% orthogonal decomposition’, Journal of Fluid Mechanics, 792, pp. 798–828. 
% doi: 10.1017/jfm.2016.103.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nsnap,Npod] = size(a);
lambda = diag(a'*a);
lambda = lambda/sum(lambda);

[~,ag] = gradient(a);

% find combined modes from DMD 
A = a(2:Nsnap,:)'/a(1:Nsnap-1,:)';
[V,mu] = eig(A);
mu = diag(mu);

% Compute frequency
omega = imag(log(mu));

% linked modes are only phase shifted -> correlate real and imaginary part
% of DMD eigenvectors
C = imag(V*diag(sign(omega))*V')/2;

% pick out single mode pairs from C
Nmode = ceil(Npod/2);
mode = struct([]);
for i=1:Nmode
    % pick maximum
    [tmp,ind] = max(C(:));
    mode(i).c = tmp;
    [indi,indj] = ind2sub([Npod,Npod],ind);
    % delete row and column
    C([indi indj],:)=0;
    C(:,[indi indj])=0;
    % save indices to output
    mode(i).ind = [indi indj];
    
    % compose analytical mode coefficient
    mode_sign = sign(sum(a(:,indi).*ag(:,indj)));
    mode(i).a = a(:,indi) + mode_sign*1i*a(:,indj);
    mode(i).at = ag(:,indi) + mode_sign*1i*ag(:,indj);
    
    % calculate common energy content
    mode(i).K = lambda(indi)+lambda(indj);

    % estimade mode frequency
    V_comb = abs(V(indi,:)).^2+abs(V(indj,:)).^2;
    [tmp,ind] = sort(V_comb,'descend');
    sub = 1:6; %sum first 3 frequencies (each appears twice)
    mode(i).f = sum(abs(omega(ind(sub))).*tmp(sub)')...
        /sum(tmp(sub))/(2*pi);
end


