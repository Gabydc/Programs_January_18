%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function

function [U,S]=PODbasis_norm_s(X)
%  X are the snapshots
[lx,ly] = size(X);
mx = max(X');
nmx =norm(mx);
% Normalize the snapshots
for i=1:ly
    X(:,i)=X(:,i)/nmx; 
end
%Compute mean of the snapshots
xm = mean(X,2);
% if ly == 2
% else
% for i=1:ly
%     x(:,i)=x(:,i)-xm;
% end
% end
D=(1/ly)*(X'*X);
D = sparse(D);
[V,L] = eigs(D,ly);
S=zeros(lx,ly);
g=diag(L);
g=sqrt(abs(g));
S(1:ly,1:ly)=diag(g);
S=sparse(S);
%size(S)
V=sparse(V);
%size(V)
X=sparse(X);
%size(x)
U=X*V*S';
U = U(:,1:ly);
S = diag(S);
for i=1:ly
    U(:,i)=U(:,i)/norm(U(:,i)); 
end
ru = rank(full(U))
U1 = U(:,8:10);
rank(full(U1))
% figure
% plot(log((diag(S))),'*r');
% nf =nf+1;
% f(nf) = figure(nf);
% plot((U(:,91:100)));
% figure; 
% for i = 91 : 100
%     plot(U(:,i)); pause;
%     hold on
% end
% norm(U(:,100))
% size(U)

