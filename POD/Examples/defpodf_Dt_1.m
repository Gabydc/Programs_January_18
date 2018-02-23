%This program gives the matrix for the solution of the Laplace problem, the
%linear system is solved with the CGF conjugate gradient function

function [U,S]=defpodf_Dt(x,dir,dvec,step,dTplot)
%x=zupd;
for i=1:size(x,2)
    x(:,i)=x(:,i)/norm(x(:,i));
end
lx=size(x,1);
ly=size(x,2);
D=x'*x;
[V,L] = eigs(D,ly);
%size(D)
%size(V)
%size(L)
S=zeros(lx,ly);
g=diag(L);
g1=sqrt(abs(g));
%g1=1./g;
S(1:ly,1:ly)=diag(g1);
mg=max(g1);
S=sparse(S);
%size(S)
V=sparse(V);
%size(V)
x=sparse(x);
%size(x)
U=x*V*S';
% size(U)
 xax=1:ly;

%[V,D] = eigs(X,n);
        if mod(step,dTplot)==0
            figure(3000+step)
%            hd=plot(xax(1:step),log(g1(1:step)),'ob');
            hd=plot(log(g1),'ob');
            axis('tight')
            title(['Eigenvalues R=X*X^T, t = ' num2str(step) ' days'],'FontSize',16);
            ylabel('log(Value) ','FontSize',16)
            xlabel('Eigenvalue','FontSize',16)
            
            % figure(4000)
            % hd1=plot(xax,log(g1/mg),'og');
            % axis('tight')
            % title('Eigenvalues R=Z*Z^T','FontSize',16);
            % ylabel('log(Value) ','FontSize',16)
            % xlabel('Eigenvalue','FontSize',16)
            % hold on
            
            
            file='eig_pod';
            B=[dir file num2str(step) '.fig'];
            saveas(hd,B)
            B=[dir  file num2str(step) '.jpg'];
            saveas(hd,B)
        end

