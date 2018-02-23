clear all; close all; clc
%% Globals (defined for integration GS)
global Lmat
%% Symbols
syms x t alpha beta nu T w intg;
%% Parameters
nu = sym(0.1);     % 1 / Re
T  = sym(2*pi);  % period for time average
nModes = 2;      % number of modes
eps = 1e-6;      % tolerance for Runge Kutta GS integration
%% Init
v = zeros(nModes,1); v = sym(v);  % temporal coefficients
b = zeros(nModes,1); b = sym(b);  % spatial basis vectors
up = sym(0);                      % velocity field
Q = zeros(nModes,nModes); Q = sym(Q); % matrix needed for eig

%% Compute alpha, beta
alpha = (1/(4*nu))*( sqrt(  2+2*sqrt(1+16*nu^2) ) - 2 );
beta  = (1/(4*nu))*( sqrt( -2+2*sqrt(1+16*nu^2) ) );
%% Define basis functions
% v_1, v_2
v(1) = exp(-alpha*x)*cos(beta*x);
v(2) = exp(-alpha*x)*sin(beta*x);
%% Define temporal coefficients
% b_1, b_2
b(1) = cos(t);
b(2) = sin(t);


%% Construct velocity field
for i = 1:nModes
up = up + b(i)*v(i);
end
%% Compute Q matrix
for i = 1:nModes
for j = 1:nModes
% sum over all modes
for k = 1:nModes
Q(i,j) = Q(i,j) + (1/T)*int(b(i)*b(k),t,0,T)*...
int( v(k)*v(j), x,0,inf );
end
end
end
%% Compute eigenvalues
[eigenVectorArray,eigenValues] = eig(Q);
eigenValues = diag(eigenValues);
[eigenValuesTemp,I] = sort(eval(eigenValues),'descend');
eigenValues = eigenValues(I);
eigenVectorArray = eigenVectorArray(:,I);
%% Assign weights c, u_i = sum_j c_ji v_j
c = eigenVectorArray;
%% Compute POD modes
u = zeros(nModes,1); u = sym(u);
for i = 1:nModes
% sum
for j = 1:nModes
u(i) = u(i) + c(j,i)*v(j);
end
end
% Compute weights and normalized modes
for i = 1:nModes
w = int( u(i)*u(i), x,0,inf );
w = sqrt(w);
u(i) = u(i) / w;
c(:,i) = c(:,i) / w;
end
%% Compute Fourier coefficients
a = zeros(nModes,1); a = sym(a);
for i = 1:nModes
a(i) = int( up*u(i), x,0,inf );
end
%% Compute Galerkin system
ug = sym(0);  % temp. storage
qij = zeros(nModes,nModes); qij = sym(qij);  % convection term
lij = zeros(nModes,nModes); lij = sym(lij);  % dissipation term
for i = 1:nModes
for j = 1:nModes
% Dissipation term: Int dx ( (d_xx u_j) * u_i )
%  Compute 2nd derivative
ug = diff(u(j),'x',2);
%  integrand: (d_xx u_j) * u_i
intg = ug*u(i); % note i-index: only POD modes!
lij(i,j) = int( intg, x,0,inf );
% Convection term: Int dx ( (d_x u_j) . u_i )
ug = diff(u(j),'x',1);
% integrand: (d_x u_j) . u_i
intg = ug*u(i); % note i-index: only POD modes!
qij(i,j) = -int( intg, x,0,inf );
end
end
Lmat = double(simplify(nu*lij+qij));
%% Integrate Galerkin system
a0 = zeros(nModes,1);
for i = 1:nModes
a0(i) = double(subs(a(i),0)); % Get initial condition
end
% Integrate with Runge Kutta
options = odeset('RelTol',eps,'AbsTol',eps*ones(nModes,1));
[tInt,aInt] = ode45(@cde_gs,[0 40],a0,options);
%% Figures
close all
% Plot eigenvalues
figure
plot([1:nModes],eval(eigenValues./sum(eigenValues)),'--o')
grid on
xlabel('i'); ylabel('\lambda_i / \Sigma \lambda_i')
set(gca, 'XTick', [1:nModes]);
title('Eigenvalues')
% Plot Fourier coefficients
figure
p1 = ezplot(a(1),[0 20]);
set(p1,'Color','blue','LineStyle','-')
hold on
p2 = ezplot(a(2),[0 20]);
set(p2,'Color','red','LineStyle','--')
axis([0 20 -1.8 1.8])
title('')
xlabel('t'); ylabel('a_i');
legend('a_1','a_2')
title('Fourier coefficients')
% Plot POD modes
figure
p1 = ezplot(u(1),[0 40]);
set(p1,'Color','blue','LineStyle','-')
hold on
p2 = ezplot(u(2),[0 40]);
set(p2,'Color','red','LineStyle','--')
axis([0 40 -0.9 0.9])
title('')
xlabel('x'); ylabel('u_i');
legend('u_1','u_2')
title('POD modes')
% Plot comparison projected and integrated Fourier coeffic
ients
figure
plot(tInt,aInt,'--x')
hold on
p1 = ezplot(a(1),tInt);
set(p1,'LineStyle','-')
p2 = ezplot(a(2),tInt);
set(p2,'LineStyle','-')
xlabel('t'); ylabel('a_i')
title('Comparison of projected a_i (-) and GS integrated a_i (--x)')




