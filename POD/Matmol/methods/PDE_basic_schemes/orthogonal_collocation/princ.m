%...  The Matmol group (2016)
%...
tic
close all
clear all
global mu
global n dz z0 zL D2 alp0 alp1

% Spatial grid
z0  = 0;
zL  = 1;
n   = 101;
nel = n-1;
dz  = (zL - z0)/(n-1);
z   = (z0:dz:zL)';

% Differentiation matrices
D2 = colorthD2(n,dz);

% Problem constants
mu = 0.001;
ne = 1;

% Initial conditions
u  = zeros(1,n);
uz = zeros(1,n);
for i = 1:n
    u(i)  = burgers_exact(z(i),0);
    uz(i) = derburgers_exact(z(i),0);
end
x(1:2:2*n-1) = u;
x(2:2:2*n)   = dz*uz/2;

% Visualization instants
dt   = 0.1;
time = (0:dt:1);
nt   = length(time);

% Time integration
options        = odeset('Mass',mcolorth(n,ne));
[timeout,yout] = ode15s(@burgerspdes,time,x,options);
err = 0;
figure
hold
ztheor = (0:.001:1);
ntheor = length(ztheor);
for k = 1:length(timeout)
    for i = 1:ntheor
        yexact(i) = burgers_exact(ztheor(i),timeout(k));
    end
    plot(ztheor,yexact,'r')
end

% Plot the resutls
for i=1:n
    yexactn(i)= burgers_exact(z(i),timeout(nt));
end
err = norm(yout(nt,1:2:2*n-1)-yexactn);


axis([0 1 0 1.4])
plot(z,yout(:,1:2:2*n-1),'g')
figure
hold
plot(z,yout(:,1:2:2*n-1),'.k')
dzvis = 0.4;
zvis  = (-1:dzvis:1);
nvis  = length(zvis);
for i = 1:length(timeout)
    for k = 1:nel
        for j = 1:nvis
            absc(j)       = (k-1)*dz+(j-1)*dzvis*dz/2;
            [N1 N2 N3 N4] = trialfunctionsher_1(zvis(j));
            ordo(j)       = yout(i,2*k-1)*N1 +...
                            yout(i,2*k)*N2 + ...
                            yout(i,2*k+1)*N3 + ...
                            yout(i,2*k+2)*N4;
        end
        plot(absc,ordo,'-k');
    end
end
ztheor = (0:.001:1);
ntheor = length(ztheor);
for k = 1:length(timeout)
    for i = 1:ntheor
        yexact(i) = burgers_exact(ztheor(i),timeout(k));
    end
    plot(ztheor,yexact,'r')
end
axis([0 1 0 1.4])
toc