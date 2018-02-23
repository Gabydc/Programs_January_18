%... The MatMol Group (2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implement the boundary conditions in Dirichlet form and add the boundary nodes % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
     function [u z] = Bcintroduct(t,x)
     
     global n ne zL zR
     global v0 eps Deff rhog cpg lamb 
%...
%... separate dependent variables and node positions
%...
     for j=1:ne,
         u(2:n+1,j)=x(j:ne+1:n*(ne+1));
     end
     z=[zL x(ne+1:ne+1:n*(ne+1))' zR];
%...
%... implement the BCs
%...
     u(1,1)=(u(2,1)+v0*(z(2)-z(1))*cin(t)/Deff)/(1+v0*(z(2)-z(1))/Deff);
     u(1,2)=(u(2,2)+eps*v0*rhog*cpg*(z(2)-z(1))*Tin(t)/lamb)/(1+eps*v0*rhog*cpg*(z(2)-z(1))/lamb);
     u(n+2,1)= u(n+1,1)+(z(n+2)-z(n+1))*(u(n+1,1)-u(n,1))/(z(n+1)-z(n)) ;
     u(n+2,2)= u(n+1,2)+(z(n+2)-z(n+1))*(u(n+1,2)-u(n,2))/(z(n+1)-z(n));
