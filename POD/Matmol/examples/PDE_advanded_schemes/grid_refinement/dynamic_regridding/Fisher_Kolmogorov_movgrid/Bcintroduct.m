%... The MatMol Group (2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implement the boundary conditions in Dirichlet form and add the boundary nodes % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
     function [u z] = Bcintroduct(t,x)
     
     global n ne zL zR
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
     u(1,1)= 1 ;
     u(n+2,1)= -1 ;
