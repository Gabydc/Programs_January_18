%... The MatMol Group (2016)
    function ut = kdv3_pde(t,u)

%... Set global variables
     global a b c
     global z0 zL z n D1 D3
%...
%... Spatial derivatives and boundary conditions
%...
%... stagewise differentiation
     uz = D1*u;
     uzz = D1*uz;
     uzzz = D1*uzz;
%... direct differentiation
%    uz = D1*u;
% 	 uzzz = D3*u;
%...
%... Partial differential equations
%...
     ut = -a*u.*uz-b*uzzz;
