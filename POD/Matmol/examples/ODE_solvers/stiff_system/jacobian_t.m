%... The MatMol Group (2016)
     function Ft = jacobian_t(t,x)
%...
%... Set global variables
     global a b
%...
%... Stiff system is autonomous
	  Ft = [0;0];
