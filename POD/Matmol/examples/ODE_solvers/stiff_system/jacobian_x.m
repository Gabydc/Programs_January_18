%... The MatMol Group (2016)
     function Jac = jacobian_x(t,x)
%...
%... Set global variables
     global a b
%...
%... Jacobian matrix
	  Jac = [-a   b;
              b  -a];
