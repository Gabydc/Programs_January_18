     function Jac = jacobian(t,x)
%...
%... Set global variables
     global r0 Clim alpha Cdeath K0 Klim beta H Q delta delta0 eta eta0
%...
%... Transfer dependent variables  
     N = x(1);
     C = x(2);
     E = x(3);
%...
%... Jacobian matrix
%...
%     Nt  = r*N - r0*N^2/K - H;
%     Ct  = Q - delta*C - delta0*E;
%     Et = eta*(C-Clim) - eta0*E;
%...
%...
    if C < Clim 
        r = r0;
        K = K0;
	    Jac = [r-2*r0*N/K    0         0       ;
               0             -delta    -delta0 ;
               0             eta       -eta0  ];
    elseif Clim <= C < Cdeath
	    r = r0-alpha*(C-Clim);
        K = K0-beta*(C-Clim);
        Jac = [r-2*r0*N/K    -alpha*N+beta*r0*N^2/K^2    0       ;
               0             -delta                      -delta0 ;
               0             eta                         -eta0  ];
     else         
        r = 0;
        K = Klim;
	    Jac = [r-2*r0*N/K    0         0       ;
               0             -delta    -delta0 ;
               0             eta       -eta0  ];
     end      
