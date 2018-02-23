%... The MatMol Group (2016)
     function xt = fish_odes(t,x)
%...
%... Set global variables
     global r0 Clim alpha Cdeath K0 Klim beta H Q delta delta0 eta eta0
%...
%... x has three columns corresponding to three different column solution
%... vectors (which can be used for numerical evalution of the Jacobian):
     ncols = size(x,2);
%...
     for j=1:ncols   
%...
%... Transfer dependent variables  
        N = x(1,j);
        C = x(2,j);
        E = x(3,j);
%...      
%... Temporal derivatives
%...
        if C < Clim
            r = r0;
        elseif Clim <= C < Cdeath
            r = r0-alpha*(C-Clim);
        else
            r = 0;
        end   
%...
        if C < Clim
            K = K0;
        elseif Clim <= C < Cdeath
            K = K0-beta*(C-Clim);
        else
            K = Klim;
        end   
%...
        Nt  = r*N - r0*N^2/K - H;
        Ct  = Q - delta*C - delta0*E;
        Et = eta*(C-Clim) - eta0*E;
%...
%... Transfer temporal derivatives
%... (One column for each column of x)
        xt(:,j) = [Nt Ct Et]';
%...           
     end;