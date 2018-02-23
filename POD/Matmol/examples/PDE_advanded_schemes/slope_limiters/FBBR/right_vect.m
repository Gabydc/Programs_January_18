%...  The MatMol Group (2016)
    function xt = right_vect(t,x)

%... Compute the right vector of the global ODEs system
%...  The following code implements the right vector b(x,t) of the system
%
%...            A(x,t)xdot = b(x,t)
%
%...  with
%...          1   2         ne      1   2         ne            1         ne
%...    x = [u   u   ...   u   z   u   u   ...   u   z   ...   u   ...   u   z ]'
%...          1   1         1   1   2   2         2   2         n         n   n
%
%... and where the meshpoints 1, 2, ... n are the interiors nodes. Nevertheless, b(x,t) uses also
%... the boundary points and values z , z   , u  , u   .
%...                                 0   n+1   0    n+1
%
%... The code first implements the boundary conditions and computes the monitor function.
%
%... Then, following elements of b(x,t) are then implemented :
%
%...    a) elements from the moving grid equation :
%
%...            1         mu        1+2mu       mu          1         mu        1+2mu       mu
%...    g(i) = --- [ - --------- + ------- - -------- ] - ----- [ - ------- + --------- - -------- ]
%...           M(i)    dltz(i+1)   dltz(i)   dltz(i-1)    M(i-1)    dltz(i)   dltz(i-1)   dltz(i-2)
%
%...    b) elements from the PDEs :
%
%...    udot(i,j)   i = 1, ..., n       j = 1, ..., ne.
%
%... Structure of b(x,t) : b(x,t) is (ne+1)*n vector :
%
%...          | udot(1,1)  |
%...          | udot(1,2)  |
%...          |   . . .    |
%...          | udot(1,ne) |
%...          |   g(1)     |
%...          | udot(2,1)  |
%...          | udot(2,2)  |
%...          |  . . .     |
%... b(x,t) = | udot(2,ne) |
%...          |   g(2)     |
%...          |  . . .     |
%...          | udot(n,1)  |
%...          | udot(n,2)  |
%...          |   . . .    |
%...          | udot(n,ne) |
%...          |   g(n)     |

%... Global variables
    global n ne
    global mu
    global u z dltz mon
    global choice
    global D1 mu0 k0 kd KS KI

%... Separate dependent variables and node positions and implement the BCs
    [u z] = Bcintroduct(t,x);

%... Compute the monitoring function mon(i)
    mon = monitor(u,z,choice,t);

%... Compute the right vector g(n) of the moving grid equation
    g = rightmon(mon,dltz,mu);

%... Spatial derivatives
    Dzz = matfd (z, 'non_uni', '2nd', 3, 'centered');
    fz  = kurg_centred_slope_limiter_fz(ne,n+2,z',t,u,'flux',...
                                        'dflux_dx');

%... Independent variables
    S = u(:,1);
    X = u(:,2);

%... Reaction rates
    muA = mu0*S./(KS*X + S + 1/KI*S.^2);

%... Temporal derivatives
    St = D1*Dzz*S - fz(:,1) - k0*muA.*X;
    Xt =                    - kd*X + muA.*X;

    udot(:,1) = St;
    udot(:,2) = Xt;     

%... Reorder the global right vector
    for jj = 1:ne
        xt(jj:ne+1:n*(ne+1),1)=udot(2:n+1,jj);
    end
    xt(ne+1:ne+1:n*(ne+1),1)=g;

    end 