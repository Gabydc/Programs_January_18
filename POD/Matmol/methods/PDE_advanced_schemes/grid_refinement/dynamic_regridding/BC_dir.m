%... The Matmol group (2016)
    function ubc = BC_dir(bc_form, u, uin, D1, pos, v, Dif)
% This function transforms Dirichlet, Neumann and Robin boundary condition
% types into a Dirichlet one.
%
% Input parameters:
%   bc_form : Type of the boundary condition 'dir' for Dirichlet boundary
%             type, 'neu' for Neumann type and 'rob' for Robin type
%   u       : Value of the field in the spatial domain
%   uin     : Value of the non homogeneous part of the boundary. For
%             homogeneous boundary conditions, set this value to 0.
%   D1      : Finite differences matrix for the first derivative
%   pos     : Position of the boundary condition: 'begin' if the boundary
%             condition refers to the first point of the spatial domain,
%             'end' if it refers to the last point.
%   v       : Fluid velocity. This parameter is only necessary in the case
%             of Robin boundary conditions, otherwise it can be set to 0
%   Dif     : Difussivity. This parameter is only necessary in the case of
%             Robin boundary conditions, otherwise it can be set to 0 
%
% Output parameters:
%   ubc     : Value of the field in the boundary (Dirichlet form)
%
% Example: Tranform the following Dankwerts boundary conditions into the
%          Dirichlet form:
%          D*uz = v(u-uin)   in z=0
%            uz = 0          in z=L=0.1
%          with uin = 5, v=3, D=0.001
%
%          global n
%          n  = 19;  uin = 5;      z = linspace(0,0.1,n+2);
%          v  = 3;   D   = 0.001;  u = sin(z)';  
%          D1 = matfd (z,'non_uni','1st',5,'centered');
%          u0 = BC_dir('rob', u, uin, D1, 'begin', v, D);
%          uL = BC_dir('neu', u, 0, D1, 'end', 0, 0);

% Global variables
    global n

% Number of discretization points
    ndisc = n + 2;


    switch bc_form
        case {'dir'}
            ubc = uin;
        case {'neu'}
            switch pos
                case {'begin'}
                    ubc = ( uin - D1(1,2:ndisc-1)*u(2:ndisc-1,1) )/D1(1,1);
                case {'end'}
                    ubc = ( uin - D1(ndisc,2:ndisc-1)*u(2:ndisc-1,1) )/D1(ndisc,ndisc);
                otherwise
                    'Error: Incorrect boundary point'
            end
        case {'rob'}
            switch pos
                case {'begin'}
                    ubc = ( v/Dif*uin + D1(1,2:ndisc-1)*u(2:ndisc-1,1) )/(v/Dif-D1(1,1));
                case {'end'}
                    ubc = ( v/Dif*uin + D1(ndisc,2:ndisc-1)*u(2:ndisc-1,1) )/(v/Dif-D1(ndisc,ndisc));
                otherwise
                    'Error: Incorrect boundary point'
            end
    end