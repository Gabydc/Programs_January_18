%... The MatMol Group (2016)
    function Tt = heat_equation_pde(t,T)
%...
%... set global variables
    global k TM 
    global nx ny D2x D2y
%...
%... boundary conditions
    T(1:nx,1) = 0;                    % BC along y = 0
    T(1:nx:nx*(ny-1)+1,1) = 0;        % BC along x = 0
    T(nx:nx:nx*ny,1) = TM;            % BC along x = xL
    T(nx*(ny-1)+1:nx*ny,1) = TM;      % BC along y = yL
%...
%... spatial derivatives
    [Txx Tyy] = second_order_derivatives_2D(T,nx,ny,D2x,D2y);
%...
%... temporal derivatives
    Tt = k*(Txx+Tyy);
