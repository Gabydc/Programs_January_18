%... The MatMol Group (2016)
    function Tt = Graetz_pde(t,T)
%...
%... set global variables
    global Pe v r_inv nv nz nr dr D1z D1r D2z D2r
%...
%... boundary conditions
    T(1:nz,1) = T(nz+1:2*nz,1);        % BC along r = 0
    T(nv-nz+1:nv,1) = 1;               % BC along r = R
    T(1:nz:nv-nz+1,1) = 0;             % BC along z = 0
    T(nz:nz:nv,1) = T(nz-1:nz:nv-1,1); % BC along z = zL
%...
%... spatial derivatives
    [Tz Tr] = first_order_derivatives_2D(T,nz,nr,D1z,D1r);
    [Tzz Trr] = second_order_derivatives_2D(T,nz,nr,D2z,D2r);
%...
%... temporal derivatives
    Tt = -Pe*v.*Tz + Tzz + Trr;
    Tt(nz+1:nz*nr,1) = Tt(nz+1:nz*nr,1) + r_inv(nz+1:nz*nr,1).*Tr(nz+1:nz*nr,1);
    Tt(1:nz,1) = Tt(1:nz,1) + Trr(1:nz,1);