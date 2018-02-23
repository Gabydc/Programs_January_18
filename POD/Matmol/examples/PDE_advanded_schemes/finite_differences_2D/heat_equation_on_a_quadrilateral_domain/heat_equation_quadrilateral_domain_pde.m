%... The MatMol Group (2016)
    function Tt = heat_equation_quadrilateral_domain_pde(t,T)
%...
%... set global variables
    global k TM 
    global nksi neta D1_ksi D1_eta D2_ksi D2_eta 
    global ksi_x ksi_y eta_x eta_y ksi_xx ksi_yy eta_xx eta_yy
    t
%...
%... boundary conditions
    T(1:nksi,1) = 0;
    T(1:nksi:nksi*(neta-1)+1,1) = 0;
    T(nksi:nksi:nksi*neta,1) = TM;
    T(nksi*(neta-1)+1:nksi*neta,1) = TM;
%...
%... spatial derivatives
    [Tksi Teta] = first_order_derivatives_2D(T,nksi,neta,D1_ksi,D1_eta);
    [Tksiksi Tetaeta] = second_order_derivatives_2D(T,nksi,neta,D2_ksi,D2_eta);
    Tksieta = mixed_second_order_derivatives_2D(T,nksi,neta,D1_ksi,D1_eta);
%...
    Txx = Tksiksi.*(ksi_x.^2) + 2*Tksieta.*ksi_x.*eta_x + Tetaeta.*(eta_x.^2) + Tksi.*ksi_xx + Teta.*eta_xx;
    Tyy = Tksiksi.*(ksi_y.^2) + 2*Tksieta.*ksi_y.*eta_y + Tetaeta.*(eta_y.^2) + Tksi.*ksi_yy + Teta.*eta_yy;
%...
%... temporal derivatives
    Tt = k*(Txx+Tyy);
