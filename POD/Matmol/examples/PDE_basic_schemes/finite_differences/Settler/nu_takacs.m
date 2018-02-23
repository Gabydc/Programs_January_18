%... The MatMol Group (2016)
    function nu = nu_takacs(C)

%... Computation of the settling velocity following Takacs' law
    global nu0 rp rh nuprim0 Cmin

    for i=1:length(C)
        if C(i) <= Cmin
            nu(i) = 0;
        else
            nu(i) = nu0*( exp(-rh*(C(i)-Cmin)) - exp(-rp*(C(i)-Cmin)) );
        end
    end

    nu(nu>nuprim0) = nuprim0;
    nu(nu<0)       = 0;
