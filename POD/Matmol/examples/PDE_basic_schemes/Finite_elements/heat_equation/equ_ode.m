%... The MatMol Group (2016)
    function dm = equ_ode(t,m,AA, alpha)

    dm = -alpha*AA*m;