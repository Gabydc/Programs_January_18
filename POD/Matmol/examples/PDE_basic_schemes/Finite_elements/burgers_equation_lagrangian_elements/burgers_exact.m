%... The MatMol Group (2016)
    function [u] = burgers_exact(z,t)
%... This function computes an exact solution to Burgers' 
%... equation

%... Global variables
    global mu

%... Analytical (exact) solution
    a  = (0.05/mu)*(z-0.5+4.95*t);
    b  = (0.25/mu)*(z-0.5+0.75*t);
    c  = (0.5/mu)*(z-0.375);
    ea = exp(-a);
    eb = exp(-b);
    ec = exp(-c);
    u  = (0.1*ea+0.5*eb+ec)/(ea+eb+ec);
