%...  The MatMol Group (2016)
%...  Example GNI problem file for solving the 2D Lotka-Volterra writting in
%...  Hamiltonian formulation
% 
%...  Use it as follows:
%
%...      gni_fgauss1('F_lotka');
%
    function [out,out2,out3]=F_lotka(~,y,flag,varargin) 
    if (nargin < 3) || isempty(flag)
        out(1,:) = exp(y(2,:)) - 2;
        out(2,:) = -exp(y(1,:)) + 1;
    else
        switch flag
            case 'init',
                out = [0 50];
                out2 = [log(2) log(2)];
                out3 = gni_set('Precision',0.02,'StepSize',0.02,'Vectorized','on');
        end
    end