%...  The MatMol Group (2016)
%   Example GNI problem file. Two bodies problem with Lennard-Jones potential
%   Use it as follows:
%
%   gni_method('G_lennard');    
% 
%   With gni_method replaced by either gni_fstrvl2, gni_vstrvl2 or gni_flmm2
    function [out,out2,out3] = G_lennard(~,q,flag,varargin)   
    if (nargin < 3) || isempty(flag)
        out = 12*q^-13 - 6*q^-7;
    else
        switch flag
            case 'init',
                out = [0 100];
                out2 = [2 0];
                out3 = gni_set('Precision',0.01,'StepSize',0.01);
        end
    end