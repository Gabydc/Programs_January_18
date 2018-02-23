%...  The MatMol Group (2016)
%...   Example GNI problem file for solving the 2D Perturbed Kepler problem.
%...   Use it as follows:
%
%...       gni_method('G_kepler_perturb',[],[],[],0.6);
%... 
%...   With gni_method replaced by either gni_fstrvl2, gni_vstrvl2 or
%...   gni_flmm2 to draw the solution with eccentricity 0.6.
    function [out,out2,out3]=G_kepler_perturb(~,q,flag,varargin)
    mu=0.015;
    
    if (nargin < 3) || isempty(flag)
        rad = q(1,:).*q(1,:) + q(2,:).*q(2,:);
        rad = sqrt(rad);
        rad1 = rad.^3;
        rad2 = rad.^5;
        out(1,:) = -q(1,:)./rad1 - mu*q(1,:)./rad2;
        out(2,:) = -q(2,:)./rad1 - mu*q(2,:)./rad2;
    else
        switch flag
            case 'init',
                if (nargin < 4)
                    ecc = 0.6;
                else
                    ecc = varargin{1};
                end
                if (ecc < 0) || (ecc >= 1)
                    error('The eccentricity must lie between 0 and 1');
                end
                out = [0 120];
                out2 = [1-ecc,0,0,sqrt((1+ecc)/(1-ecc))];
                out3 = gni_set('Precision',0.1,'StepSize',0.01,'Vectorized','on');
        end
    end