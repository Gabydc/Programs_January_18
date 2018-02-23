%...  The MatMol Group (2016)
%...  Example GNI problem file for solving the 2D Kepler problem.
%...  Use it as follows:
%
%...      gni_fgauss1('F_kepler',[],[],[],0.5); or
%...      gni_vgauss1('F_kepler',[],[],[],0.5);
%
%...  to compute and show the solution with eccentricity 0.5.
    function [out,out2,out3]=F_kepler(~,y,flag,varargin) 
    
    if (nargin < 3) || isempty(flag)
        rad = y(1,:).*y(1,:)+y(2,:).*y(2,:);
        rad = rad .* sqrt(rad);
        out(1,:) = y(3,:);
        out(2,:) = y(4,:);
        out(3,:) = -y(1,:)./rad;
        out(4,:) = -y(2,:)./rad;
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
            out = [0 100];
            out2 = [1-ecc 0 0 sqrt((1+ecc)/(1-ecc))];
            out3 = gni_set('Precision',0.1,'StepSize',0.05,'Vectorized','on');
        end
    end