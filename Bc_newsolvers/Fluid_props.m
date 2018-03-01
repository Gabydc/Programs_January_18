vx = linspace(0, 1, 11) .';
vy = linspace(1, 0, 11) .';

pc_form = 'nonwetting';
[muw, muo] = deal(1,10);
[rhow, rhoo] = deal(1000,700);
[krw, kro] = deal(2,2);
cap_scale = 10;
props = constantProperties([   muw,  muo] .* centi*poise, ...
                           [rhow, rhoo] .* kilogram/meter^3);  
[kr, pc]  = tabulatedSatFunc([vx, vx.^krw, vy.^kro, vy.*cap_scale*barsa]);

if(~use_cp)
    fluid = initSimpleFluid('mu' , [   muw,    muo] .* centi*poise     , ...
        'rho', [rhow, rhoo] .* kilogram/meter^3, ...
        'n'  , [   krw,    kro]);
else
    
    fluid = struct('properties', props                  , ...
        'saturation', @(x, varargin)    x.s  , ...
        'relperm'   , kr                     , ...
        'pc'        , @(x, varargin) pc(x.s));
    xDummy   = initState(G, [], [0, 1]);
    xDummy.s = linspace(0, 1, numel(xDummy.s))'; ...
        pc = convertTo(fluid.pc(xDummy), barsa);

% title = 'Capillary pressure';
% xlab = 'S_w'; 
% ylab = 'Pc [bar]';
% nf = 2;
%     Plot_2(pc,'x',xDummy.s,'dnames',{{title}},'linestyles','-','o_ax',...
%     2, 'x_lab',xlab,'y_lab', ylab,'figure',nf-1)
end

% Plot relative permeability
x_kr = linspace(0, 1, 1001).'; kr = fluid.relperm(x_kr);
