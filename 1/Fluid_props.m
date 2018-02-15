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
    nf = nf + 1;
    hcp=figure(nf);
    plot(xDummy.s, pc);
    xlabel('s_w'); ylabel('pc [bar]');
    title('Capillary pressure curve')
    file{nf} = ['Capillary pressure'];
end

% Plot relative permeability
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);