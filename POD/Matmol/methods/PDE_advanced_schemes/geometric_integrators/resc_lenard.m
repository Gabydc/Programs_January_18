%...  The MatMol Group (2016)
    function out = resc_lenard(~,q,~,varargin)
%...  Example for changing the rescaling function in gni_fstrvl2 or
%...  gni_vstrvl2. The function by default is out=norm(f(y)).
%...  
%...  Use it as follow:
%...      options = gni_set('RescalFun','resc_lenard');
%...      gni_method('F',tspan,y0,options);
%...      with gni_method remplaced by gni_fstrvl2 or gni_vstrvl2.

    out = 1 + norm(q)^-2;