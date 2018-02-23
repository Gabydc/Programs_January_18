%... The MatMol Group (2016)
    function xt = bstr_odes(t,x)
%...
%... Set global variables
     global rhocp Cd Cr Ed Er DH R
%...
%... Transfer dependent variables  
     A = x(1);
     B = x(2);
     T = x(3);
%...      
%... Temporal derivatives
%...
     kd = Cd*exp(-Ed/(R*T));
     kr = Cr*exp(-Er/(R*T));
%...
     At  = - kd*A + kr*B;
     Bt  = + kd*A - kr*B;
     Tt  = + (kd*A - kr*B)*DH/(rhocp);
%...
%... Transfer temporal derivatives
     xt = [At Bt Tt]';
