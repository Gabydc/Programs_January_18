%... The MatMol Group (2016)
    function xt = bacteria_odes(t,x)
%...
%... Set global variables
     global nu mumax Km Ki Kc kinetics
%...
%... Transfer dependent variables  
     X = x(1);
     S = x(2);
%...      
%... Temporal derivatives
%...
     switch kinetics    
%...
%... constant
       case('constant')
           mu = mumax;
%... monod
       case('monod')
           mu = mumax*S/(Km+S);
%... haldane
       case('haldane')
           mu = mumax*S/(Km+S+(S^2)/Ki);
%... contois
       case('contois')
           mu = mumax*S/(Kc*X+S);
%...
     end
%...
     phi = mu*X;
%...
     Xt  = phi;
     St  = -nu*phi;
%...
%... Transfer temporal derivatives
     xt = [Xt St]';
