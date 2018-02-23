%... The MatMol Group (2016)
     function [x] = burgers_exact(z,t)

%...
%... this function computes an exact solution to Burgers' equation
%...
     global visc
     a(1)=-(0.05/visc)*(z-0.5+4.95*t);
     a(2)=-(0.25/visc)*(z-0.5+0.75*t);
     a(3)=-(0.5/visc)*(z-0.375);
%...
     ea=0;
     eb=0;
     ec=0;
     temp=max(a);
     if a(1)-temp >= -35
         ea=exp(a(1)-temp);
     end
     if a(2)-temp >= -35
         eb=exp(a(2)-temp);
     end
     if a(3)-temp >= -35
         ec=exp(a(3)-temp);
     end
%...     
     x=(0.1*ea+0.5*eb+ec)/(ea+eb+ec);
