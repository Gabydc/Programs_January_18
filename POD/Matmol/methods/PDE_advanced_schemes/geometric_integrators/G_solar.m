%...  The MatMol Group (2016)
    function [out,out2,out3] = G_solar(~,q,flag,varargin)

%   Outer Solar System

%   REFERENCE:
%   E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration. 
%   Structure-Preserving Algorithms for Ordinary Diferential Equations,
%   springer-verlag, berlin, second edition Edition, Vol. 31, Springer
%   Series in Computational Mathematics 31, 2006.
%
    if (nargin < 3) || isempty(flag)
        G  = 2.95912208286e-4;
        m0 = 1.00000597682;
        m1 = 9.54786104043e-4;
        m2 = 2.85583733151e-4;
        m3 = 4.37273164546e-5;
        m4 = 5.17759138449e-5;
        m5 = 1/(1.3e8);
        
        q0 = q(1:3);
        q1 = q(4:6);
        q2 = q(7:9);
        q3 = q(10:12);
        q4 = q(13:15);
        q5 = q(16:18);
        
        Q01 = norm(q0-q1)^3;
        Q02 = norm(q0-q2)^3;
        Q03 = norm(q0-q3)^3;
        Q04 = norm(q0-q4)^3;
        Q05 = norm(q0-q5)^3;
        Q12 = norm(q1-q2)^3;
        Q13 = norm(q1-q3)^3;
        Q14 = norm(q1-q4)^3;
        Q15 = norm(q1-q5)^3;
        Q23 = norm(q2-q3)^3;
        Q24 = norm(q2-q4)^3;
        Q25 = norm(q2-q5)^3;
        Q34 = norm(q3-q4)^3;
        Q35 = norm(q3-q5)^3;
        Q45 = norm(q4-q5)^3;
        
        out = zeros(18,1);
        out(1:3)   = - G*(m1/Q01*(q0-q1)+m2/Q02*(q0-q2)+m3/Q03*(q0-q3)+m4/Q04*(q0-q4)+m5/Q05*(q0-q5));
        out(4:6)   = - G*(m0/Q01*(q1-q0)+m2/Q12*(q1-q2)+m3/Q13*(q1-q3)+m4/Q14*(q1-q4)+m5/Q15*(q1-q5));
        out(7:9)   = - G*(m0/Q02*(q2-q0)+m1/Q12*(q2-q1)+m3/Q23*(q2-q3)+m4/Q24*(q2-q4)+m5/Q25*(q2-q5));
        out(10:12) = - G*(m0/Q03*(q3-q0)+m1/Q13*(q3-q1)+m2/Q23*(q3-q2)+m4/Q34*(q3-q4)+m5/Q35*(q3-q5));
        out(13:15) = - G*(m0/Q04*(q4-q0)+m1/Q14*(q4-q1)+m2/Q24*(q4-q2)+m3/Q34*(q4-q3)+m5/Q45*(q4-q5));
        out(16:18) = - G*(m0/Q05*(q5-q0)+m1/Q15*(q5-q1)+m2/Q25*(q5-q2)+m3/Q35*(q5-q3)+m4/Q45*(q5-q4));
    else
        switch flag
            case 'init',
                out = [0 150000];
                
                q0 = [0 0 0];
                q1 = [-3.5023653 -3.8169847 -1.5507963];
                q2 = [9.0755314 -3.0458353 -1.6483708];
                q3 = [8.3101420 -16.2901086 -7.2521278];
                q4 = [11.4707666 -25.7294829 -10.8169456];
                q5 = [-15.5387357 -25.2225594 -3.1902382];
                p0 = [0 0 0];
                p1 = [0.00565429 -0.00412490 -0.00190589];
                p2 = [0.00168318 0.00483525 0.00192462];
                p3 = [0.00354178 0.00137102 0.00055029];
                p4 = [0.00288930 0.00114527 0.00039677];
                p5 = [0.00276725 -0.00170702 -0.00136504];
                out2 = [q0 q1 q2 q3 q4 q5 p0 p1 p2 p3 p4 p5]';
                
                out3 = gni_set('Precision',1,'StepSize',50);
        end
    end