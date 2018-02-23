%...  The MatMol Group (2016)
%... 
%... Example GNI problem file for solving the 2D Kepler problem. To compute
%... and show the solution with eccentricity 0.5. 

%... Solve non-stiff differential equations
     gni_fgauss1('F_kepler',[],[],[],0.5);
     
%... Solve non-stiff autonomous differential equations
     gni_vgauss1('F_kepler',[],[],[],0.5);

   