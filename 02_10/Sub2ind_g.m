 function [I]=Sub2ind(i,j,z,Nx,Ny,Nz)
ncell = 0;
ni = numel(i);
nj = numel(j);
nz = numel(z);


 for zz = 1: nz
for jj = 1 : nj
for ii = 1 : ni
ncell = ncell + 1;
I(ncell) = i(ii)+(j(jj)-1)*Nx+(z(zz)-1)*(Nx*Ny);
end
end
end
end
