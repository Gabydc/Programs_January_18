function [i,j]=Ind2sub(I,Nx)
i = mod(I-1,Nx)+1;
j = (I - i ) / Nx +1;
end
