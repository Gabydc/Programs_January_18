clear all
close all
clc
x = 3;
y = 3;
l = 3;
s0 = 1;
s  = 1;
n  = x*y*l-2*x;
[a,b,z] = matrixf(x,y,l,s0,s);
xi(1:size(b,1)) = rand;
maxIter = 500;
tol = 10^-8;





