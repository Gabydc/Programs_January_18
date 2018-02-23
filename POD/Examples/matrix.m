%%  Matrix_POD
%Create matrix to solve with POD


x = linspace(0,1,25);
t = linspace(0,2,50);
[X, T] = meshgrid(x,t);
Z = exp(-abs(X - 0.5) .* (T-1)) + sin(X.* T);
