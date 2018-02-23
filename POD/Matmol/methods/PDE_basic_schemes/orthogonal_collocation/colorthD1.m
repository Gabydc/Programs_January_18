%...  The Matmol group (2016)
     function D1 = colorthD1(n,h)
%...
%... Computation of the finite element matrix D1 with the
%... orthogonal collocation technique
    a = sqrt(3);
    b = [-1 1/a 1 -1/a ; -1 -1/a 1 1/a]/2; 
    base(1:2,:) = horzcat(b, zeros(2,2*n-4));
    for i = 1:n-2
        base(2*i+1:2*i+2,:) = circshift(base(2*i-1:2*i,:),[0 2]);
    end

    D1 = [zeros(1,2*n) ; base ; zeros(1,2*n)];
    D1 = (2/h)*sparse(D1);