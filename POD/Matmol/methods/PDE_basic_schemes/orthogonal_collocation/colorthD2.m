%...  The Matmol group (2016)
     function D2 = colorthD2(n,h)
%...    
%... Computation of the finite element matrix D2 with the
%... orthogonal collocation technique
    a = sqrt(3);
    b = [-a -a-1 a 1-a ; a a-1 -a 1+a]/2; 
    base(1:2,:) = horzcat(b, zeros(2,2*n-4));
    for i = 1:n-2
        base(2*i+1:2*i+2,:) = circshift(base(2*i-1:2*i,:),[0 2]);
    end

    D2 = [zeros(1,2*n) ; base ; zeros(1,2*n)];
    D2 = (4/(h^2))*sparse(D2);