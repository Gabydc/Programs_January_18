%...  The Matmol group (2016)
     function D0 = colorthD0(n)
%...    
%... Computation of the finite element matrix D0 using the 
%... orthogonal collocation technique
    a = sqrt(3);
    b = [3*a+4 a+1 3*a-4 1-a ; 3*a-4 a-1 3*a+4 -a-1]/(6*a); 
    base(1:2,:) = horzcat(b, zeros(2,2*n-4));
    for i = 1:n-2
        base(2*i+1:2*i+2,:) = circshift(base(2*i-1:2*i,:),[0 2]);
    end

    D0 = [zeros(1,2*n) ; base ; zeros(1,2*n)];
    D0 = sparse(D0);