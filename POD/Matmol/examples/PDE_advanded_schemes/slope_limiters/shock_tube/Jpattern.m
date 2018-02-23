%... The MatMol Group (2016)
    function Jp = Jpattern(n,ne)
%...
    global D2
%...
    Jp = zeros(n*ne,n*ne);
    Jp(1:n,1:n) = eye(n) + spones(D2);
    Jp(1:n,n+1:2*n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2);
    Jp(1:n,2*n+1:3*n) = eye(n);
    Jp(n+1:2*n,1:n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2);
    Jp(n+1:2*n,n+1:2*n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2) + spones(D2);
    Jp(n+1:2*n,2*n+1:3*n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2);
    Jp(2*n+1:3*n,1:n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2);
    Jp(2*n+1:3*n,n+1:2*n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2);
    Jp(2*n+1:3*n,2*n+1:3*n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2) + spones(D2);
    Jp = spones(Jp);
    Jp = sparse(Jp);
