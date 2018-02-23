%... The MatMol Group (2016)
    function out = Jpattern(n,ne,flindex)
%...
    global D1 D2
%...
    out = zeros(n*ne,n*ne);
%...
    if flindex == 0
        out(1:n,1:n) = eye(n) + spones(D1+D2);
        out(1:n,n+1:2*n) = eye(n);
        out(n+1:2*n,1:n) = eye(n);
        out(n+1:2*n,n+1:2*n) = eye(n) + spones(D1+D2);
    else
        out(1:n,1:n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2) + spones(D2);
        out(1:n,n+1:2*n) = eye(n);
        out(n+1:2*n,1:n) = eye(n);
        out(n+1:2*n,n+1:2*n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-2),2) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),-2) + spones(D2);
    end
    out = spones(out);
