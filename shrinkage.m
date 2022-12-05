function [X] = shrinkage(Z, tau)
[m, n] = size(Z);
if 2*m < n
        [s,v,d] = svd(Z*Z');
        x = diag(v);
        x = sqrt(x);
        t = tau;
        x = max(x-t,0)./x;
        x = diag(x);
        X = s*x*s'*Z;
    return;
end
if m > 2*n
        [s,v,d] = svd(Z'*Z);
        x = diag(v);
        x = sqrt(x);
        t = tau;
        x = max(x-t,0)./x;
        x = diag(x);
        X = Z*s'*x*s;
    return;
end
[S,V,D] = svd(Z);
Sigma2 = diag(V).^2;
n = sum(diag(V) > tau);
X = S(:, 1:n) * max(V(1:n,1:n)-tau, 0) * D(:, 1:n)';