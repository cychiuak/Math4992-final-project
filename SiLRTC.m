function [X, difference_S] = SiLRTC(T,Omega,alpha,beta,maxIteration,epsilon)
X = T;
X(logical(1-Omega)) = mean(T(Omega));
%[X] = equatematrix(X,T,Omega);
difference_S = zeros(maxIteration,1);
n = ndims(X);
M = cell(n,1);
tau = alpha./beta;
betasum = sum(beta);
normT = norm(T(:));
for k = 1:maxIteration
    theX = X;
    Xsoln = 0;
    if mod(k, 20) == 0
        fprintf('SiLRTC: iterations = %d   difference=%f\n', k, difference_S(k-1));
    end
    for i = 1:n
        M{i} =  Fold(shrinkage(Unfold(theX, size(theX), i), tau(i)), size(theX), i);
        Xsoln = Xsoln + beta(i) * M{i};
    end
            Xlast = X;
    X = Xsoln/betasum;
        X(Omega) = T(Omega);
        difference_S(k) = norm(X(:)-Xlast(:)) / normT;
    if (difference_S(k) < epsilon)
        difference_S = difference_S(1:k);
        break;
    end
end
fprintf('SiLRTC ends: total iterations = %d   difference=%f\n\n', k, difference_S(k));
end

