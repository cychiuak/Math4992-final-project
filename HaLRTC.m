function [X, difference_H] = HaLRTC(T, Omega, alpha, row, maxIter, epsilon)
X = T;
X(logical(1-Omega)) = mean(T(Omega));
difference_H = zeros(maxIter, 1);
dim = size(T);
Y = cell(ndims(T), 1);
M = Y;
normT = norm(T(:));
n = ndims(T);
for i = 1:ndims(T) 
    Y{i} = zeros(dim);
    M{i} = zeros(dim);
end
for k = 1: maxIter
    row = row * 1.15;
    

    Msum = 0;
    Ysum = 0;
    for i = 1:n
        result = X + Y{i} / row;
        M{i} = Fold(shrinkage(Unfold(result, dim, i), alpha(i)/row), dim, i);
        Msum = Msum + M{i} - Y{i}/row;
    end
    
    prevX = X;
    X = (Msum) / (ndims(T));
    X(Omega) = T(Omega);
    for i = 1:n
        Y{i} = Y{i} - row*(M{i}-X);
    end

    difference_H(k) = norm(X(:)-prevX(:)) / normT;
    if difference_H(k) < epsilon
        break;
    end
end

difference_H = difference_H(1:k);
fprintf('HaLRTC ends: total iterations = %d   difference=%f\n\n', k, difference_H(k));
end

