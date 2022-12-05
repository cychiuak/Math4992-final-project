function [X,difference_F] = FaLRTC(T,Omega,alpha,mu,L,C,maxIteration,epsilon)
X= T;
X(logical(1-Omega)) = mean(T(Omega));
W = X;
Z = X;
B = 0;
n = ndims(T);
dim = size(T);
Gradient = 0;
normT = norm(T(:));
difference_F = zeros(maxIteration,1);
LP = L;
mu0 = mu;
for k = 1:maxIteration
    if mod(k, 20) == 0
         fprintf('FaLRTC: iterations = %d   difference=%f\n', k, difference_F(k-1));
    end
    gradcoef = alpha.^2 ./ mu;
    Xprev = X;
    while true
        updated = false;
        theta = L*(1+sqrt(1+4*L*B)) / (2*L);
        W = (theta/L)/(B+theta/L) * Z + B/(B+theta/L) * Xprev;
        Gradient = Gradient * 0;
        fmiu = 0;
        fmiuX = 0;
        for i = 1:n
            [temp,sigma] = Truncate(Unfold(W, dim, i), 1);
            temp = Fold(temp, dim, i);
            Gradient = Gradient + gradcoef(i)*temp;
            fmiu = fmiu + gradcoef(i)*(sum(sigma)-n);
        end
        for i = 1:n
            [temp,sigma] = Truncate(Unfold(X, dim, i), 1);
            fmiuX = fmiuX + gradcoef(i)*(sum(sigma)-n);
        end
        
        Gradient(Omega) = 0;
        X = W - Gradient/L;
        fmiux = 0;
        for i = 1 : n
            [sigma] = SingularValue(Unfold(X, dim, i));
            fmiux = fmiux + gradcoef(i)*(sum(sigma.^2) - n);
        end
        if(fmiu-fmiux)*L<sum(Gradient(:).^2)/2
            break;
    
        else
            L = L/C;
            break;
        end
    end
        difference_F(k) =  norm(X(:)-Xprev(:)) / normT;
        if difference_F(k) < epsilon
            break;
        end
    Z = Z - (theta/L)*Gradient;
    B = B+theta/L;
end
difference_F = difference_F(1:k);
fprintf('FaLRTC ends: total iterations = %d   difference=%f\n\n', k, difference_F(k));
end

