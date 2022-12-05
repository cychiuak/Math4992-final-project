function [X] = equatematrix(X,T,Omega)
dim = size(X);
x1 = Unfold(X,dim,1);
t1 = Unfold(T,dim,1);
o1 = Unfold(Omega,dim,1);
[m,n] = size(x1);
for i = 1:m
    for j = 1:n
        if o1(m,n) == 1
            x1(m,n) = t1(m,n);
        end
    end
end
X = Fold(x1,dim,1);
end

