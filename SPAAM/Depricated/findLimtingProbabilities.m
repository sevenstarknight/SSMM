function[est] = findLimtingProbabilities(A)

X = A;
[n,m] = size(X);

for(i = 1:1:n)
    
    X(i,i) = A(i,i) - 1;
    
end

X(n + 1,:) = ones(1,m);
Y = [zeros(1,n), 1];

est = (X'*X)^(-1)*X'*Y';

end