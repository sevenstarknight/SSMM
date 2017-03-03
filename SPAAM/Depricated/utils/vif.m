function [vifArray] = vif(X)

aX = (X{1})';

[n,m] = size(aX);
vifArray = zeros(1,m);

indexSet = 1:1:m;

for(i = 1:1:m)

    idxX = indexSet(indexSet ~= i);
    
    yTmp = aX(:,i);
    xTemp = horzcat(ones(n, 1), aX(:,idxX));
    
    [b,bint,r,rint,stats]  = regress(yTmp,xTemp);
    
    vifArray(i) = 1/(1 - stats(1));

end


S = sparse(A)

end