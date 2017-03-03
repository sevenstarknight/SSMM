function [Xmean,meanX] = sub_mc(X)
n = size(X,1);
meanX = mean(X);
Xmean = (X-meanX(ones(n,1),:));

end