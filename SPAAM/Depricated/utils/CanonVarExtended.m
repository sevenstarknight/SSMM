function [t_cv] = CanonVarExtended(X,K)
%CANONVAR  Canonical Variates (no post rotation of variates)
% W=CanonVar(X,K,<shrinkage>)
%
% Inputs:
% X : cell array of data matrices in which each column contains a
% datapoint. X{c} contains the data for class c
% K : dimension of the projection
%
% Output:
% W : the eigenvector (or approximation)
% 
% see demoCanonVarDigits.m

%% Initialize
nClasses = length(X); 
D = size(X{1},1); % dimension of data

%% For each class
XX=[];
n = zeros(1,nClasses);
x_bar = zeros(D, 1);
x_bar_i = zeros(D,nClasses);
for c = 1: nClasses
   % n_i
   n(c) = size(X{c},2);
   
   %x_bar_i
   x_bar_i(:,c) = mean(X{c},2);
   
   % sum(n*x_bar_i)
   x_bar = x_bar + n(c)*x_bar_i(:,c);
   
   % total population
   XX=[XX X{c}];
end
x_bar = x_bar/sum(n);

% basis to represent the solution
Q=orth(XX);

[sWithin, sBetween, Y] = scatterEstimate(X, n, x_bar_i, x_bar);

[~,p] = chol(sWithin);
if(p == 0)
    % [V,D]
    [w, ~]=eigs(sWithin\sBetween,K);
else
    ncomp = min(nClasses - 1, D);
    [~,~,~,~,w] = plsregress(sBetween,Y,ncomp);
end

t_cv = Q*(w(:,1:nClasses - 1));

end

function[sWithin, sBetween, Y] = scatterEstimate(X, n, x_bar_i, x_bar)

nClasses = length(X); 
sWithin = zeros(size(cov((X{1})')));
sBetween = zeros(size(cov((X{1})')));

for c = 1:1:nClasses
   
    %% Within Scatter
    for i = 1:1:n(c)
        halfCov = X{c}(:,i) - x_bar_i(:,c);
        sWithin = sWithin + halfCov*halfCov';
    end
    
    %% Between Scatter
    delta = (x_bar_i(c) - x_bar);
    sBetween = sBetween + n(c)*(delta*delta');
    
    Y(:,c) = delta;
end

sWithin = sWithin/(sum(n) - nClasses);
sBetween = sBetween/(nClasses - 1);

end