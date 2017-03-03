function [structLDAParam] = GenerateLDAParams(fltPatternArray_Training, grpSource_Training, switch_expr)
% function [structLDAParam] = GenerateLDAParams(fltPatternArray_Training, grpSource_Training,
% switch_expr)
%
% Author: Kyle Johnston     14 Jan 2010
%
% Usage: This function generates class estimates based on the LDA
%   algorithm. Ties are represented by the "missValue" string in the
%   output labels. 
% 
% Input: 
%       fltPatternArray_Training: the training dataset
%       grpSource_Training: the labels for the training dataset
%       crit_alpha: int representing the kernel selected

% Output:
%       structLDAParam: distribution parameters  

%% Initialize variables
uniqueGrp = unique(grpSource_Training);
intNumberGroups = length(uniqueGrp);

structLDAParam = struct([]);

% Select method for generation of covariance matrix
if(switch_expr == 1)
    [covMatrix] = GenerateGeneralCaseLDA(fltPatternArray_Training, grpSource_Training, uniqueGrp, intNumberGroups);
elseif(switch_expr == 2)
    [covMatrix] = GenerateNaiveCaseLDA(fltPatternArray_Training, grpSource_Training, uniqueGrp, intNumberGroups);
elseif(switch_expr == 3)
    [covMatrix] = GenerateIsotropicLDA(fltPatternArray_Training, grpSource_Training, uniqueGrp, intNumberGroups);
else
    covMatrix = [];
end


% Estimate and store parameters
for i = 1:1:intNumberGroups
    boolDecision =  strcmp(uniqueGrp{i},grpSource_Training);
    
    fltReducedSet = fltPatternArray_Training(boolDecision,:);
    
    structLDAParam(i).mean = mean(fltReducedSet, 1);
    structLDAParam(i).n = sum(boolDecision);
    structLDAParam(i).cov = covMatrix;
    structLDAParam(i).invCov = inv(structLDAParam(i).cov);
    structLDAParam(i).logDet = log(det(structLDAParam(i).cov));
    
    structLDAParam(i).constant = log(structLDAParam(i).n/length(grpSource_Training)) - 0.5*structLDAParam(i).logDet;
    
    structLDAParam(i).type = uniqueGrp(i);
end

end


%% Subfunction General LDA
function [covMatrix] = GenerateGeneralCaseLDA(fltPatternArray_Training, grpSource_Training, uniqueGrp, intNumberGroups)

intLengthData = length(fltPatternArray_Training(:,1));
dimen = length(fltPatternArray_Training(1,:));

covMatrix = zeros(dimen);
    
for j = 1:1:intNumberGroups
    boolDecision =  strcmp(uniqueGrp{j},grpSource_Training);
    fltReducedSet = fltPatternArray_Training(boolDecision,:);
    intLengthSet = length(fltReducedSet(:,1));
    
    covMatrix = covMatrix + cov(fltReducedSet)*intLengthSet/intLengthData;
end

end


%% Subfunction Naive LDA
function [covMatrix] = GenerateNaiveCaseLDA(fltPatternArray_Training, grpSource_Training, uniqueGrp, intNumberGroups)

intLengthData = length(fltPatternArray_Training(:,1));
dimen = length(fltPatternArray_Training(1,:));

covMatrix = zeros(dimen);
    
for j = 1:1:intNumberGroups
    boolDecision =  strcmp(uniqueGrp{j},grpSource_Training);
    fltReducedSet = fltPatternArray_Training(boolDecision,:);
    intLengthSet = length(fltReducedSet(:,1));
    
    for k = 1:1:dimen
        covMatrix(k,k) = covMatrix(k,k) + var(fltReducedSet(:,k))*intLengthSet/intLengthData;
    end
end

end

%% Subfunction Isotropic LDA
function [covMatrix] = GenerateIsotropicLDA(fltPatternArray_Training, grpSource_Training, uniqueGrp, intNumberGroups)

    intLengthData = length(fltPatternArray_Training(:,1));
    dimen = length(fltPatternArray_Training(1,:));

    varMatrix = zeros(dimen);

    for j = 1:1:intNumberGroups
        boolDecision =  strcmp(uniqueGrp{j},grpSource_Training);
        fltReducedSet = fltPatternArray_Training(boolDecision,:);
        intLengthSet = length(fltReducedSet(:,1));

        for k = 1:1:dimen
            varMatrix(k,k) = varMatrix(k,k) + var(fltReducedSet(:,k))*intLengthSet/intLengthData;
        end
    end

    varPool = 0;
    
    for j = 1:1:dimen
        varPool = varMatrix(j,j) + varPool;
    end
    
    covMatrix = eye(dimen)*(varPool/dimen);
end