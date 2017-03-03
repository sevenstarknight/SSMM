function [errorProb, classEstimate, structParaEst] = LDA_QDA_Classifier(...
    fltPatternArray_Training, grpSource_Training, fltPatternArray_CrossVal, grpSource_CrossVal, switch_expr)
% function [errorProb, classEstimate, structParaEst] = LDA_QDA_Classifier(...
%     fltPatternArray_Training, grpSource_Training, fltPatternArray_CrossVal, grpSource_CrossVal, ...
%      switch_expr)
%
% Author: Kyle Johnston     14 Jan 2010
%
% Usage: This function generates class estimates based on the LDA/QDA
%   algorithm. This function acts as a "traffic guard"
% 
% Input: 
%       fltPatternArray_Training: the training dataset
%       grpSource_Training: the labels for the training dataset
%       fltPatternArray_CrossVal: used to estimate error
%       grpSource_CrossVal: used to estimate error
%       switch_expr: select the LDA/QDA flavor of your choice

% Output:
%       errorProb: estimated probabiltiy of misclassification
%       classEstimate: estimated labels of fltPatternArray_CrossVal
%       structParaEst: paramters estimated by the LDA/QDA routine
%% LDA QDA Classification Routines
[structParaEst] = Train_LDA_QDA_Classification(fltPatternArray_Training, grpSource_Training, switch_expr);

[classEstimate] = Use_LDA_QDA_Classifer(fltPatternArray_CrossVal, structParaEst);

TF = strcmp(classEstimate,grpSource_CrossVal);
errorProb = 1 - sum(TF)/length(classEstimate);

end

%% Traffic Cop for LDA/QDA Classification
function [structParaEst] = Train_LDA_QDA_Classification(fltPatternArray_Training, grpSource_Training, switch_expr)

    if(switch_expr == 1)
        structParaEst = GenerateQDAParams(fltPatternArray_Training, grpSource_Training, 1); %General QDA
    elseif(switch_expr == 2)
        structParaEst = GenerateQDAParams(fltPatternArray_Training, grpSource_Training, 2); %Naive QDA
    elseif(switch_expr == 3)
        structParaEst = GenerateQDAParams(fltPatternArray_Training, grpSource_Training, 3); %Iso QDA
    elseif(switch_expr == 4)
        structParaEst = GenerateLDAParams(fltPatternArray_Training, grpSource_Training, 1); %General QDA
    elseif(switch_expr == 5)
        structParaEst = GenerateLDAParams(fltPatternArray_Training, grpSource_Training, 2); %Naive QDA
    elseif(switch_expr == 6)
        structParaEst = GenerateLDAParams(fltPatternArray_Training, grpSource_Training, 3); %Iso QDA
    else
        structParaEst = struct([]);
    end

end