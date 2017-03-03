function [structLRC, nodes, spread, errorEst, classEstimate] = RadialBasisFunctionClassifier_KMeans(fltPatternArray_Training, grpSource_Training, fltPatternArray_CrossVal, grpSource_CrossVal, spread)

%% Initialize Parameters
intLengthTraining = length(fltPatternArray_Training(:,1));
uniqueGrp = unique(grpSource_Training);

structPsi_All = [];

for i = 1:1:length(uniqueGrp)
    TF = strcmp(grpSource_Training, uniqueGrp(i));
    
    fltLimitedSet = fltPatternArray_Training(TF,:);
    
    fltUniquePatterns = unique(fltLimitedSet, 'rows');
    
    if(length(fltUniquePatterns(:,1)) > 5)
        [structPsi] = K_Means(fltPatternArray_Training(TF,:)', 5);
    else
        [structPsi] = K_Means(fltPatternArray_Training(TF,:)', length(fltUniquePatterns(:,1)));
    end
    
    structPsi_All(i).structPsi = structPsi;
    
end

% Make Nodes in Hyper-Dimensions
[nodes] = MakeNodes(structPsi_All);
fltRejectionValue = 0;

%% RBF Method
% Translate into Radial Basis Function

% Training
fltRBData = zeros(intLengthTraining, length(nodes));

for i = 1:1:intLengthTraining
    fltPositionVector = fltPatternArray_Training(i,:);
    [fDistance] = GaussianKernel(fltPositionVector, nodes, spread);
    fltRBData(i, :) = fDistance;
end

method = 2;

[structLRC, uniqueGrps] = LogisticRegressionClassifier(fltRBData, grpSource_Training, method);

[errorEst, percentReject, classEstimate] = RBFCrossVal(fltPatternArray_CrossVal, grpSource_CrossVal, structLRC, uniqueGrps , nodes, spread, fltRejectionValue);

end

function [nodes] = MakeNodes(structPsi_All)

intClasses = length(structPsi_All);
count = 0;

for i = 1:1:intClasses
    intNodes = length(structPsi_All(i).structPsi);
    structPsi = structPsi_All(i).structPsi;
    for j = 1:1:intNodes
        count = count + 1;
        nodes{count} = structPsi(j).mean;
    end
end


end