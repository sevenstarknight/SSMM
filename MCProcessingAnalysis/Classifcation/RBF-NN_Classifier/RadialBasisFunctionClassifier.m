function [structLRC, nodes, spread, errorEst, classEstimate] = ...
    RadialBasisFunctionClassifier(fltPatternArray_Training, ...
    grpSource_Training, fltPatternArray_CrossVal, grpSource_CrossVal, ...
    spread)

%% Initialize Parameters
intLengthTraining = length(fltPatternArray_Training(:,1));

% Make Nodes in Hyper-Dimensions
[nodes] = MakeNodes(fltPatternArray_Training);
fltRejectionValue = 0;

%% RBF Method
% Training
fltRBData = zeros(intLengthTraining, length(nodes));

for i = 1:1:intLengthTraining
    fltPositionVector = fltPatternArray_Training(i,:);
    
    for j = 1:1:length(nodes)
        centerEstimate = nodes{j};
        distance = fltPositionVector - centerEstimate;
        fDistance = exp(-spread*norm(distance));
    
%     [fDistance] = GaussianKernel(fltPositionVector, nodes, spread);
        fltRBData(i, j) = fDistance;
    end 
end

method = 2;

[structLRC, uniqueGrps] = LogisticRegressionClassifier(fltRBData, grpSource_Training, method);

[errorEst, percentReject, classEstimate] = RBFCrossVal(fltPatternArray_CrossVal, grpSource_CrossVal, structLRC, uniqueGrps , nodes, spread, fltRejectionValue);


end

function [nodes] = MakeNodes(fltPatternArray_Training)

intLengthTraining = length(fltPatternArray_Training(:,1));
nodes = cell(1,intLengthTraining);

for i = 1:1:intLengthTraining
    nodes{i} = fltPatternArray_Training(i,:);
end

end
