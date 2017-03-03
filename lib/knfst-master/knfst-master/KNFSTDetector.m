function [model, nodes, errorEst, transformedScores] = ...
    KNFSTDetector(fltPatternArray_Training, fltPatternArray_CrossVal, spread)

%% Initialize Parameters & Make Nodes in Hyper-Dimensions
intLengthTraining = length(fltPatternArray_Training(:,1));

nodes = cell(1,intLengthTraining);
for i = 1:1:intLengthTraining
    nodes{i} = fltPatternArray_Training(i,:);
end


%% RBF Method
% Training
[fltRBData] = GaussianKernelKNFST(fltPatternArray_Training, nodes, spread);
model = learn_oneClassNovelty_knfst(fltRBData);


%% Cross Val
[fltRBCrossData] = GaussianKernelKNFST(fltPatternArray_CrossVal, nodes, spread);
scores = test_oneClassNovelty_knfst(model, fltRBCrossData');


transformedScores = log10(scores);
transformedScores = (transformedScores - mean(transformedScores))/std(transformedScores);

errorEst = length(transformedScores(transformedScores > 3.0))/length(transformedScores);


end

function [fltRBData] = GaussianKernelKNFST(fltPatternArray, nodes, spread)

fDistance = zeros(1,length(nodes));
intLengthPatternArray = length(fltPatternArray);
fltRBData = zeros(intLengthPatternArray, length(nodes));

for j = 1:1:intLengthPatternArray
    fltPositionVector = fltPatternArray(j,:);
    for i = 1:1:length(nodes)
        distance = norm(fltPositionVector - nodes{i});
        fDistance(i) = exp(-(distance)/spread);
    end
    
    fltRBData(j, :) = fDistance;
end

end
