function [structLRC, uniqueGrps] = LogisticRegressionClassifier(fltPatternArray_Training, grpSource_Training, method)
%% Initialize Parameters
uniqueGrps = unique(grpSource_Training);
intUniqueGroups = length(uniqueGrps);

[fltNewTrainingData] = InitializeTrainingData(fltPatternArray_Training);
intDimensions = length(fltNewTrainingData(1,:));
intLengthTraining = length(fltNewTrainingData(:,1));

fltThetaParams = zeros(intUniqueGroups, intDimensions) + 0.5;

structLRC = struct(...
    'fltNewTrainingData', fltNewTrainingData, ...
    'intDimensions', intDimensions, ...
    'intLengthTraining', intLengthTraining, ...
    'intUniqueGroups', intUniqueGroups, ...
    'fltThetaParams', fltThetaParams);

%% LRC Method
if(method == 1)
    [structLRC] = ComputeGradDescent(structLRC, grpSource_Training);
else
    [structLRC] = ComputeNewtonIteration(structLRC, grpSource_Training);
end

end