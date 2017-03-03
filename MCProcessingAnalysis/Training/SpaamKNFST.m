function [model, nodes, errorEst] = ...
    SpaamKNFST(fltPatternArrayMerge, grpSource_TrainCross, fltPatternArray_Testing, structPatternArray_TrainingCross)

spreadArray = 10.^(0:0.05:1);
errorArray = zeros(1,length(spreadArray));

for i = 1:1:length(spreadArray)
    errorProbMean = 0;
    
    for j = 1:1:length(structPatternArray_TrainingCross)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArrayMerge,...
            grpSource_TrainCross, structPatternArray_TrainingCross, j);

        [model, nodes, errorEst, scores] = ...
            KNFSTDetector(fltPatternArray_Training, fltPatternArray_CrossVal, spreadArray(i));


        errorProbMean = errorProbMean + errorEst/5;
    end
    
    errorArray(i) = errorProbMean;
end

semilogx(spreadArray, errorArray, '.-r')
xlabel('KNFST Kernel Width')
ylabel('Misclassification Rate 5-Fold Cross-validation')
grid on


[model, nodes, errorEst, transformedScores] = ...
    KNFSTDetector(fltPatternArray_Training, fltPatternArray_CrossVal, 1.0);

[fltRBTestingData] = GaussianKernelKNFST(fltPatternArray_Testing, nodes, 1.0);
scores = test_oneClassNovelty_knfst(model, fltRBTestingData');

transformedScores = log10(scores);
transformedScores = (transformedScores - mean(transformedScores))/std(transformedScores);

errorEst = length(transformedScores(transformedScores > 3.0))/length(transformedScores);


end