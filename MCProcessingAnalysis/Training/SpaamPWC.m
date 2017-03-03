function [confusionMatrix, fltSpread] = SpaamPWC(fltPatternArrayMerge, grpSource_TrainCross,...
    fltPatternArray_Testing, grpSource_Testing, structPatternArray_TrainingCross)

%% ===================================================================
%PWC

spreadArray = 10.^(-2:0.2:2);
errorArray = zeros(1,length(spreadArray));

for i = 1:1:length(spreadArray)
    errorProbMean = 0;
    
    for j = 1:1:length(structPatternArray_TrainingCross)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArrayMerge, grpSource_TrainCross, structPatternArray_TrainingCross, j);

        [fltResponse, classEstimate] = Parzen_Window_Classifier(...
            fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, ...
            3, spreadArray(i), 'Missed');

        count = 0;
        for k = 1:1:length(classEstimate)
            if(strcmp(classEstimate{k}, grpSource_CrossVal{k}))
                count = count + 1;
            end
        end
        errorProb = (length(classEstimate) - count)/length(classEstimate);


        errorProbMean = errorProbMean + errorProb/5;
    end
    
    errorArray(i) = errorProbMean;
end

figure()
semilogx(spreadArray, errorArray, '.-r')
xlabel('PWC Gaussian Kernel Spread')
ylabel('Misclassification Rate 5-Fold Cross-validation')
grid on

[fltResponse, classEstimate] = Parzen_Window_Classifier(...
    fltPatternArrayMerge, grpSource_TrainCross', fltPatternArray_Testing, ...
    3, 0.4, 'Missed');

fltSpread = 0.4;

%% TPR 

tp = 0;

for i = 1:1:length(classEstimate)
    if(strcmp(classEstimate{i}, grpSource_Testing{i}))
        tp = tp + 1;
    end
end
     
tpr = tp/length(classEstimate);


%% Confusion Matrix Generation
uniqueClasses = unique(grpSource_TrainCross);

missedVector = zeros(1, length(uniqueClasses));
confusionMatrix = zeros(length(uniqueClasses), length(uniqueClasses));

for i = 1:1:length(classEstimate)
    if(strcmp(classEstimate{i}, 'Missed'))
        indexi = find(strcmp(grpSource_Testing{i}, uniqueClasses));
        missedVector(indexi) = missedVector(indexi) + 1;
    else
        indexj = find(strcmp(classEstimate{i}, uniqueClasses));
        indexi = find(strcmp(grpSource_Testing{i}, uniqueClasses));
        confusionMatrix(indexi, indexj) = confusionMatrix(indexi, indexj) + 1;
    end
end

rowSums = sum(confusionMatrix,2);
for i = 1:1:length(rowSums)
    confusionMatrix(i,:) = confusionMatrix(i,:)./rowSums(i);
end



end