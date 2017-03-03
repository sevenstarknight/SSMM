function [confusionMatrix, fltSpread, tpr] = SpaamKNN(fltPatternArrayMerge, grpSource_TrainCross,...
    fltPatternArray_Testing, grpSource_Testing, structPatternArray_TrainingCross)

%% ===================================================================
%kNN
fltSpread = 1.0;
kNeigh = 1:1:10;
errorArray = zeros(1,length(kNeigh));

for i = 1:1:length(kNeigh)
    errorProbMean = 0;
    
    for j = 1:1:length(structPatternArray_TrainingCross)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArrayMerge, grpSource_TrainCross, structPatternArray_TrainingCross, j);
        
        [fltResponse, classEstimate] = Naive_K_Nearest(...
            fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, ...
            kNeigh(i), 2, 1.0,  'Missed');

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
plot(kNeigh, errorArray, '.-r')
xlabel('kNN Number of k')
ylabel('Misclassification Rate 5-Fold Cross-validation')
grid on


[fltResponse, classEstimate] = Naive_K_Nearest(...
            fltPatternArrayMerge, grpSource_TrainCross', fltPatternArray_Testing, ...
            1, 2, 1.0,  'Missed');

%% Misclassification Rate

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