function[confusionMatrix, structLRC, nodes, spread, tpr] = SpaamRBFNN(fltPatternArrayMerge, grpSource_TrainCross,...
    fltPatternArray_Testing, grpSource_Testing, structPatternArray_TrainingCross)

%% ==================================================
% RBF-NN

spreadArray = 10.^(-1:0.2:1);
errorArray = zeros(1,length(spreadArray));

for i = 1:1:length(spreadArray)
    errorProbMean = 0;
    
    for j = 1:1:length(structPatternArray_TrainingCross)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArrayMerge(:,1:5),...
            grpSource_TrainCross, structPatternArray_TrainingCross, j);

        [structLRC, nodes, spread, errorEst, classEstimate] = ...
            RadialBasisFunctionClassifier(fltPatternArray_Training, ...
            grpSource_Training, fltPatternArray_CrossVal, grpSource_CrossVal, ...
            spreadArray(i));

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
xlabel('RBF-NN Gaussian Kernel Spread')
ylabel('Misclassification Rate 5-Fold Cross-validation')
grid on

[structLRC, nodes, spread, errorEst, classEstimate] = ...
            RadialBasisFunctionClassifier(fltPatternArrayMerge, grpSource_TrainCross',... 
        fltPatternArray_Testing, grpSource_Testing', ...
            1);
        
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
