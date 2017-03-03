centroids = 4:2:26;
errorArray = zeros(size(centroids));   
 
uniqueTraining = unique(grpSource);
    
for i = 1:1:length(centroids)

%% CV and Error Estimates 
errorProbMean = 0;
for j = 1:1:length(structPatternArray_TrainingCV)

    [fltPatternArray_Training, grpSource_Training, ...
        fltPatternArray_CrossVal, grpSource_CrossVal] = ...
        PullTrainingAndCrossFromStruct(fltPatternArray, grpSource, structPatternArray_TrainingCV, j);

%         [structLRC, nodes, spread, errorEst, classEstimate] = ...
%             RadialBasisFunctionClassifier(fltPatternArray_Training, ...
%             grpSource_Training', fltPatternArray_CrossVal, grpSource_CrossVal', spread);

    idxClass = zeros(1,length(grpSource_Training));
    for(k = 1:1:length(grpSource_Training))
       idxClass(k) = find(ismember(uniqueTraining,grpSource_Training{k}));
    end

    [Centers, betas, Theta] = trainRBFN(fltPatternArray_Training, idxClass', centroids(i), false);

    counter = 0;
    for(k = 1:1:length(grpSource_CrossVal))
        z = evaluateRBFN(Centers, betas, Theta, fltPatternArray_CrossVal(k,:));
        zest = find(z == max(z));
        zprime = find(ismember(uniqueTraining,grpSource_CrossVal{k}));
        if(zest == zprime)
            counter = counter + 1;
        end
    end
    errorEst = 1 - counter/length(grpSource_CrossVal);
    errorProbMean = errorProbMean + errorEst/5;
end

errorArray(i) = errorProbMean;

end


plot(centroids, errorArray, '.-r')
xlabel('# centroids')
ylabel('Classification Error')