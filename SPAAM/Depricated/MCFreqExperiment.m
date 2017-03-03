clc
clear

%% ====================================================
% Make feature space and load data
resArray = 0.02:0.02:0.3;
errorArray = zeros(1,length(resArray));

[structMCReduct, grpSource] = ReadInTimeSeries();

%% Add in flat lines

[structMCReduct, grpSource] = ConstructNonVariableStars(1000, structMCReduct, grpSource);

[structPatternArray_TrainingCV, indexTesting] = ...
    Generate_5FoldCrossVal(grpSource);
uniqueLabels = unique(grpSource);
tic

%% Cycle over resolutions
for index = 1:1:length(resArray) 
    % construct state space, for given resolution.
    states = -2:resArray(index):2;
    
    structFreqSet = struct('freqSet', []);
    for j = 1:1:length(uniqueLabels)
        structFreqSet(j).freqSet = zeros(length(states), length(states));
    end
    
    
    freqProbs = zeros(length(structMCReduct) ,length(states)*length(states));
    spaceTotal = zeros(length(structMCReduct) ,length(states)*length(states) + 2);
    for i = 1:1:length(structMCReduct)
        [markovChain, ~, ~] = ConstructMCAmp(structMCReduct(i).timeSeries, states);
        [transitionOccurance, meanAmplitude, stdAmplitude] = ConstructMCFreqAmp(structMCReduct(i).timeSeries, states);
        
        structMCReduct(i).MC = markovChain;
        structMCReduct(i).unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
        spaceTotal(i,:) = [structMCReduct(i).unpackMC', meanAmplitude, stdAmplitude];
        
        for j = 1:1:length(uniqueLabels)
            if(strcmp(grpSource{i}, uniqueLabels{j}))
                structFreqSet(j).freqSet = structFreqSet(j).freqSet + transitionOccurance;
            end
        end
        
    end
    
    %%
    for j = 1:1:length(uniqueLabels)
        markovChainClass = structFreqSet(j).freqSet;
        
        %% turn counts into probablities
        totalInState = sum(markovChainClass, 2);
        for(i = 1:1:length(totalInState))
            if(totalInState(i) ~= 0.0)
                markovChainClass(i,:) = (markovChainClass(i,:)./totalInState(i));
            end
        end  
        
        structFreqSet(j).markovChain = markovChainClass;
        structFreqSet(j).markovChainShape = reshape(markovChainClass, [length(states)*length(states), 1]);
    end
    
    deltaProbs = zeros(length(structMCReduct) ,length(uniqueLabels));
    for i = 1:1:length(grpSource)
        for j = 1:1:length(uniqueLabels)
           deltaProbs(i,j) = sum((structFreqSet(j).markovChainShape - structMCReduct(i).unpackMC).^2);
        end
    end
    
    spaceTotal = horzcat(spaceTotal,deltaProbs);
    
    %% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
    ECVAm = ecva('model', spaceTotal, grpSource, 20, 'none', 'syst123', 6);
    fltPatternArray = ECVAm.CanonicalVariates;
    
    %% CV and Error Estimates 
    errorProbMean = 0;
    for j = 1:1:length(structPatternArray_TrainingCV)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArray, grpSource, structPatternArray_TrainingCV, j);

        [errorProb, classEstimate, structParaEst] = LDA_QDA_Classifier(...
            fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, grpSource_CrossVal', 4); 

        errorProbMean = errorProbMean + errorProb/5;
    end

    errorArray(index) = errorProbMean;

end

toc

plot(resArray, errorArray, '.-r')
ylabel('Misclassification Rate')
xlabel('State Space Resolution')

%% ======================================================================
% Generate Feature Space JUST on Variablity
states = -2:0.08:2;
    
    structFreqSet = struct('freqSet', []);
    for j = 1:1:length(uniqueLabels)
        structFreqSet(j).freqSet = zeros(length(states), length(states));
    end
    
    
    freqProbs = zeros(length(structMCReduct) ,length(states)*length(states));
    spaceTotal = zeros(length(structMCReduct) ,length(states)*length(states));
    for i = 1:1:length(structMCReduct)
        [markovChain, ~, ~] = ConstructMCAmp(structMCReduct(i).timeSeries, states);
        [transitionOccurance, meanAmplitude, stdAmplitude] = ConstructMCFreqAmp(structMCReduct(i).timeSeries, states);
        
        structMCReduct(i).MC = markovChain;
        structMCReduct(i).unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
        spaceTotal(i,:) = structMCReduct(i).unpackMC';
        
        for j = 1:1:length(uniqueLabels)
            if(strcmp(grpSource{i}, uniqueLabels{j}))
                structFreqSet(j).freqSet = structFreqSet(j).freqSet + transitionOccurance;
            end
        end
        
    end
    
    %%
    for j = 1:1:length(uniqueLabels)
        markovChainClass = structFreqSet(j).freqSet;
        
        %% turn counts into probablities
        totalInState = sum(markovChainClass, 2);
        for(i = 1:1:length(totalInState))
            if(totalInState(i) ~= 0.0)
                markovChainClass(i,:) = (markovChainClass(i,:)./totalInState(i));
            end
        end  
        
        structFreqSet(j).markovChain = markovChainClass;
        structFreqSet(j).markovChainShape = reshape(markovChainClass, [length(states)*length(states), 1]);
    end
    
%     deltaProbs = zeros(length(structMCReduct) ,length(uniqueLabels));
%     for i = 1:1:length(grpSource)
%         for j = 1:1:length(uniqueLabels)
%            deltaProbs(i,j) = sum((structFreqSet(j).markovChainShape - structMCReduct(i).unpackMC).^2);
%         end
%     end
    
%     spaceTotal = horzcat(spaceTotal,deltaProbs);

    ECVAm = ecva('model', spaceTotal, grpSource, 20, 'none', 'syst123', 6);
    
    mx = ECVAm.Detail.mx; % To be used for predictions
    canonicalWeights = ECVAm.CanonicalWeights; % To be used for predictions

    
    % transform testing data
for i = 1:1:length(uniqueLabels)
    ECVATesting(i, :) = (structFreqSet(i).markovChainShape - mx');
end
fltPatternArray_Testing_Canonical = ECVATesting*canonicalWeights;
    



    
idxTrainCross = horzcat(structPatternArray_TrainingCV(1).indexSet, structPatternArray_TrainingCV(2).indexSet,...
    structPatternArray_TrainingCV(3).indexSet, structPatternArray_TrainingCV(4).indexSet, ...
    structPatternArray_TrainingCV(5).indexSet);

fltPatternArray_TrainCross = spaceTotal(idxTrainCross, :); 
grpSource_TrainCross = grpSource(idxTrainCross);

%% =====================================
% Extended CVA Applied to Training and Crossval (Total, with 5-Fold Validation)      
ECVAm = ecva('model', fltPatternArray_TrainCross, grpSource_TrainCross, 20, 'none', 'syst123', 5);
fltPatternArrayMC = ECVAm.CanonicalVariates;
   
mx = ECVAm.Detail.mx; % To be used for predictions
canonicalWeights = ECVAm.CanonicalWeights; % To be used for predictions
% Test for Equal Outputs
% for i = 1:1:length(grpSource_TrainCross)
%     ECVAExample(i, :) = (fltPatternArray_TrainCross(i, :) - mx);
% end
%
% canonicalVariates = ECVAExample*ECVAm.CanonicalWeights;
% realVariates = ECVAm.CanonicalVariates;

% Merge with photometric data
counter = 1;
fltPatternPhotometric = [];
for i = idxTrainCross

    fltPatternPhotometric(counter,:) = structMCReduct(i).parameters;
    counter = counter + 1;

end
fltPatternArrayMerge = horzcat(fltPatternArrayMC, fltPatternPhotometric);

% standardize the merged dataset
vectorMean = mean(fltPatternArrayMerge,1);
vectorStd = std(fltPatternArrayMerge, 1);
for i = 1:1:length(fltPatternArrayMerge)
    
    fltPatternArrayMerge(i,:) = (fltPatternArrayMerge(i,:) - vectorMean)./vectorStd;
    
end

%% Grab Testing Data
fltPatternArray_Testing = spaceTotal(indexTesting, :);
grpSource_Testing = grpSource(indexTesting);

% transform testing data
for i = 1:1:length(grpSource_Testing)
    ECVATesting(i, :) = (fltPatternArray_Testing(i, :) - mx);
end
fltPatternArray_Testing_Canonical = ECVATesting*canonicalWeights;

counter = 1;
fltPatternPhotometricTesting = [];
for i = indexTesting

    fltPatternPhotometricTesting(counter,:) = structMCReduct(i).parameters;
    counter = counter + 1;

end

fltPatternArray_Testing = horzcat(fltPatternArray_Testing_Canonical, fltPatternPhotometricTesting);

for i = 1:1:length(fltPatternArray_Testing)
    fltPatternArray_Testing(i,:) = (fltPatternArray_Testing(i,:) - vectorMean)./vectorStd;
end
grpSource_Testing = grpSource(indexTesting);

%% PWC
spreadArray = 10.^(-4:0.1:1);
errorArray = zeros(1,length(spreadArray));

for i = 1:1:length(spreadArray)
    errorProbMean = 0;
    
    for j = 1:1:length(structPatternArray_TrainingCV)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArrayMerge, grpSource_TrainCross, structPatternArray_TrainingCV, j);

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


[fltResponse, classEstimate] = Parzen_Window_Classifier(...
    fltPatternArray_Training, grpSource_Training', fltPatternArray_Testing, ...
    3, 0.125, 'Missed');

% Confusion Matrix Generation
uniqueClasses = unique(classEstimate);

confusionMatrix = zeros(length(uniqueClasses), length(uniqueClasses));

for i = 1:1:length(classEstimate)
   indexj = find(strcmp(classEstimate{i}, uniqueClasses));
   indexi = find(strcmp(grpSource_Testing{i}, uniqueClasses));
    
    confusionMatrix(indexi, indexj) = confusionMatrix(indexi, indexj) + 1;
end

rowSums = sum(confusionMatrix,2);
for i = 1:1:length(rowSums)
    confusionMatrix(i,:) = confusionMatrix(i,:)./rowSums(i);
end

