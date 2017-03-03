clc
clear

%% ====================================================
% Make feature space and load data
resArray = 0.02:0.01:0.3;
errorArray = zeros(1,length(resArray));

[structMCReduct, grpSource] = ReadInTimeSeries();

%% Add in flat lines

[structMCReduct, grpSource] = ConstructNonVariableStars(1000, structMCReduct, grpSource);

[structPatternArray_TrainingCV, indexTesting] = ...
    Generate_5FoldCrossVal(grpSource); 

% construct state space, for given resolution.
    states = -2:resArray(7):2;
    
    spaceTotal = zeros(length(structMCReduct) ,length(states)*length(states));
    limitTotal = zeros(length(structMCReduct) ,length(states));
    
    for i = 1:1:length(structMCReduct)
        [markovChain] = ConstructMCDiff(structMCReduct(i).timeSeries, states);
        structMCReduct(i).MC = markovChain;
        
        
        [est] = findLimtingProbabilities(markovChain);
        
        structMCReduct(i).limitingProb = est;
        structMCReduct(i).unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
        spaceTotal(i,:) = structMCReduct(i).unpackMC';
        limitTotal(i,:) = est;
    end

    %% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
    ECVAm = ecva('model', spaceTotal, grpSource, 20, 'none', 'syst123', 6);
    
    fltPatternArray = ECVAm.CanonicalVariates;
    
    
    %% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
    ECVAlimit = ecva('model', limitTotal, grpSource, 20, 'none', 'syst123', 6);
    
    spreadArray = 10.^(-3.25:0.25:-3.25);
    
    uniqueTraining = unique(grpSource);
    
    for i = 1:1:length(spreadArray)
    
    %% CV and Error Estimates 
    errorProbMean = 0;
    spread = spreadArray(i);
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
        
        [Centers, betas, Theta] = trainRBFN(fltPatternArray_Training, idxClass, false);
        
        counter = 0;
        for(k = 1:1:length(grpSource_CrossVal))
            z = evaluateRBFN(Centers, betas, Theta, fltPatternArray_CrossVal(k,:));
        
            zprime = find(grpSource_CrossVal{k} == uniqueTraining);
            if(z == zprime)
                counter = counter + 1;
            end
        end
        errorEst = 1 - counter/length(grpSource_CrossVal);
        errorProbMean = errorProbMean + errorEst/5;
    end

    errorArray(i) = errorProbMean;
    
    end