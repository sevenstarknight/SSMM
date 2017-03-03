function [errorArray, resArray] = SpaamAnalysis(structMCReduct, grpSource, structPatternArray_TrainingCV)
%  tic

resArray = 0.03:0.005:0.12;
% resArray = 5:2:30;
errorArray = zeros(length(resArray),2);

%% Cycle over resolutions, determine optimal resolution

for indexi = 1:1:length(resArray) 
    %======================================================
    %% construct state space on the training set, for given resolution.
    
    states = -2:resArray(indexi):2;
    
%     logSpace = logspace(log10(0.1),log10(3),resArray(indexi));
%     states = horzcat(-fliplr(logSpace), logSpace);
    
%     spaceTotal = zeros(length(structMCReduct) ,length(states)*length(states) + 2 + 10);
    spaceTotal = zeros(length(structMCReduct) ,length(states)*length(states) + 2 );
    for i = 1:1:length(structMCReduct)
        [markovChain, ~, meanAmplitude, stdAmplitude] = ConstructMCAmp_Kernel(structMCReduct(i).timeSeries, states);

%         [fltSet] = ComputeTextureMetrics(markovChain);
        
        structMCReduct(i).MC = markovChain;
        structMCReduct(i).unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
%         spaceTotal(i,:) = [structMCReduct(i).unpackMC', meanAmplitude, stdAmplitude, fltSet];
        spaceTotal(i,:) = [structMCReduct(i).unpackMC',meanAmplitude, stdAmplitude];
    end

    %% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
    ECVAm = ecva('model', spaceTotal, grpSource, 10, 'none', 'syst123', 5, [], 0);
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

    errorArray(indexi,1) = resArray(indexi);
    errorArray(indexi,2) = errorProbMean;
end





end

