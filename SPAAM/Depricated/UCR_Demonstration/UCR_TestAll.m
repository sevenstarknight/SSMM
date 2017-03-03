clc
clear

stringDirectory = 'C:\Users\kyle.johnston\Desktop\Dissertation\starvars\matlab\data\dataset';
listing = dir(stringDirectory);

errorProb_NN = length(listing);
errorProb_QDA = length(listing);

resArray = 0.01:0.01:0.03;
kernelArray = 10.^(2:0.25:4);

for i = 3:1:length(listing)

    i
    
    strLabel = listing(i).name;
    %% ====================================================
    % Make feature space and load data
    [fltTrainSet, grpSource_TrainCross, fltPatternArray_Testing, grpSource_Testing] = ReadInTestTimeSeries_UCR(strLabel);

    %% ====================================================
    % Randomly distribute into five sets
    p = randperm(length(grpSource_TrainCross));
    sampleSize = round(length(p)/5) - 1;

    structPatternArray_TrainingCV = [];
    count = 1;
    for j = 1:1:5
        indexSet = p(count:count + sampleSize - 1);
        count = count + sampleSize;
        structPatternArray_TrainingCV(j).indexSet = indexSet;
    end


    %% ====================================================
    % Cycle over resolutions, determine optimal resolution
    errorArray = zeros(length(kernelArray),length(resArray));
    
    for indexj = 1:1:length(kernelArray)
        for indexi = 1:1:length(resArray) 
            % construct state space, for given resolution.
            states = -2:resArray(indexi):2;

            spaceTotal = zeros(length(grpSource_TrainCross) ,length(states)*length(states) + 2);
            for indexk = 1:1:length(grpSource_TrainCross)
                [markovChain, meanAmplitude, stdAmplitude] = ConstructMCAmp_Kernel(fltTrainSet(indexk,:),...
                    states, 0.05, 4.0, kernelArray(indexj));

                unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
                spaceTotal(indexk,:) = [unpackMC', meanAmplitude, stdAmplitude];
            end

            %% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
            ECVAm = ecva('model', spaceTotal, grpSource_TrainCross, 60, 'none', 'syst123', 6, [], 0);
            fltPatternArray = ECVAm.CanonicalVariates;

            % CV and Error Estimates 
            errorProbMean = 0;
            for j = 1:1:length(structPatternArray_TrainingCV)

                [fltPatternArray_Training, grpSource_Training, ...
                    fltPatternArray_CrossVal, grpSource_CrossVal] = ...
                    PullTrainingAndCrossFromStruct(fltPatternArray, grpSource_TrainCross, structPatternArray_TrainingCV, j);

                [errorProb, classEstimate, structParaEst] = LDA_QDA_Classifier(...
                    fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, grpSource_CrossVal', 4); 

                errorProbMean = errorProbMean + errorProb/5;
            end

            errorArray(indexj,indexi) = errorProbMean;

        end
    end

    %% ======================================================================
    % Generate Feature Space, Find Optimal State Resolution
    [~,iResolution] = min(min(errorArray));
    [~,iKernel] = min(errorArray);
    
    states = -2:resArray(iResolution):2;
    spaceTotal = zeros(length(grpSource_TrainCross) ,length(states)*length(states));
    for indexk = 1:1:length(grpSource_TrainCross)
        [markovChain, meanAmplitude, stdAmplitude] = ConstructMCAmp_Kernel(fltTrainSet(indexk,:),...
                states, 0.05, 4.0, kernelArray(iKernel(iResolution)));

        unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
        spaceTotal(indexk,:) = [unpackMC', meanAmplitude, stdAmplitude];
    end


    fltPatternArray_TrainCross = spaceTotal; 
    grpSource_TrainCross = grpSource(idxTrainCross);

    %% =====================================
    % Extended CVA Applied to Training and Crossval 
    ECVAm = ecva('model', fltPatternArray_TrainCross, grpSource_TrainCross, 20, 'none', 'syst123', 5);
    fltPatternArrayMC = ECVAm.CanonicalVariates;
   
    mx = ECVAm.Detail.mx; % To be used for predictions
    canonicalWeights = ECVAm.CanonicalWeights; % To be used for predictions

    fltPatternArrayMerge = fltPatternArrayMC;

    % standardize the merged dataset
    vectorMean = mean(fltPatternArrayMerge,1);
    vectorStd = std(fltPatternArrayMerge, 1);
    for indexk = 1:1:length(fltPatternArrayMerge)
        fltPatternArrayMerge(indexk,:) = (fltPatternArrayMerge(indexk,:) - vectorMean)./vectorStd;
    end

    %% ==================================================================
    % transform testing data based on optimization
    for indexk = 1:1:length(grpSource_Testing)
        ECVATesting(indexk, :) = (fltPatternArray_Testing(indexk, :) - mx);
    end
    fltPatternArray_Testing_Canonical = ECVATesting*canonicalWeights;
    fltPatternArray_Testing = fltPatternArray_Testing_Canonical;

    for indexk = 1:1:length(fltPatternArray_Testing)
        fltPatternArray_Testing(indexk,:) = (fltPatternArray_Testing(indexk,:) - vectorMean)./vectorStd;
    end
    grpSource_Testing = grpSource(indexTesting);

    %% ===================================================================
    % 1-NN Testing
    [fltResponse, classEstimate] = Naive_K_Nearest(...
            fltPatternArrayMerge, grpSource_TrainCross', fltPatternArray_Testing, ...
            1.0, 2, 1.0,  'Missed');
    
    count = 0;
    for k = 1:1:length(classEstimate)
        if(strcmp(classEstimate{k}, grpSource_Testing{k}))
            count = count + 1;
        end
    end
    errorProb_NN(i) = (length(classEstimate) - count)/length(classEstimate);
        
    
        %% ===================================================================
    % QDA Testing        
    [errorProb, classEstimate, structParaEst] = LDA_QDA_Classifier(...
                fltPatternArrayMerge, grpSource_TrainCross', fltPatternArray_Testing, grpSource_Testing', 4); 
    
    errorProb_QDA(i) = errorProb;
 
end