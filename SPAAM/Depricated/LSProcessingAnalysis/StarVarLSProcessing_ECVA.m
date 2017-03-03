clc
clear

%% ====================================================
% Make feature space and load data
errorArray = zeros(1,length(resArray));

[structLSReduct, grpSource] = ReadInTimeSeries();

[structPatternArray_TrainingCV, indexTesting] = ...
    Generate_5FoldCrossVal(grpSource);

for index = 1:1:length(structLSReduct)
   
    meanDeltaT(index) = mean(diff(structLSReduct(index).timeSeries(:,1)));
    minDeltaT(index) = min(diff(structLSReduct(index).timeSeries(:,1)));
    maxDeltaT(index) = max(diff(structLSReduct(index).timeSeries(:,1)));
    sizeSample(index) = length(structLSReduct(index).timeSeries(:,1));
end

median(sizeSample)
fmax = 1/(2*median(meanDeltaT));
fmin = 1/(4*median(sizeSample)*median(meanDeltaT));

%% Cycle over resolutions
% construct frequencies to be evaluated, for given resolution.
freqSet = fmin:fmin:fmax;

spaceTotal = [];
for i = 1:1:length(structLSReduct)
    [lombVector,~] = plomb(structLSReduct(i).timeSeries(:,2), structLSReduct(i).timeSeries(:,1),freqSet);
    
    structLSReduct(i).lombVector = lombVector;
    spaceTotal(i,:) = [lombVector; diff(lombVector)];
end


%% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
ECVAm = ecva('model', spaceTotal, grpSource, 20, 'none', 'syst123', 5);

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

