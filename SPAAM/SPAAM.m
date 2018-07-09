% ==================================================== 
%  Copyright (C) 2016 Kyle Johnston
%  
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ====================================================
function [errorProb_NN] = SPAAM(stringDirectory, strLabel)
%% ====================================================
% Make feature space and load data
% 
% [structTraining, grpSource_Training, structTesting, grpSource_Testing] = ReadInLINEARTimeSeries(stringDirectory);

[structTraining, grpSource, structTesting, grpSource_Testing] = ReadInUCRTimeSeries(stringDirectory, strLabel);

%% ====================================================
uniqueClasses = unique(grpSource);
numArray = zeros(1,length(uniqueClasses));

for i = 1:1:length(uniqueClasses)
    label = uniqueClasses{i};
    numArray(i) = length(find(ismember(grpSource,label)));
end

if(min(numArray) < 6)
    errorProb_NN = 1;
    return;
end
%% ====================================================
%% Randomly distribute into five sets
[structPatternArray_TrainingCross] = Make5FoldCrossVal(grpSource);

%% ====================================================
% Optimize
resArray = 0.03:0.01:0.2;
% kernelArray = 10.^(2:0.2:4);
kernelArray = 10;
%% =====================================================
% Slot
% windowWidth = 0.05; windowOverlap = 4;
% tic
% for idx = 1:1:length(structTraining)
%     timeSeries = structTraining(idx).timeSeries;
%     structSlots(length(kernelArray)) = struct('structSampleAreas', []);
%     for jdx = 1:1:length(kernelArray)
%         [structSampleAreas, meanAmplitude, stdAmplitude] = Slotting(timeSeries, windowWidth, windowOverlap, kernelArray(jdx));
%         structSlots(jdx).structSampleAreas = structSampleAreas;
%     end
%     structTraining(idx).structSlots = structSlots;
%     structTraining(idx).meanAmplitude = meanAmplitude;
%     structTraining(idx).stdAmplitude = stdAmplitude;
% end

for idx = 1:1:length(structTraining)
    
    ampSet = structTraining(idx).timeSeries(:,2);

    %% Standardize the amplitude
    meanAmplitude = mean(ampSet);
    stdAmplitude = std(ampSet);

    ampSet = (ampSet - meanAmplitude)/stdAmplitude;
    
    
    structSampleAreas(1).set = ampSet;
    structSlots(1).structSampleAreas = structSampleAreas;
    structTraining(idx).structSlots = structSlots;
    structTraining(idx).meanAmplitude = meanAmplitude;
    structTraining(idx).stdAmplitude = stdAmplitude;
end

%% ====================================================
% construct state space, for given resolution.
errorArray = zeros(length(kernelArray),length(resArray));

for idx = 1:1:length(resArray)
    states = -2:resArray(idx):2;
    for jdx = 1:1:length(kernelArray)
        
        spaceTotal = zeros(length(structTraining) ,length(states)*length(states) + 2);
        
        for kdx = 1:1:length(structTraining)
            structSlots = structTraining(kdx).structSlots;
            [markovChain] = MarkovChain(structSlots(jdx).structSampleAreas, states);
            [unpackMC] = reshape(markovChain, [length(states)*length(states), 1]);
            
            spaceTotal(kdx,:) = [unpackMC', structTraining(kdx).meanAmplitude, structTraining(kdx).stdAmplitude];
        end
        
         %% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
        ECVAm = ecva('model', spaceTotal, grpSource, 20, 'none', 'syst123', 5, [], 0);
        fltPatternArray = ECVAm.CanonicalVariates;

        %% CV and Error Estimates 
        errorProbMean = 0;
        for kdx = 1:1:length(structPatternArray_TrainingCross)

            [fltPatternArray_Training, grpSource_Training, ...
                fltPatternArray_CrossVal, grpSource_CrossVal] = ...
                PullTrainingAndCrossFromStruct(fltPatternArray, grpSource, structPatternArray_TrainingCross, kdx);

%             [errorProb, ~, ~] = LDA_QDA_Classifier(...
%                 fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, grpSource_CrossVal', 4); 
%             
            
            [~, classEstimate] = Naive_K_Nearest(...
                    fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, ...
                    1.0, 2, 1.0,  'Missed');

            count = 0;
            for k = 1:1:length(classEstimate)
                if(strcmp(classEstimate{k}, grpSource_CrossVal{k}))
                    count = count + 1;
                end
            end

            errorProb = (length(classEstimate) - count)/length(classEstimate);
            

            errorProbMean = errorProbMean + errorProb/5;
        end

        errorArray(jdx,idx) = errorProbMean;
        
    end
end

 %% Generate Feature Space, Find Optimal State Resolution
[~,iResolution] = min(min(errorArray));
% [~,iKernelSet] = min(errorArray);
% iKernel = iKernelSet(iResolution);

states = -2:resArray(iResolution):2;
% kernelWidth = kernelArray(iKernel);

spaceTotal = zeros(length(structTraining) ,length(states)*length(states) + 2);

% for idx = 1:1:length(structTraining)
%     timeSeries = structTraining(idx).timeSeries;
%     
%     [structSampleAreas, meanAmplitude, stdAmplitude] = ...
%         Slotting(timeSeries, windowWidth, windowOverlap, kernelWidth);
% 
%     structTraining(idx).structSampleAreas = structSampleAreas;
%     structTraining(idx).meanAmplitude = meanAmplitude;
%     structTraining(idx).stdAmplitude = stdAmplitude;
% end


for idx = 1:1:length(structTraining)
    
    ampSet = structTraining(idx).timeSeries(:,2);

    %% Standardize the amplitude
    meanAmplitude = mean(ampSet);
    stdAmplitude = std(ampSet);

    ampSet = (ampSet - meanAmplitude)/stdAmplitude;
    
    
    structSampleAreas(1).set = ampSet;
    structTraining(idx).structSampleAreas = structSampleAreas;
    structTraining(idx).meanAmplitude = meanAmplitude;
    structTraining(idx).stdAmplitude = stdAmplitude;
end


for kdx = 1:1:length(structTraining)
    structSampleAreas = structTraining(kdx).structSampleAreas;
    [markovChain] = MarkovChain(structSampleAreas, states);
    [unpackMC] = reshape(markovChain, [length(states)*length(states), 1]);
    spaceTotal(kdx,:) = [unpackMC', structTraining(kdx).meanAmplitude, structTraining(kdx).stdAmplitude];
end

%% =================================================
% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
ECVAm = ecva('model', spaceTotal, grpSource, 20, 'none', 'syst123', 6, [], 0);
fltPatternArray = ECVAm.CanonicalVariates;

mx = ECVAm.Detail.mx; % To be used for predictions
canonicalWeights = ECVAm.CanonicalWeights; % To be used for predictions

fltPatternArrayMerge = fltPatternArray;

% standardize the merged dataset
vectorMean = mean(fltPatternArrayMerge,1);
vectorStd = std(fltPatternArrayMerge, 1);
for indexk = 1:1:length(fltPatternArrayMerge)
    fltPatternArrayMerge(indexk,:) = (fltPatternArrayMerge(indexk,:) - vectorMean)./vectorStd;
end


%% =================================================
% Transform Test Data
% for idx = 1:1:length(structTesting)
%     timeSeries = structTesting(idx).timeSeries;
%     
%     [structSampleAreas, meanAmplitude, stdAmplitude] = Slotting(timeSeries, windowWidth, windowOverlap, kernelWidth);
%     
%     structTesting(idx).structSampleAreas = structSampleAreas;
%     structTesting(idx).meanAmplitude = meanAmplitude;
%     structTesting(idx).stdAmplitude = stdAmplitude;
% end

for idx = 1:1:length(structTesting)
    
    ampSet = structTesting(idx).timeSeries(:,2);

    %% Standardize the amplitude
    meanAmplitude = mean(ampSet);
    stdAmplitude = std(ampSet);

    ampSet = (ampSet - meanAmplitude)/stdAmplitude;
    
    
    structSampleAreas(1).set = ampSet;
    structTesting(idx).structSampleAreas = structSampleAreas;
    structTesting(idx).meanAmplitude = meanAmplitude;
    structTesting(idx).stdAmplitude = stdAmplitude;
end


fltPatternArray_Testing = zeros(length(structTesting), length(states)*length(states) + 2);

for kdx = 1:1:length(structTesting)
    structSampleAreas = structTesting(kdx).structSampleAreas;
    [markovChain] = MarkovChain(structSampleAreas, states);
    [unpackMC] = reshape(markovChain, [length(states)*length(states), 1]);
    fltPatternArray_Testing(kdx,:) = [unpackMC', structTesting(kdx).meanAmplitude, structTesting(kdx).stdAmplitude];
end

ECVATesting = zeros(size(fltPatternArray_Testing));
for indexk = 1:1:length(grpSource_Testing)
    ECVATesting(indexk, :) = (fltPatternArray_Testing(indexk, :) - mx);
end
fltPatternArray_Testing_Canonical = ECVATesting*canonicalWeights;
fltPatternArray_Testing = fltPatternArray_Testing_Canonical;

for indexk = 1:1:length(fltPatternArray_Testing)
    fltPatternArray_Testing(indexk,:) = (fltPatternArray_Testing(indexk,:) - vectorMean)./vectorStd;
end

toc
 %% ===================================================================
% 1-NN Testing
[~, classEstimate] = Naive_K_Nearest(...
        fltPatternArrayMerge, grpSource', fltPatternArray_Testing, ...
        1.0, 2, 1.0,  'Missed');

count = 0;
for k = 1:1:length(classEstimate)
    if(strcmp(classEstimate{k}, grpSource_Testing{k}))
        count = count + 1;
    end
end

errorProb_NN = (length(classEstimate) - count)/length(classEstimate);


end