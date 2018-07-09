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
clc
clear

location = pwd;
stringDirectory = strcat(location, '\SSMM\matlab\data\dataset2');

stringDirectory = '/Users/kjohnston/Google Drive/VarStarData/UCR';

%% ====================================================
% % Make feature space and load data
[structMCReduct, grpSource] = ReadInTestTimeSeries(stringDirectory);

%% ====================================================
% Make Cross Validation
[structPatternArray_TrainingCV, indexTesting] =  Generate_5FoldCrossVal(grpSource);

%% ====================================================
% Cycle over resolutions, determine optimal resolution
% [errorArray, resArray] = SpaamAnalysis_Metric(structMCReduct, grpSource, structPatternArray_TrainingCV);

%% ======================================================================
% Generate Feature Space JUST on Variablity
states = -2:0.18:2;
spaceTotal = zeros(length(structMCReduct) ,length(states)*length(states) + 2);
for i = 1:1:length(structMCReduct)
    [markovChain, stateTransitionMatrix, meanAmplitude, stdAmplitude] = ConstructMCAmp_Kernel(...
        structMCReduct(i).timeSeries,states);
        
    structMCReduct(i).MC = markovChain;
    structMCReduct(i).unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
    spaceTotal(i,:) = [structMCReduct(i).unpackMC',meanAmplitude, stdAmplitude];
end

idxTrainCross = horzcat(structPatternArray_TrainingCV(1).indexSet, structPatternArray_TrainingCV(2).indexSet,...
    structPatternArray_TrainingCV(3).indexSet, structPatternArray_TrainingCV(4).indexSet, ...
    structPatternArray_TrainingCV(5).indexSet);

% changing from index on everything to index on the sub-sample
equiSet = length(structPatternArray_TrainingCV(1).indexSet);
structPatternArray_TrainingCross(1).indexSet = 1:equiSet*1;
structPatternArray_TrainingCross(2).indexSet = equiSet + 1:equiSet*2;
structPatternArray_TrainingCross(3).indexSet = equiSet*2 + 1:equiSet*3;
structPatternArray_TrainingCross(4).indexSet = equiSet*3 + 1:equiSet*4;
structPatternArray_TrainingCross(5).indexSet = equiSet*4 + 1:equiSet*5;

fltPatternArray_TrainCross = spaceTotal(idxTrainCross, :); 
grpSource_TrainCross = grpSource(idxTrainCross);

%% =====================================
% Extended CVA Applied to Training and Crossval (Total, with 5-Fold Validation)      
ECVAm = ecva('model', fltPatternArray_TrainCross, grpSource_TrainCross, 20, 'none', 'syst123', 5);
fltPatternArrayMC = ECVAm.CanonicalVariates;
   
mx = ECVAm.Detail.mx; % To be used for predictions
canonicalWeights = ECVAm.CanonicalWeights; % To be used for predictions


fltPatternArrayMerge = fltPatternArrayMC;

% standardize the merged dataset
vectorMean = mean(fltPatternArrayMerge,1);
vectorStd = std(fltPatternArrayMerge, 1);
for i = 1:1:length(fltPatternArrayMerge)
    
    fltPatternArrayMerge(i,:) = (fltPatternArrayMerge(i,:) - vectorMean)./vectorStd;
    
end


%% ======================================================================
% Graphics (EDA)

categories = {'ECV1', 'ECV2'};

parallelcoords(fltPatternArrayMerge, 'Group', grpSource_TrainCross, ...
    'Labels', categories)
grid on


gplotmatrix(fltPatternArrayMerge,fltPatternArrayMerge,grpSource_TrainCross',...
    'brgkym','.',[],'on','',categories',...
    categories');

%% ==================================================================
% Grab Testing Data
fltPatternArray_Testing = spaceTotal(indexTesting, :);
grpSource_Testing = grpSource(indexTesting);

% transform testing data
for i = 1:1:length(grpSource_Testing)
    ECVATesting(i, :) = (fltPatternArray_Testing(i, :) - mx);
end
fltPatternArray_Testing_Canonical = ECVATesting*canonicalWeights;


fltPatternArray_Testing = fltPatternArray_Testing_Canonical;

for i = 1:1:length(fltPatternArray_Testing)
    fltPatternArray_Testing(i,:) = (fltPatternArray_Testing(i,:) - vectorMean)./vectorStd;
end
grpSource_Testing = grpSource(indexTesting);

%% ===================================================================
[confusionMatrix, fltSpread] = SpaamKNN(fltPatternArrayMerge, grpSource_TrainCross,...
    fltPatternArray_Testing, grpSource_Testing, structPatternArray_TrainingCross);

%% ===================================================================
%PWC
[confusionMatrix, fltSpread] = SpaamPWC(fltPatternArrayMerge, grpSource_TrainCross,...
    fltPatternArray_Testing, grpSource_Testing, structPatternArray_TrainingCross);

%% ==================================================
% RBF-NN
[confusionMatrix, structLRC, nodes, spread] = SpaamRBFNN(fltPatternArrayMerge, grpSource_TrainCross,...
    fltPatternArray_Testing, grpSource_Testing, structPatternArray_TrainingCross);


%% ==================================================
% Bagged decision trees
[confusionMatrix, B] = SpaamBaggedTree(fltPatternArrayMerge, grpSource_TrainCross,...
    fltPatternArray_Testing, grpSource_Testing, structPatternArray_TrainingCross);


%% ==================================================
%  KNFST
[model, nodes, errorEst] = ...
    SpaamKNFST(fltPatternArrayMerge, grpSource_TrainCross, ...
    fltPatternArray_Testing, structPatternArray_TrainingCross);


%% ==================================================
% Experiment
[fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArrayMC,...
            grpSource_TrainCross, structPatternArray_TrainingCross, 1);

[model, nodes, errorEst, scores] = ...
    KNFSTDetector(fltPatternArray_Training, fltPatternArray_Testing_Canonical, 1.3);

[model, nodes, errorEst, scores] = ...
    KNFSTDetector(fltPatternArray_Testing_Canonical, fltPatternArray_Training, 1.3);

[model, nodes, errorEst, scores] = ...
    KNFSTDetector(fltPatternArray_Training, fltPatternArray_CrossVal, 1.3);


%% ===================================================================
%OC-PWC
[w, errorEst] = ...
    SpaamOCSVM(fltPatternArrayMerge, grpSource_TrainCross, ...
    fltPatternArray_Testing, structPatternArray_TrainingCross);

