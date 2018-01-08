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

stringDirectory = 'C:\Users\kyle.johnston\Desktop\Dissertation\starvars\matlab\data\trainingsets';

%% ====================================================
% Make feature space and load data
[structMCReduct, grpSource] = ReadInTimeSeries(stringDirectory);

%% ====================================================
% Add in flat lines
[structMCReduct, grpSource] = ConstructNonVariableStars(1000, structMCReduct, grpSource);

%% ====================================================
% Make Cross Validation
[structPatternArray_TrainingCV, indexTesting] =  Generate_5FoldCrossVal(grpSource);


windowWidthArray = 0.005:0.005:0.025;
resArray = 0.02:0.01:0.14;
% resArray = 0.03;
% resArray = 5:2:30;
errorArray = zeros(length(resArray),length(windowWidthArray));

%% Cycle over resolutions, determine optimal resolution

for indexj = 1:1:length(windowWidthArray)

for indexi = 1:1:length(resArray) 
    %======================================================
    %% construct state space on the training set, for given resolution.
    
    states = -2:resArray(indexi):2;
    
%     logSpace = logspace(log10(0.1),log10(3),resArray(indexi));
%     states = horzcat(-fliplr(logSpace), logSpace);
    
%     spaceTotal = zeros(length(structMCReduct) ,length(states)*length(states) + 2 + 10);
    spaceTotal = zeros(length(structMCReduct) ,length(states)*length(states) + 2 );
    for i = 1:1:length(structMCReduct)
        [markovChain, meanAmplitude, stdAmplitude] = ConstructMCAmp_Kernel(structMCReduct(i).timeSeries, states, windowWidthArray(indexj));

%         [fltSet] = ComputeTextureMetrics(markovChain);
        
        structMCReduct(i).MC = markovChain;
        structMCReduct(i).unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
%         spaceTotal(i,:) = [structMCReduct(i).unpackMC', meanAmplitude, stdAmplitude, fltSet];
        spaceTotal(i,:) = [structMCReduct(i).unpackMC',meanAmplitude, stdAmplitude];
    end

    y = randsample(length(structMCReduct),floor(length(structMCReduct)*0.80));
    
    %% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
    ECVAm = ecva('model', spaceTotal(y,:), grpSource(y), 10, 'none', 'syst123', 5, [], 0);

    mx = ECVAm.Detail.mx; % To be used for predictions
    canonicalWeights = ECVAm.CanonicalWeights; % To be used for predictions

    % transform testing data
    ECVATesting = zeros(length(grpSource), length(mx));
    for i = 1:1:length(grpSource)
        ECVATesting(i, :) = (spaceTotal(i, :) - mx);
    end
    fltPatternArray_Testing_Canonical = ECVATesting*canonicalWeights;
    
    %% CV and Error Estimates 
    errorProbMean = 0;
    for j = 1:1:length(structPatternArray_TrainingCV)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArray_Testing_Canonical, grpSource, structPatternArray_TrainingCV, j);

        [errorProb, classEstimate, structParaEst] = LDA_QDA_Classifier(...
            fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, grpSource_CrossVal', 4); 

        errorProbMean = errorProbMean + errorProb/5;
    end

    errorArray(indexi,indexj) = errorProbMean;
end

end