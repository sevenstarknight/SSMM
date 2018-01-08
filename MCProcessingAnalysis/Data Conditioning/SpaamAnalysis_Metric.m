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
function [errorArray, resArray] = SpaamAnalysis_Metric(structMCReduct, grpSource, structPatternArray_TrainingCV)

resArray = 0.02:0.005:0.2;
resArray = fliplr(resArray);
errorArray = zeros(length(resArray),3);
uniqueSources = unique(grpSource);
%% Cycle over resolutions, determine optimal resolution
for indexi = 1:1:length(resArray) 
    
    %======================================================
    %% construct state space on the training set, for given resolution.
    states = -2:resArray(indexi):2;    
    stMatrixGroup = zeros(length(states) ,length(states), length(uniqueSources) );
    mcMatrix = zeros(length(states) ,length(states), length(structMCReduct) ); 
    for indxi = 1:1:length(structMCReduct)
        [markovChain, stateTransitionMatrix, meanAmplitude, stdAmplitude] = ...
            ConstructMCAmp_Kernel(structMCReduct(indxi).timeSeries, states);
        
        structMCReduct(indxi).MC = markovChain;
        structMCReduct(indxi).timeDomain = [meanAmplitude, stdAmplitude];
        
        idx = strcmp(grpSource(indxi), uniqueSources);
        stMatrixGroup(:,:,idx) = stMatrixGroup(:,:,idx) + stateTransitionMatrix;
        
        mcMatrix(:,:,indxi) = markovChain;
    end
        
    
    %% turn counts into probablities
    markovChainGroup = zeros(length(states), length(states), length(uniqueSources));
    for indxi = 1:1:length(uniqueSources)
        totalInState = sum(stMatrixGroup(:,:,indxi), 2);
        markovChain = zeros(length(states), length(states));
        for indxj = 1:1:length(totalInState)
            if(totalInState(indxj) ~= 0.0)
                markovChain(indxj,:) = (stMatrixGroup(indxj,:,indxi)./totalInState(indxj));
            end
        end
        markovChainGroup(:,:,indxi) = markovChain;
    end
    
    fltPatternArray = zeros(length(structMCReduct), length(markovChainGroup(1,1,:)));
    for indxi = 1:1:length(structMCReduct)
        for indxj = 1:1:length(markovChainGroup(1,1,:))
            delta = structMCReduct(indxi).MC - markovChainGroup(:,:,indxj);
            fltPatternArray(indxi, indxj) = norm(delta,'fro');
        end
    end
    

    %% CV and Error Estimates 
    errorProbMean = 0;
    errorProbMeanMetric = 0;
    
    meanTimeLDA = 0;
    meanPullTime = 0;
    meanPushTime = 0;
    meanKNNTime = 0;
    
    for j = 1:1:length(structPatternArray_TrainingCV)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArray, grpSource, structPatternArray_TrainingCV, j);

        tic
        [errorProb, ~, ~] = LDA_QDA_Classifier(...
            fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, grpSource_CrossVal', 4); 
        ldaTime = toc;
        
        errorProbMean = errorProbMean + errorProb/5;
        
        %% ===================================
        
        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossMatrixFromStruct(mcMatrix, ...
            grpSource, structPatternArray_TrainingCV, j);
        
        [M, pullTime, pushTime] = pushPullMethod(fltPatternArray_Training,...
            grpSource_Training, 1.0, 0.5);
        
        options = struct('kValue', 1, 'M', M);
        
        tic
        [~, accuracy] = DF_KNN (...
            fltPatternArray_Training, grpSource_Training,...
            fltPatternArray_CrossVal, grpSource_CrossVal, options);
        
        dfknnTime = toc;
        errorProbMeanMetric = errorProbMeanMetric + (1-accuracy)/5;
        
        meanTimeLDA = meanTimeLDA + ldaTime/5;
        meanPullTime = meanPullTime + pullTime/5;
        meanPushTime = meanPushTime + pushTime/5;
        meanKNNTime = meanKNNTime + dfknnTime/5;
        
    end

    errorArray(indexi,1) = resArray(indexi);
    errorArray(indexi,2) = errorProbMean;
    errorArray(indexi,3) = errorProbMeanMetric;
    
    disp(strcat(num2str(meanTimeLDA) , ' ; ' , num2str(meanPullTime) ...
            , '  ; ' , num2str(meanPushTime) , '  ; ' ,+ num2str(meanKNNTime)));
        
    disp(strcat(num2str(resArray(indexi)) ,...
        '  ; ' , num2str(errorProbMean) , ...
        '  ; ' , num2str(errorProbMeanMetric)));
end


end

