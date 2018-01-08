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
function [errorArray, resArray] = SpaamAnalysis_Big(structMCReduct, grpSource,  structPatternArray_TrainingCV)

%% Sort By Classes
uniqueClasses = unique(grpSource);
structUniqueClassesIdx = [];

for indexk = 1:1:length(uniqueClasses)
    idxSet = [];
    
    for indexj = 1:1:length(grpSource)
        if(strcmp(grpSource{indexj}, uniqueClasses{indexk}))
            idxSet = vertcat(idxSet, indexj);
        end
    end
    
    structUniqueClassesIdx(indexk).idxSet = idxSet;
end



%%
resArray = 0.02:0.01:0.13;
% resArray = 5:1:35;
% resArray = 0.09;
errorArray = zeros(length(resArray), 2);

%% Cycle over resolutions, determine optimal resolution

for indexi = 1:1:length(resArray) 

    %% construct state space on the training set, for given resolution.
    states = -2:resArray(indexi):2;
    
%     logSpace = logspace(log10(0.01),log10(2),resArray(indexi));
%     states = horzcat(-fliplr(logSpace), logSpace);
%     
    intCluster = 15;
    
    spaceTotal = zeros(1,length(states)*length(states) + 2);
    grpSource_Partial = cell(1,1);
    counter = 1;
    %% simplify for given set using k-means to outline divisions
    for indexh = 1:1:length(uniqueClasses)
        idxSet = structUniqueClassesIdx(indexh).idxSet;

%         spacePartial = zeros(length(idxSet) ,length(states)*length(states)  + 2 + 10);
        spacePartial = zeros(length(idxSet) ,length(states)*length(states)  + 2);

        structMCReduct_Temp = structMCReduct(idxSet);

        for i = 1:1:length(structMCReduct_Temp)
            [markovChain, ~, meanAmplitude, stdAmplitude] = ...
                ConstructMCAmp_Kernel(structMCReduct_Temp(i).timeSeries, states);

%             [fltSet] = ComputeTextureMetrics(markovChain);
            
            structMCReduct_Temp(i).MC = markovChain;
            structMCReduct_Temp(i).unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
            spacePartial(i,:) = [structMCReduct_Temp(i).unpackMC', meanAmplitude, stdAmplitude];
%             spacePartial(i,:) = [structMCReduct_Temp(i).unpackMC', meanAmplitude, stdAmplitude, fltSet];
        end

         structPsi = K_Means(spacePartial', intCluster);
%          [structPsi] = Sampling(spacePartial, intCluster);

        for i = 1:1:length(structPsi)
            spaceTotal(counter,:) = structPsi(i).mean;
            grpSource_Partial{counter} = uniqueClasses{indexh};
            counter = counter + 1;
        end

    end


    %% Extended CVA (performs the 5-Fold Cross Validation as part of the analysis
    ECVAm = ecva('model', spaceTotal, grpSource_Partial, 10, 'none', 'syst123', 5, [], 0);
    fltPatternArray = ECVAm.CanonicalVariates;

    [structPatternArray_TrainingCV_Reduced, indexTesting] = ...
        Generate_5FoldCrossVal(grpSource_Partial);


    %% CV and Error Estimates 
    errorProbMean = 0;
    for j = 1:1:length(structPatternArray_TrainingCV_Reduced)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArray, grpSource_Partial, structPatternArray_TrainingCV_Reduced, j);

        [errorProb, classEstimate, structParaEst] = LDA_QDA_Classifier(...
            fltPatternArray_Training, grpSource_Training', fltPatternArray_CrossVal, grpSource_CrossVal', 4); 

        errorProbMean = errorProbMean + errorProb/5;
    end

    errorArray(indexi,1) = resArray(indexi);
    errorArray(indexi,2) = errorProbMean;
end



end

function [structPsi] = Sampling(spacePartial, intCluster)

    intLength = length(spacePartial(:,1));
    
    structPsi = [];
    
    intSelectSample = unique(randi(intLength,[1,intCluster])); 
    
    for i = 1:1:length(intSelectSample)
        structPsi(i).mean = spacePartial(intSelectSample(i),:);
    end
    

end



