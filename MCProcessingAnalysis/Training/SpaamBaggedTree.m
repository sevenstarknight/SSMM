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
function [confusionMatrix, B, tpr] = SpaamBaggedTree(fltPatternArrayMerge, grpSource_TrainCross,...
    fltPatternArray_Testing, grpSource_Testing, structPatternArray_TrainingCross)


%% ==================================================
% Bagged decision trees

nTrees = 1:3:50;
errorArray = zeros(1,length(nTrees));

for i = 1:1:length(nTrees)
    errorProbMean = 0;
    
    for j = 1:1:length(structPatternArray_TrainingCross)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArrayMerge,...
            grpSource_TrainCross, structPatternArray_TrainingCross, j);

        B = TreeBagger(nTrees(i),fltPatternArray_Training, grpSource_Training, 'Method', 'classification',...
        'NVarToSample', length(fltPatternArray_CrossVal(1,:)) - 1);

        
        count = 0;
        for k = 1:1:length(grpSource_CrossVal)
            
            predChar1 = B.predict(fltPatternArray_CrossVal(k,:));
            
            if(strcmp(predChar1, grpSource_CrossVal{k}))
                count = count + 1;
            end
        end
        errorProb = (length(grpSource_CrossVal) - count)/length(grpSource_CrossVal);


        errorProbMean = errorProbMean + errorProb/5;
    end
    
    errorArray(i) = errorProbMean;
end

figure()
plot(nTrees, errorArray, '.-r')
xlabel('Number of Trees Generated')
ylabel('Misclassification Rate 5-Fold Cross-validation')
grid on


B = TreeBagger(43,fltPatternArrayMerge, grpSource_TrainCross', 'Method', 'classification',...
        'NVarToSample', length(fltPatternArrayMerge(1,:)) - 1);

classEstimate = {};
count = 0;
for k = 1:1:length(grpSource_Testing)
    
    classEstimate{k} = B.predict(fltPatternArray_Testing(k,:));

end

%% TPR 

tp = 0;

for i = 1:1:length(classEstimate)
    if(strcmp(classEstimate{i}, grpSource_Testing{i}))
        tp = tp + 1;
    end
end
     
tpr = tp/length(classEstimate);

%% Confusion Matrix Generation
uniqueClasses = unique(grpSource_TrainCross);

missedVector = zeros(1, length(uniqueClasses));
confusionMatrix = zeros(length(uniqueClasses), length(uniqueClasses));

for i = 1:1:length(classEstimate)
    if(strcmp(classEstimate{i}, 'Missed'))
        indexi = find(strcmp(grpSource_Testing{i}, uniqueClasses));
        missedVector(indexi) = missedVector(indexi) + 1;
    else
        indexj = find(strcmp(classEstimate{i}, uniqueClasses));
        indexi = find(strcmp(grpSource_Testing{i}, uniqueClasses));
        confusionMatrix(indexi, indexj) = confusionMatrix(indexi, indexj) + 1;
    end
end

rowSums = sum(confusionMatrix,2);
for i = 1:1:length(rowSums)
    confusionMatrix(i,:) = confusionMatrix(i,:)./rowSums(i);
end


