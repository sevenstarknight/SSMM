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
function[w, errorEst] = SpaamOCSVM(fltPatternArrayMerge, ...
    grpSource_TrainCross, fltPatternArray_Testing, structPatternArray_TrainingCross)


spreadArray = 10.^(-1:0.05:1);
errorArray = zeros(1,length(spreadArray));
uniqueSet = unique(grpSource_TrainCross);

for i = 1:1:length(spreadArray)
    errorProbMean = 0;
    
    for j = 1:1:length(structPatternArray_TrainingCross)

        [fltPatternArray_Training, grpSource_Training, ...
            fltPatternArray_CrossVal, grpSource_CrossVal] = ...
            PullTrainingAndCrossFromStruct(fltPatternArrayMerge,...
            grpSource_TrainCross, structPatternArray_TrainingCross, j);

        training_oc = oc_set(fltPatternArray_Training);
        crossval_oc = oc_set(fltPatternArray_CrossVal);
      
        w = parzen_dd(training_oc, 0.000, spreadArray(i));

        a = crossval_oc*w*labeld;
          
        counter = 0;
        for k = 1:1:length(a)
            if(strcmp(a(k,:),'outlier'))
                counter = counter + 1;
            end
        end
        errorEst = counter/length(a);
        
        errorProbMean = errorProbMean + errorEst/5;
    end
    
    errorArray(i) = errorProbMean;
end

figure()
semilogx(spreadArray, errorArray, '.-r')
xlabel('OC-PWC Kernel Width')
ylabel('Misclassification Rate 5-Fold Cross-validation')
grid on

testing_oc = oc_set(fltPatternArray_Testing);
w = parzen_dd(training_oc, 0.001, 2.5);
a = testing_oc*w*labeld;

counter = 0;
for k = 1:1:length(a)
    if(strcmp(a(k,:),'outlier'))
        counter = counter + 1;
    end
end
errorEst = counter/length(a);



end