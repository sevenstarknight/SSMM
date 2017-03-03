function [fltPatternArray_Training, grpSource_Training, ...
    fltPatternArray_CrossVal, grpSource_CrossVal] = ...
    PullTrainingAndCrossFromStruct(fltPatternArray, grpSource, structPatternArray_TrainingCV, index)

fltPatternArray_Training = [];
grpSource_Training = {};
fltPatternArray_CrossVal = [];
grpSource_CrossVal = {};

for i = 1:1:length(structPatternArray_TrainingCV) 
    
   if(i == index)
       cvSet = structPatternArray_TrainingCV(i).indexSet;
       fltPatternArray_CrossVal = fltPatternArray(cvSet, :);
       grpSource_CrossVal = grpSource(cvSet);
   else
       trainSet = structPatternArray_TrainingCV(i).indexSet;
       newSet = fltPatternArray(trainSet,:);
       grpNewSet = grpSource(trainSet);
       fltPatternArray_Training = vertcat(fltPatternArray_Training, newSet);
       grpSource_Training = horzcat(grpSource_Training, grpNewSet);
   end
    
    
end



end