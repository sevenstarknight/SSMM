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
function [fltPatternArray_Training, grpSource_Training, ...
    fltPatternArray_CrossVal, grpSource_CrossVal] = ...
    PullTrainingAndCrossMatrixFromStruct(fltPatternArray, grpSource, ...
    structPatternArray_TrainingCV, index)

fltPatternArray_Training = [];
grpSource_Training = {};
fltPatternArray_CrossVal = [];
grpSource_CrossVal = {};

for i = 1:1:length(structPatternArray_TrainingCV) 
    
   if(i == index)
       cvSet = structPatternArray_TrainingCV(i).indexSet;
       fltPatternArray_CrossVal = fltPatternArray(:,:, cvSet);
       grpSource_CrossVal = grpSource(cvSet);
   else
       trainSet = structPatternArray_TrainingCV(i).indexSet;
       newSet = fltPatternArray(:,:, trainSet);
       grpNewSet = grpSource(trainSet);
       fltPatternArray_Training = cat(3,fltPatternArray_Training, newSet);
       grpSource_Training = horzcat(grpSource_Training, grpNewSet);
   end
    
    
end

end