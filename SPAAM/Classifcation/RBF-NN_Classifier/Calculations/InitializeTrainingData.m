function [fltNewTrainingData] = InitializeTrainingData(fltPatternArray_Training)
    
intLengthTraining = length(fltPatternArray_Training(:,1));
interceptArray = zeros(1,intLengthTraining) + 1;

fltNewTrainingData = cat(2, interceptArray', fltPatternArray_Training);

end