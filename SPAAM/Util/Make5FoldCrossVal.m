function [structPatternArray_TrainingCV] = Make5FoldCrossVal(grpSource_Training)

p = randperm(length(grpSource_Training));
sampleSize = round(length(p)/5) - 1;

structPatternArray_TrainingCV = [];
count = 1;
for j = 1:1:5
    indexSet = p(count:count + sampleSize - 1);
    count = count + sampleSize;
    structPatternArray_TrainingCV(j).indexSet = indexSet;
end


end