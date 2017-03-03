function [structPatternArray_TrainingCV, indexTesting] =   ...
    Generate_5FoldCrossVal(grpSource)

intLength = length(grpSource);
uniqueTypes = unique(grpSource);
intNumTypes = length(uniqueTypes);

indexTrainingCrossVal = [];
indexTesting = [];

%% Pipe in each class into a separate 
for i = 1:1:intLength
    X = rand;
    
    for j = 1:1:intNumTypes
        if(strcmp(grpSource{i}, uniqueTypes{j}))
            %% Random sample into both training/cross-val and testing
            if(X > 0.5)
                if(isempty(indexTrainingCrossVal))
                    indexTrainingCrossVal(1) = i; 
                else
                    n = length(indexTrainingCrossVal);
                    indexTrainingCrossVal(n + 1) = i; 
                end
            else
                if(isempty(indexTesting))
                    indexTesting(1) = i;
                else
                    n = length(indexTesting);
                    indexTesting(n + 1) = i;
                end
            end
        end
    end
end

%% Randomly distribute into five sets
p = randperm(length(indexTrainingCrossVal));
sampleSize = round(length(p)/5) - 1;

structPatternArray_TrainingCV = [];
count = 1;
for j = 1:1:5
    indexSet = p(count:count + sampleSize - 1);
    count = count + sampleSize;
    structPatternArray_TrainingCV(j).indexSet = indexSet;
end


end