function [intLabelArray] = Use_LDA_QDA_Classifer(fltPatternArray, structParaEst)
%% Use Classifier
intNumberTesting = length(fltPatternArray(:,1));
intLabelArray = cell(intNumberTesting,1);

for i = 1:1:intNumberTesting
    [structProb] = DetermineDistanceEstimates(fltPatternArray(i,:), structParaEst);
    
    maxProb = structProb(1).prob;
    
    for j = 1:1:length(structProb)
        if(maxProb <= structProb(j).prob)
            maxProb = structProb(j).prob;
            currentType = structProb(j).type;
        end
    end
    
    
    intLabelArray{i} = cell2mat(currentType);
end


end
