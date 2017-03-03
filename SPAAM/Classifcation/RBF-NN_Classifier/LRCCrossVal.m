function [errorEst, percentReject, classEstimate] = LRCCrossVal(fltPatternArray_CrossVal, grpSource_CrossVal, uniqueGrps, structLRC, fltRejectionValue, testing)

%% Cross Validation
[fltNewCrossValData] = InitializeTrainingData(fltPatternArray_CrossVal);

[classEstimate] = GenerateClassEstimates(fltNewCrossValData, uniqueGrps, structLRC, fltRejectionValue);

correct = 0;
intReject = 0;

for i = 1:1:length(classEstimate)
    if(strcmp(classEstimate{i},grpSource_CrossVal{i}))
        correct = correct + 1;
    elseif(strcmp(classEstimate{i},'Reject'))
        intReject = intReject + 1;
    else
        
    end
end

errorEst = 1 - (correct)/length(grpSource_CrossVal);

percentReject = intReject/length(grpSource_CrossVal);

if(testing == 1)
    %% Testing
    maxSizeA1  = max(fltPatternArray_CrossVal);
    minSizeA1  = min(fltPatternArray_CrossVal);
    [X,Y] = meshgrid(linspace(minSizeA1(1),maxSizeA1(1), 100),linspace(minSizeA1(2),maxSizeA1(2), 100));
    X = X(:); Y = Y(:);
    fltTesting = cat(2, X, Y);

    [fltNewTesting] = InitializeTrainingData(fltTesting);

    [classEstimate] = GenerateClassEstimates(fltNewTesting, uniqueGrps, structLRC);

    figure
    gscatter(fltTesting(:,1), fltTesting(:,2), classEstimate, 'brgk','.', 2);
    hold on
    gscatter(fltPatternArray_CrossVal(:,1), fltPatternArray_CrossVal(:,2), grpSource_CrossVal, 'gkrb','o', 6)
    hold off
end


end