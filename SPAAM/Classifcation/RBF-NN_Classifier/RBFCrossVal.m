function [errorEst, percentReject, classEstimate] = RBFCrossVal(fltPatternArray_CrossVal, grpSource_CrossVal, structLRC, uniqueGrps, nodes, spread, fltRejectionValue)
%% Cross Val
    testing = 0;
    intLengthCrossVal = length(fltPatternArray_CrossVal(:,1));
    fltRBCrossData = zeros(intLengthCrossVal, length(nodes));

    for i = 1:1:intLengthCrossVal
        fltPositionVector = fltPatternArray_CrossVal(i,:);
        [fDistance] = GaussianKernel(fltPositionVector, nodes, spread);
        fltRBCrossData(i, :) = fDistance;
    end
    
    [errorEst, percentReject, classEstimate] = LRCCrossVal(fltRBCrossData, grpSource_CrossVal, uniqueGrps, structLRC, fltRejectionValue, 0);
    
    
    if(testing == 1)
        %% Testing
        maxSizeA1  = max(fltPatternArray_CrossVal);
        minSizeA1  = min(fltPatternArray_CrossVal);
        [X,Y] = meshgrid(linspace(minSizeA1(1),maxSizeA1(1)),linspace(minSizeA1(2),maxSizeA1(2)));
        X = X(:); Y = Y(:);
        fltTesting = cat(2, X, Y);

        for i = 1:1:length(fltTesting)
            fltPositionVector = fltTesting(i,:);
            [fDistance] = GaussianKernel(fltPositionVector, nodes, spread);
            fltRBTestData(i, :) = fDistance;
        end
       
        [fltNewTesting] = InitializeTrainingData(fltRBTestData);
        
        [classEstimate] = GenerateClassEstimates(fltNewTesting, uniqueGrps, structLRC);

        figure
        gscatter(fltTesting(:,1), fltTesting(:,2), classEstimate);
        hold on
        gscatter(fltPatternArray_CrossVal(:,1), fltPatternArray_CrossVal(:,2), grpSource_CrossVal)
        hold off
    end
    
    
end