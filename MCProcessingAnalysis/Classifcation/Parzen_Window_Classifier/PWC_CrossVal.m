function [classEstimate] = PWC_CrossVal(fltPatternArray_CrossVal, grpSource_CrossVal, spread)

    %% Testing
    maxSizeA1  = max(fltPatternArray_CrossVal);
    minSizeA1  = min(fltPatternArray_CrossVal);
    [X,Y] = meshgrid(linspace(minSizeA1(1),maxSizeA1(1), 100),linspace(minSizeA1(2),maxSizeA1(2), 100));
    X = X(:); Y = Y(:);
    fltTesting = cat(2, X, Y);

    [~, classEstimate] = Parzen_Window_Classifier(...
        fltPatternArray_CrossVal, grpSource_CrossVal, fltTesting, ...
        3, spread, 'Missed');
    
    
    figure
    gscatter(fltTesting(:,1), fltTesting(:,2), classEstimate', 'br','.', 2);
    hold on
    gscatter(fltPatternArray_CrossVal(:,1), fltPatternArray_CrossVal(:,2), grpSource_CrossVal, 'br','xo', 6)
    hold off

end