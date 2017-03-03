function [structSampleAreas, meanAmplitude, stdAmplitude] = Slotting(timeSeries, windowWidth, windowOverlap, kernelSpread)

%% Diff
timeSet = timeSeries(:,1);
ampSet = timeSeries(:,2);

%% Standardize the amplitude
meanAmplitude = mean(ampSet);
stdAmplitude = std(ampSet);

ampSet = (ampSet - meanAmplitude)/stdAmplitude;

%% Parse Model Into Continous Segements
x = ampSet;
t = timeSet - min(timeSet);

% 75% Overlap
kernelCenterArray = 0:windowWidth/windowOverlap:(max(t) + windowWidth); 
structSampleAreas = [];
gaussianMeanSet = [];
count = 1;

i = 1;
while(i < length(kernelCenterArray))
    % Gaussian Window
    idx = t > kernelCenterArray(i) - windowWidth & t < kernelCenterArray(i) + windowWidth;
    inKernelX = x(idx);
    inKernelT = t(idx);

    if(isempty(inKernelX));
        if(isempty(gaussianMeanSet))
            currentPt = find(t > kernelCenterArray(i) + windowWidth, 1);

            % new start (hop points)
            i = find(t(currentPt) > kernelCenterArray, 1, 'last');

        else
            structSampleAreas(count).set = gaussianMeanSet;
            count = count + 1;
            gaussianMeanSet = [];
        end
    else  

        weights = exp(-(((inKernelT - kernelCenterArray(i)).^2)*kernelSpread));
        gaussianMean = sum(weights.*inKernelX)/sum(weights);
        gaussianMeanSet = vertcat(gaussianMeanSet, gaussianMean);
    end

    i = i + 1;
end
    
structSampleAreas(count).set = gaussianMeanSet;
end