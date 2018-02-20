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