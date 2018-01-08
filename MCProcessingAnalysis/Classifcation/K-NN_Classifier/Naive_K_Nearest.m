% =======================================================
% Copyright (c) 2005, Kyle Johnston
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% The views and conclusions contained in the software and documentation are those
% of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of the FreeBSD Project.
% =======================================================
function [fltResponse, classEstimate] = Naive_K_Nearest(...
    fltPatternArray_Training, grpSource_Training, fltPatternArray_CrossVal, ...
    intKN, p, fltSpread,  strMissValue)
% function [fltResponse, classEstimate] = Naive_K_Nearest(...
%     fltPatternArray_Training, grpSource_Training, fltPatternArray_CrossVal, ...
%     intKN, p, fltSpread,  strMissValue)
%
% Author: Kyle Johnston     14 Jan 2010
%
% Usage: This function generates class estimates based on the K-NN
%   algorithm. Ties are represented by the "strMissValue" string in the
%   output labels. This particular implimentation is a weighted K-NN
%   algorithm based on the distance from the input pattern to other
%   training pattterns.
% 
% Input: 
%       fltPatternArray_Training: the training dataset
%       grpSource_Training: the labels for the training dataset
%       fltPatternArray_CrossVal: the test dataset
%       intKN: value of k
%       p: distance type (2 is quadratic distance)
%       fltSpread: for the distance measurement
%       strMissValue: strnig label for rejected points (ties)
% Output:
%       fltResponse: the set of posterior probabilities
%       classEstimate: estimated labels for fltPatternArray_CrossVal


% Initialize Variables and Arrays
intNumberOfDimensionsTraining = length(fltPatternArray_Training(1,:));
intNumberOfTrainingPatterns = length(fltPatternArray_Training(:,1));

intNumberOfDimensionsTest = length(fltPatternArray_CrossVal(1,:));
intNumberOfTestPatterns = length(fltPatternArray_CrossVal(:,1));

classEstimate = cell(1,intNumberOfTestPatterns);

grpUniqueTypes = unique(grpSource_Training);
intNumberOfUnique = length(grpUniqueTypes);

fltResponse = zeros(intNumberOfTestPatterns, intNumberOfUnique);

% Check to make sure two datasets operate in the same vector space
if(intNumberOfDimensionsTest == intNumberOfDimensionsTraining)

    % Initialize distance arrays
    fltDistance = zeros(1,intNumberOfTrainingPatterns);

    % Loop through test patterns provided
    for i = 1:1:intNumberOfTestPatterns

        % Compute distance between training and test
        for j = 1:1:intNumberOfTrainingPatterns

            fltDistanceSet = zeros(1,intNumberOfDimensionsTest);
            for k = 1:1:intNumberOfDimensionsTest
                fltDistanceSet(k) = (fltPatternArray_CrossVal(i,k) - fltPatternArray_Training(j,k))^p;
            end

            fltDistance(j) = (sum(fltDistanceSet))^(1/p);
        end

        %Find the K closest points
        for j = 1:1:intKN
            indexCurrent = find(min(fltDistance) == fltDistance);

            %Weighted scheme of associated based on distance
            for k = 1:1:intNumberOfUnique
                if(strcmp(grpUniqueTypes(k),grpSource_Training(indexCurrent(1))))
                    fltResponse(i,k) = fltResponse(i,k) + 1;%% - normcdf(min(fltDistance), 0 , fltSpread);
                    break;
                end
            end

            fltDistance = fltDistance(min(fltDistance) ~= fltDistance);
        end

        % Compute posterior probabilities
        fltResponse(i,:) = fltResponse(i,:)./sum(fltResponse(i,:));
    end

end

% Label the datasets basedon posterior probabilities

for i = 1:1:length(fltResponse(:,1))

    fltMaxProb = max(fltResponse(i,:));
    intSingle = sum(fltResponse(i,:) == fltMaxProb);
    index = find(fltMaxProb == fltResponse(i,:));

    if(intSingle ~= 1)
        classEstimate{i} = strMissValue;
    else
        classEstimate{i} = grpUniqueTypes{index};
    end
end

end