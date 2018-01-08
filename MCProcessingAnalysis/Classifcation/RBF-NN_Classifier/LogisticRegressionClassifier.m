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
function [structLRC, uniqueGrps] = LogisticRegressionClassifier(fltPatternArray_Training, grpSource_Training, method)
%% Initialize Parameters
uniqueGrps = unique(grpSource_Training);
intUniqueGroups = length(uniqueGrps);

[fltNewTrainingData] = InitializeTrainingData(fltPatternArray_Training);
intDimensions = length(fltNewTrainingData(1,:));
intLengthTraining = length(fltNewTrainingData(:,1));

fltThetaParams = zeros(intUniqueGroups, intDimensions) + 0.5;

structLRC = struct(...
    'fltNewTrainingData', fltNewTrainingData, ...
    'intDimensions', intDimensions, ...
    'intLengthTraining', intLengthTraining, ...
    'intUniqueGroups', intUniqueGroups, ...
    'fltThetaParams', fltThetaParams);

%% LRC Method
if(method == 1)
    [structLRC] = ComputeGradDescent(structLRC, grpSource_Training);
else
    [structLRC] = ComputeNewtonIteration(structLRC, grpSource_Training);
end

end