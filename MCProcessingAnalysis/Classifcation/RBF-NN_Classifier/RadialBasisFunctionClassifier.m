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
function [structLRC, nodes, spread, errorEst, classEstimate] = ...
    RadialBasisFunctionClassifier(fltPatternArray_Training, ...
    grpSource_Training, fltPatternArray_CrossVal, grpSource_CrossVal, ...
    spread)

%% Initialize Parameters
intLengthTraining = length(fltPatternArray_Training(:,1));

% Make Nodes in Hyper-Dimensions
[nodes] = MakeNodes(fltPatternArray_Training);
fltRejectionValue = 0;

%% RBF Method
% Training
fltRBData = zeros(intLengthTraining, length(nodes));

for i = 1:1:intLengthTraining
    fltPositionVector = fltPatternArray_Training(i,:);
    
    for j = 1:1:length(nodes)
        centerEstimate = nodes{j};
        distance = fltPositionVector - centerEstimate;
        fDistance = exp(-spread*norm(distance));
    
%     [fDistance] = GaussianKernel(fltPositionVector, nodes, spread);
        fltRBData(i, j) = fDistance;
    end 
end

method = 2;

[structLRC, uniqueGrps] = LogisticRegressionClassifier(fltRBData, grpSource_Training, method);

[errorEst, percentReject, classEstimate] = RBFCrossVal(fltPatternArray_CrossVal, grpSource_CrossVal, structLRC, uniqueGrps , nodes, spread, fltRejectionValue);


end

function [nodes] = MakeNodes(fltPatternArray_Training)

intLengthTraining = length(fltPatternArray_Training(:,1));
nodes = cell(1,intLengthTraining);

for i = 1:1:intLengthTraining
    nodes{i} = fltPatternArray_Training(i,:);
end

end
