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
function [structLRC, nodes, spread, errorEst, classEstimate] = RadialBasisFunctionClassifier_KMeans(fltPatternArray_Training, grpSource_Training, fltPatternArray_CrossVal, grpSource_CrossVal, spread)

%% Initialize Parameters
intLengthTraining = length(fltPatternArray_Training(:,1));
uniqueGrp = unique(grpSource_Training);

structPsi_All = [];

for i = 1:1:length(uniqueGrp)
    TF = strcmp(grpSource_Training, uniqueGrp(i));
    
    fltLimitedSet = fltPatternArray_Training(TF,:);
    
    fltUniquePatterns = unique(fltLimitedSet, 'rows');
    
    if(length(fltUniquePatterns(:,1)) > 5)
        [structPsi] = K_Means(fltPatternArray_Training(TF,:)', 5);
    else
        [structPsi] = K_Means(fltPatternArray_Training(TF,:)', length(fltUniquePatterns(:,1)));
    end
    
    structPsi_All(i).structPsi = structPsi;
    
end

% Make Nodes in Hyper-Dimensions
[nodes] = MakeNodes(structPsi_All);
fltRejectionValue = 0;

%% RBF Method
% Translate into Radial Basis Function

% Training
fltRBData = zeros(intLengthTraining, length(nodes));

for i = 1:1:intLengthTraining
    fltPositionVector = fltPatternArray_Training(i,:);
    [fDistance] = GaussianKernel(fltPositionVector, nodes, spread);
    fltRBData(i, :) = fDistance;
end

method = 2;

[structLRC, uniqueGrps] = LogisticRegressionClassifier(fltRBData, grpSource_Training, method);

[errorEst, percentReject, classEstimate] = RBFCrossVal(fltPatternArray_CrossVal, grpSource_CrossVal, structLRC, uniqueGrps , nodes, spread, fltRejectionValue);

end

function [nodes] = MakeNodes(structPsi_All)

intClasses = length(structPsi_All);
count = 0;

for i = 1:1:intClasses
    intNodes = length(structPsi_All(i).structPsi);
    structPsi = structPsi_All(i).structPsi;
    for j = 1:1:intNodes
        count = count + 1;
        nodes{count} = structPsi(j).mean;
    end
end


end