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
function [structLRC, nodes, spread, errorEst, percentReject] = RadialBasisFunctionClassifier_EM(fltTraining, grpTraining, structPsi_All)

%% Initialize Parameters
intLengthTraining = length(fltTraining(:,1));
fltRejectionValue = 0; % if zero, everything is classified

% Make Nodes in Hyper-Dimensions
[nodes, spread] = MakeNodes(structPsi_All);

%% RBF Method
% Translate into Radial Basis Function

% Training
fltRBData = zeros(intLengthTraining, length(nodes));

for i = 1:1:intLengthTraining
    fltPositionVector = fltTraining(i,:);
    [fDistance] = GaussianKernel(fltPositionVector, nodes, spread);
    fltRBData(i, :) = fDistance;
end

method = 1;

gplotmatrix(fltRBData,fltRBData,grpTraining);

[structLRC] = LogisticRegressionClassifier(fltRBData, grpTraining, method);

errorEst = 0;
percentReject = 0;

% %% LOOCV
% errorArray = zeros(1,length(fltTraining));
% h = waitbar(0,'LOOCV in Progress') ;
% 
% for i = 1:1:length(fltTraining)
%     
%     ratio = double(i)/double(length(fltTraining));
%     waitbar(ratio)
% 
%     set = [];
%     for j = 1:1:length(fltTraining)
%         if(j ~= i)
%             set = cat(2,set,j);
%         end
%     end
%     
% [errorArray(i), percentReject] = RBFCrossVal(fltTraining(set,:), grpTraining(set), structLRC, nodes, spread, fltRejectionValue);
% 
% end
% 
% close(h) 
% 
% errorEst = mean(errorArray);
end

function [nodes, spread] = MakeNodes(structPsi_All)

intClasses = length(structPsi_All);
count = 0;

for i = 1:1:intClasses
    intNodes = length(structPsi_All(i).structPsi);
    structPsi = structPsi_All(i).structPsi;
    for j = 1:1:intNodes
        count = count + 1;
        nodes{count} = structPsi(j).mean;
        spread{count} = structPsi(j).cov;
    end
end


end
