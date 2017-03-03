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
