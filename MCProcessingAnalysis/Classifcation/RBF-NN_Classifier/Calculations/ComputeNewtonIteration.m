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
function [structLRC] = ComputeNewtonIteration(structLRC, grpSource_Training)

%% Initialize Variables
tol = 1e-6;

fltNewTrainingData = structLRC.fltNewTrainingData;

intUniqueGroups = structLRC.intUniqueGroups;
intLengthTraining = structLRC.intLengthTraining;

pEst = zeros(intLengthTraining, intUniqueGroups);

x = structLRC.fltNewTrainingData;
y = InitializeYArray(structLRC,grpSource_Training);
w = eye(intLengthTraining,intLengthTraining);

counter = 0;
iter = 0;

convergenceRatio = length(grpSource_Training);
convergencePrior = convergenceRatio;

while( convergenceRatio > tol)

%     convergenceRatio
    
    %% Generate P Estimates
    for i = 1:1:intLengthTraining
        pEst(i,:) = ComputeSoftMax(structLRC, fltNewTrainingData(i,:));      
        
        yEst = y(i,:)*pEst(i,:)';
        w(i,i) = (yEst)*(1 - yEst);
    end

    crossEntropyArray = zeros(intLengthTraining,1);
    
    for i = 1:1:intLengthTraining
        logP = log(pEst(i,:));
        
        oneMinuslogP = log(1 - pEst(i,:));
        
        crossEntropyArray(i) = (y(i,:)*logP' + (1 - y(i,:))*oneMinuslogP');
    end
    
    crossEntropy = -sum(crossEntropyArray)/intLengthTraining;
     

    % insure symmetry
    invC = x'*w*x;
    invCtrans = invC';
    invC = triu(invC) + tril(invCtrans) - diag(diag(invC));

%     first = (invC)^-1;
%     second = x'*(y - pEst);
    
    %% Compute Iterative Step
    deltaTheta = (invC)\(x'*(y - pEst));
%     deltaTheta = first*second;

    %% Iterate Variable
    currentW = structLRC.fltThetaParams;
    currentW = currentW + deltaTheta';
    structLRC.fltThetaParams = currentW;

    delta = 0;
    for i = 1:1:intLengthTraining
       delta = delta + sum((y(i,:) - pEst(i,:)).^2); 
    end
    
    convergenceRatio = (delta - convergencePrior)^2;
    convergencePrior = delta;
    
    %% Graphics for LRC Convergence
    counter = counter + 1;
    
    if(mod(counter,1) == 0)
    
        iter = iter + 1;
        crossEntropyVsIterations(iter, 1) = counter;
        crossEntropyVsIterations(iter, 2) = crossEntropy;

        convergenceRatioVsIterations(iter,1) = counter;
        convergenceRatioVsIterations(iter,2) = convergenceRatio;

        subplot(2,1,1)
        semilogy(convergenceRatioVsIterations(:,1),convergenceRatioVsIterations(:,2),'*r')
        title('Convergence Ratio')

        subplot(2,1,2)
        semilogy(crossEntropyVsIterations(:,1),crossEntropyVsIterations(:,2),'*r')
        title('Cross Entropy')

        pause(0.5)
    end
    
end

end
