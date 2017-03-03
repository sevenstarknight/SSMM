function [structLRC] = ComputeGradDescent(structLRC, grpTraining)

%% Initialize Variables
tol = 1e-6;
stepSize = 0.00001;

iter = 0;
counter = 0;

deltaCrossEntropy = 1;
crossEntropyPrior = 1;


while (deltaCrossEntropy > tol)
    
    %% Compute Iterative Cross-Entropy
    [deltaW, ~, crossEntropy] = ComputeDelEDelWIteration(structLRC, grpTraining);
    
    deltaCrossEntropy = abs(crossEntropyPrior - crossEntropy)/crossEntropy;
    
    crossEntropyPrior = crossEntropy;

    % In and out for Theta Params
    fltThetaParams = structLRC.fltThetaParams;
    fltThetaParams = fltThetaParams - stepSize*deltaW;
    structLRC.fltThetaParams = fltThetaParams;
    
    %% Graphics for LRC Convergence
    counter = counter + 1;
    
    if(mod(counter,100) == 0)
    
        iter = iter + 1;

        crossEntropyVsIterations(iter,1) = iter;
        crossEntropyVsIterations(iter,2) = crossEntropy;

        convergenceRatioVsIterations(iter,1) = iter;
        convergenceRatioVsIterations(iter,2) = deltaCrossEntropy;

        subplot(2,1,1)
        semilogy(convergenceRatioVsIterations(:,1),convergenceRatioVsIterations(:,2),'*r')
        title('Delta Cross-Entropy vs. Iteration')
        ylabel('Delta Cross-Entropy')

        subplot(2,1,2)
        semilogy(crossEntropyVsIterations(:,1),crossEntropyVsIterations(:,2),'*r')
        title('Cross-Entropy vs. Iteration')
        ylabel('Cross-Entropy')

        pause(0.5)
    end
    
end


end

function [deltaW, delta, crossEntropy] = ComputeDelEDelWIteration(structLRC, grpTraining)

%% Initialization
intUniqueGroups = structLRC.intUniqueGroups;
intLengthTraining = structLRC.intLengthTraining;
fltNewTrainingData = structLRC.fltNewTrainingData;

delta = 0;

y = InitializeYArray(structLRC,grpTraining);
pEst = zeros(intLengthTraining, intUniqueGroups);

%% Compute pEst
for i = 1:1:intLengthTraining
    pEst(i,:) = ComputeSoftMax(structLRC, fltNewTrainingData(i,:));
end    

crossEntropyArray = zeros(intLengthTraining,1);

for i = 1:1:intLengthTraining
    logP = log(pEst(i,:));
    crossEntropyArray(i) = y(i,:)*logP';
end

crossEntropy = -sum(crossEntropyArray);

errorTerm = pEst - y;

deltaW = errorTerm'*fltNewTrainingData;

%% Compute Error Estimate
for i = 1:1:intLengthTraining
   delta = delta + sum((y(i,:) - pEst(i,:)).^2); 
end

end
