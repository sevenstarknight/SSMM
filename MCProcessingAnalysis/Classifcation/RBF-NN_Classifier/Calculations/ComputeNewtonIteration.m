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
