function[markovChain, meanAmplitude, stdAmplitude] = ConstructMCDiff(timeSeries, states)

%% Diff
deltaTimeSet = diff(timeSeries(:,1));
deltaAmpSet = diff(timeSeries(:,2));

%% Standardize prior to the diff
meanAmplitude = mean(deltaAmpSet);
stdAmplitude = std(deltaAmpSet);
deltaAmpSet = (deltaAmpSet - meanAmplitude)/stdAmplitude;

[markovChain] = ConstructMC(deltaAmpSet, deltaTimeSet, states);

end
