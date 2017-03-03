function[markovChain, meanAmplitude, stdAmplitude] = ConstructMCAmp(timeSeries, states)

%% Diff
timeSet = diff(timeSeries(:,1));
ampSet = (timeSeries(:,2));

%% Standardize the amplitude

meanAmplitude = mean(ampSet);
stdAmplitude = std(ampSet);

ampSet = (ampSet - meanAmplitude)/stdAmplitude;

[markovChain] = ConstructMC(ampSet(2:end), timeSet, states);

end