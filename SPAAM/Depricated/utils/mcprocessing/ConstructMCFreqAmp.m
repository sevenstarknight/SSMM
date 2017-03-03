function[markovChain, meanAmplitude, stdAmplitude] = ConstructMCFreqAmp(timeSeries, states)

%% Diff
timeSet = diff(timeSeries(:,1));
ampSet = (timeSeries(:,2));

%% Standardize the amplitude

meanAmplitude = mean(ampSet);
stdAmplitude = std(ampSet);

ampSet = (ampSet - meanAmplitude)/stdAmplitude;

[markovChain] = ConstructMCFreq(ampSet(2:end), timeSet, states);

end