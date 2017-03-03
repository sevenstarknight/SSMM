function[markovChain] = ConstructMCSlope(timeSeries, states)

%% Diff
deltaTimeSet = diff(timeSeries(:,1));
deltaAmpSet = diff(timeSeries(:,2));

% Generate Slope Estimate
slope = deltaAmpSet./deltaTimeSet;

%% Standardize the slope
slope = (slope - mean(slope))/std(slope);

[markovChain] = ConstructMC(slope, deltaTimeSet, states);

end

