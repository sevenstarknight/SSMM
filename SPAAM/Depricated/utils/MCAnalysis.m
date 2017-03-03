function[structMCReduct] = MCAnalysis(structMCReduct, states)

%% States
for i = 1:1:length(structMCReduct)

    timeSeries = structMCReduct(i).timeSeries;

    [markovChain] = ConstructMCDiff(timeSeries, states);
    structMCReduct(i).MC = markovChain;
    structMCReduct(i).unpackMC = reshape(markovChain, [length(states)*length(states), 1]);
  
end


end
