function [markovChain] = ConstructMC(seqX, states)

%% Construct the MC
markovChain = zeros(length(states), length(states));

for i = 2:1:length(seqX)
    %CURRENT
    [indexStart] = findIndex(seqX(i - 1), states);

    %NEXT
    [indexStop] = findIndex(seqX(i), states);

    markovChain(indexStart, indexStop) = markovChain(indexStart, indexStop) + 1;
end

%% turn counts into probablities
totalInState = sum(markovChain, 2);
for(i = 1:1:length(totalInState))
    if(totalInState(i) ~= 0.0)
        markovChain(i,:) = (markovChain(i,:)./totalInState(i));
    end
end


end