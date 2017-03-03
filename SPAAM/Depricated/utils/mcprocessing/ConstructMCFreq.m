function [markovChain] = ConstructMCFreq(seqX, seqT, states)

%% Construct the MC
markovChain = zeros(length(states), length(states));

for(i = 2:1:length(seqX))
    
    % skip over the jumps in day
    if(seqT(i) < 40 && seqT(i-1) < 40)
        
        %CURRENT
        [indexStart] = findIndex(seqX(i - 1), states);
        
        %NEXT
        [indexStop] = findIndex(seqX(i), states);
        
        markovChain(indexStart, indexStop) = markovChain(indexStart, indexStop) + 1;
    end 
end


end