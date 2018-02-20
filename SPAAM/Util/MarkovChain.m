% ==================================================== 
%  Copyright (C) 2016 Kyle Johnston
%  
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ====================================================
function [markovChain] = MarkovChain(structSampleAreas, states)


%% Construct the MC
markovChain = zeros(length(states), length(states));
for i = 1:1:length(structSampleAreas)
    [markovChainNew] = ConstructMC(structSampleAreas(i).set, states);
    markovChain = markovChain + markovChainNew;
end

%% turn counts into probablities
totalInState = sum(markovChain, 2);
for i = 1:1:length(totalInState)
    if(totalInState(i) ~= 0.0)
        markovChain(i,:) = (markovChain(i,:)./totalInState(i));
    end
end


end


function [markovChain] = ConstructMC(seqX, states)

%% Construct the MC
markovChain = zeros(length(states), length(states));

%% Translation Array for the Time Series (Symbolic Representation)
stateSeqX = zeros(size(seqX));

%% Handle outside of state cases
stateSeqX(seqX < states(1)) = 1;
stateSeqX(seqX > max(states)) = length(states);

%% Handle inside of state cases
for i = 2:1:length(states)
    stateSeqX(seqX >= states(i - 1) & seqX < states(i)) = i - 1;
end

%% Populate Markov Model
for i = 2:1:length(seqX) 
    markovChain(stateSeqX(i - 1), stateSeqX(i)) = markovChain(stateSeqX(i - 1), stateSeqX(i)) + 1;
end

% markovChain_Old = zeros(length(states), length(states));
% for i = 2:1:length(seqX)
%     %CURRENT
%     [indexStart] = findIndex(seqX(i - 1), states);
% 
%     %NEXT
%     [indexStop] = findIndex(seqX(i), states);
% 
%     markovChain_Old(indexStart, indexStop) = markovChain_Old(indexStart, indexStop) + 1;
% end
% 
% 
% x = sum(markovChain - markovChain_Old);

end