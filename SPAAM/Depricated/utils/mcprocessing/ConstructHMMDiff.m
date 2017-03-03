function[TRANS_EST] = ConstructHMMDiff(timeSeries, states)

timeSeries(:,2) = (timeSeries(:,2) - mean(timeSeries(:,2)))/std(timeSeries(:,2));

deltaTimeSet = diff(timeSeries(:,1));
deltaAmpSet = diff(timeSeries(:,2));

seq = zeros(1,length(deltaTimeSet));

for(i = 1:1:length(deltaTimeSet))
    
    if(deltaAmpSet(i) >= states(length(states)))
            indexStart = length(states);            
    else
        indexStart = find(deltaAmpSet(i) < states, 1) - 1;

        % bound the state search
        if(indexStart <= 0) 
            indexStart = 1;
        elseif(indexStart > length(states))
            indexStart = length(states);
        end
    end
    
    seq(i) = states(indexStart);
    
    
end


[TRANS_EST, EMIS_EST] = hmmestimate(deltaAmpSet, seq);

% R = poissrnd(1,size(markovChain));
% markovChain = markovChain + R;

% %% turn counts into probablities
% totalInState = sum(markovChain, 2);
% for(i = 1:1:length(totalInState))
%     if(totalInState(i) ~= 0.0)
%         markovChain(i,:) = (markovChain(i,:)./totalInState(i));
%     end
% end

% figure()
% image(markovChain)

end

function [A, B] = GenerateLine(x, y)

N = length(x);
ssqX = sum(x.^2);
sxy = sum(x.*y);
sx = sum(x);
sy = sum(y);

delta = N*ssqX - (sx^2);

% y = A + Bx
A = (ssqX*sy - sx*sxy)/delta;
B = (N*sxy - sx*sy)/delta;

end
