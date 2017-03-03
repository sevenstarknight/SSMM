function [structMCNew, grpSourceNew] = ConstructNonVariableStars(n, structMCReduct, grpSource)

class = 6;

for indexi = 1:1:n

    timeSeriesSize = randi([150, 300]);

    
    T = 0.01;                   % Sample time
    L = timeSeriesSize;           % Length of signal
    time = (0:L-1)*T;                % Time vector
    
    
    amplitudes = randn(1, timeSeriesSize)';
    error = zeros(size(amplitudes));
    time = time + T*0.1*randn(1,length(time));
    
    
    timeSeries = [time', amplitudes, error];

    indexNew = length(structMCReduct) + 1;

    grpSource{indexNew} = 'No Variablity';
    
    structMCReduct(indexNew).ID = strcat('Synth ', num2str(indexi));
    structMCReduct(indexNew).parameters = [];
    structMCReduct(indexNew).class = class;
    structMCReduct(indexNew).timeSeries = timeSeries;

    structMCReduct(indexNew).MC = [];
    structMCReduct(indexNew).unpackMC = [];

end

% copy
structMCNew = structMCReduct;
grpSourceNew = grpSource;
% permutate
p = randperm(length(structMCReduct));

for indexi = 1:1:length(structMCReduct)
    
   structMCNew(p(indexi)) = structMCReduct(indexi); 
   grpSourceNew{p(indexi)} = grpSource{indexi}; 
end


end