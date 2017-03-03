clc
clear

%% ====================================================
% Make feature space and load data
resArray = 0.02:0.01:0.3;
errorArray = zeros(1,length(resArray));

[structMCReduct, grpSource] = ReadInTimeSeries();

[structPatternArray_TrainingCV, indexTesting] = ...
    Generate_5FoldCrossVal(grpSource);


index = 1;

timeSeries = structMCReduct(index).timeSeries;

% construct state space, for given resolution.
states = -2:resArray(index):2;

%% Diff
deltaTimeSet = diff(timeSeries(:,1));
deltaAmpSet = diff(timeSeries(:,2));

%% Standardize prior to the diff
seqT = deltaTimeSet;
seqX = (deltaAmpSet - mean(deltaAmpSet))/std(deltaAmpSet);
seq = [];

count = 1;
for(i = 1:1:length(seqX))
    % skip over the jumps in day
    if(seqT(i) < 40)
        %CURRENT
        index = findIndex(seqX(i), states);
        
        if(index == 0)
            seq(count) = 1;
            count = count + 1;
        else
            seq(count) = index;
            count = count + 1;
        end
    end 
end


[TRANS,EMIS] = hmmestimate(seq',ones(size(seq')));

