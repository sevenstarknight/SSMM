function [structTraining, grpSource_Training, structTesting, grpSource_Testing] = ReadInLINEARTimeSeries(stringDirectory)

%% Location of file
cd('C:\Users\kyle.johnston\Desktop\Dissertation\starvars\data\trainingsets');
fileIDLinear = csvread('LINEARattributesFinalApr2013.csv');

idArray = {};
for i = 1:1:length(fileIDLinear)
    idArray{i} = [num2str(fileIDLinear(i,1)), '.dat'];
end

cd('C:\Users\kyle.johnston\Desktop\Dissertation\Old\Data\datFileFromLINEAR');

listing = dir;
grpSource = {};
structMCReduct = [];
count = 1;
%% Add in the 
% LINEARobjectID, ra, dec, ug, gi, iK, JK, logP, Ampl, skew, kurt, magMed, nObs, LCtype 
for i = 1:1:length(listing)  
       
    if(~ismember(listing(i).name,idArray))
        continue
    end
    
    ind = find(strcmp(listing(i).name, idArray));
    
    structMCReduct(count).ID = idArray{ind};
    structMCReduct(count).parameters = fileIDLinear(ind, 4:12); 
    structMCReduct(count).classType =  fileIDLinear(ind, 14);
    
    timeSeries = dlmread(idArray{ind}, '\t');
    structMCReduct(count).timeSeries = timeSeries;
    
    if(structMCReduct(count).classType == 1)
        grpSource{count} = 'RRab';
    elseif (structMCReduct(count).classType == 2)
        grpSource{count} = 'RRc';
    elseif (structMCReduct(count).classType == 4)
        grpSource{count} = 'Algol';
    elseif (structMCReduct(count).classType == 5)
        grpSource{count} = 'Contact Binary';
    elseif (structMCReduct(count).classType == 6)
        grpSource{count} = 'Delta Scu / SX Phe';
    else
        grpSource{count} = 'Unknown';
    end
    
     count = count + 1;
end

%% Divide into Training and Testing
intLength = length(grpSource);
uniqueTypes = unique(grpSource);
intNumTypes = length(uniqueTypes);

indexTrainingCrossVal = [];
indexTesting = [];

for i = 1:1:intLength
    X = rand;
    
    for j = 1:1:intNumTypes
        if(strcmp(grpSource{i}, uniqueTypes{j}))
            %% Random sample into both training/cross-val and testing
            if(X > 0.25)
                if(isempty(indexTrainingCrossVal))
                    indexTrainingCrossVal(1) = i; 
                else
                    n = length(indexTrainingCrossVal);
                    indexTrainingCrossVal(n + 1) = i; 
                end
            else
                if(isempty(indexTesting))
                    indexTesting(1) = i;
                else
                    n = length(indexTesting);
                    indexTesting(n + 1) = i;
                end
            end
        end
    end
end

%% Parse and Store
structTraining = [];
grpSource_Training = {};
structTesting = [];
grpSource_Testing = {};

for idx = 1:1:length(indexTrainingCrossVal)
    
    structTraining(idx) = structMCReduct(indexTrainingCrossVal(idx));
    grpSource_Training{idx} = grpSource{indexTrainingCrossVal(idx)};
    
end

for idx = 1:1:length(indexTesting)
    
    structTesting(idx) = structMCReduct(indexTesting(idx));
    grpSource_Testing{idx} = grpSource{indexTesting(idx)};
    
end


end