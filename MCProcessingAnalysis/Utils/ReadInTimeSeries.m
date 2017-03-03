function [structMCReduct, grpSource] = ReadInTimeSeries(stringDirectory)

%% Location of file
cd(stringDirectory);
fileIDLinear = csvread('LINEARattributesFinalApr2013.csv');

idArray = {};
for i = 1:1:length(fileIDLinear)
    idArray{i} = [num2str(fileIDLinear(i,1)), '.dat'];
end


cd(strcat(stringDirectory,'\datFileFromLINEAR'));

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
    
    timeSeries = dlmread(idArray{ind});
    [n,m] = size(timeSeries);
    
    if(m ~= 3)
       m; 
    end
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


end