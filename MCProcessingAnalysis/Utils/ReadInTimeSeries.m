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