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
function [structMCReduct, grpSource] = ReadInTestTimeSeries(stringDirectory)


cd(strcat(stringDirectory,'\StarLightCurves'))

%% Location of file
TRAIN = load('StarLightCurves_TRAIN'); % Only these two lines need to be changed to test a different dataset. %
TEST  = load('StarLightCurves_TEST' ); % Only these two lines need to be changed to test a different dataset. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRAIN_class_labels = TRAIN(:,1);     % Pull out the class labels.
TEST_class_labels = TEST(:,1);       % Pull out the class labels.

intSource = vertcat(TRAIN(:,1), TEST(:,1));
fltDataSet = vertcat(TRAIN(:, 2:end), TEST(:, 2:end));

structMCReduct = [];
grpSource = {};

%% Add in the 
for i = 1:1:length(intSource)  
       
    structMCReduct(i).ID = i;
    structMCReduct(i).parameters = []; 
    structMCReduct(i).classType =  intSource(i);
    
    maxTime = 0.05*length(fltDataSet(i,:));
    
    structMCReduct(i).timeSeries = [0:0.05:(maxTime-0.05); fltDataSet(i,:)]';
    
    if(structMCReduct(i).classType == 1)
        grpSource{i} = '1';
    elseif (structMCReduct(i).classType == 2)
        grpSource{i} = '2';
    elseif (structMCReduct(i).classType == 3)
        grpSource{i} = '3';
    else
        grpSource{i} = 'Unknown';
    end
end


end