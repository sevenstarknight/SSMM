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
function [fltTrainSet, grpSourceTraining, fltTestingSet, grpSourceTesting] = ReadInTestTimeSeries_UCR(strLabel)

%% Location of file
TRAIN = load(strcat(strLabel, '_TRAIN')); % Only these two lines need to be changed to test a different dataset. %
TEST  = load(strcat(strLabel, '_TEST' )); % Only these two lines need to be changed to test a different dataset. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRAIN_class_labels = TRAIN(:,1);     % Pull out the class labels.
TEST_class_labels = TEST(:,1);       % Pull out the class labels.

fltTrainSet = TRAIN(:, 2:end);
fltTestingSet = TEST(:, 2:end);

grpSourceTraining = cell(1,length(TRAIN_class_labels));
grpSourceTesting = cell(1,length(TEST_class_labels));

%% Add in the 
for i = 1:1:length(TRAIN_class_labels)  
    grpSourceTraining{i} = num2str(TRAIN_class_labels(i));
end


for i = 1:1:length(TEST_class_labels)  
    grpSourceTesting{i} = num2str(TEST_class_labels(i));
end


end