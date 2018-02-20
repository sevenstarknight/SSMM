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
function [structTraining, grpSource_Training, structTesting, grpSource_Testing] = ReadInUCRTimeSeries(stringDirectory, strLabel)

cd(strcat(stringDirectory,'\',strLabel));

%% Location of file
TRAIN = load(strcat(strLabel,'_TRAIN')); % Only these two lines need to be changed to test a different dataset. %
TEST  = load(strcat(strLabel,'_TEST' )); % Only these two lines need to be changed to test a different dataset. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRAIN_class_labels = TRAIN(:,1);     % Pull out the class labels.
TEST_class_labels = TEST(:,1);       % Pull out the class labels.

structTraining = [];
grpSource_Training = {};

structTesting = [];
grpSource_Testing = {};

alphabet = char('a' +(1:26)-1)';
[I,J] = ndgrid(1:26,1:26);
I=I';J=J';
XX = [alphabet(I(:)), alphabet(J(:))];
% XX = strvcat(alphabet, XX);
alphabet = XX;

alphabetCell = cell(1,length(alphabet));
for(i = 1:1:length(alphabet))
    
   alphabetCell{i} = alphabet(i,:); 
    
end

disp(strcat(num2str(length(TRAIN_class_labels)) , ': ', 'Training Labels'))
%% Add in the 
for i = 1:1:length(TRAIN_class_labels)  
       
    structTraining(i).ID = i;
    structTraining(i).parameters = []; 
    structTraining(i).classType =  TRAIN_class_labels(i);
    
    maxTime = 0.05*length(TRAIN(i, 2:end));
    
    structTraining(i).timeSeries = [0:0.05:(maxTime-0.05); TRAIN(i, 2:end)]';
    
    grpSource_Training{i} = num2str(TRAIN_class_labels(i));
end
uniqTraining = unique(grpSource_Training);
for i = 1:1:length(TRAIN_class_labels)
    
    idx = find(ismember(uniqTraining,grpSource_Training{i}));
    grpSource_Training{i} = alphabetCell{idx};
end

%% Add in the 
for i = 1:1:length(TEST_class_labels)  
       
    structTesting(i).ID = i;
    structTesting(i).parameters = []; 
    structTesting(i).classType =  TEST_class_labels(i);
    
    maxTime = 0.05*length(TEST(i, 2:end));
    
    structTesting(i).timeSeries = [0:0.05:(maxTime-0.05); TEST(i, 2:end)]';
    
    grpSource_Testing{i} = num2str(TEST_class_labels(i));
end

uniqTesting = unique(grpSource_Testing);
for i = 1:1:length(TEST_class_labels)
    
    idx = find(ismember(uniqTesting,grpSource_Testing{i}));
    grpSource_Testing{i} = alphabetCell{idx};
end
end