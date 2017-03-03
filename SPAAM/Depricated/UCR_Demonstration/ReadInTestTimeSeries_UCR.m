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