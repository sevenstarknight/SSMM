%-------------------------------------------------------------------------
% Function: DF_KNN
%
% Description: This function takes in training and testing curves including
%              their respective labels. Computes the distances between
%              testing curves and all training curves. Then takes the k
%              closest ones and picks the label based on a majority vote of
%              the k nearest neighbors. Once it assigns the predicted
%              labels, it compares it to the true labels and computes an
%              accuracy measure. 
%
% Inputs:
%     trainConvDF   - B x M x N matrix of the distribution fields of the
%                     training curves. B is the number of bins and M is the
%                     length of the curves and N is the number of curves.
%                     OR it can be theN x M matrix of the training data
%                     curves.
%     trainLabels   - M x 1 vector that has the class labels for Training
%                     Curves.
%     testConvDF    - B x M x N matrix of the distribution fields of the
%                     testing curves. B is the number of bins and M is the
%                     length of the curves and N is the number of curves.
%                     OR it can be theN x M matrix of the testing data
%                     curves.
%     testLabels    - M x 1 vector that has the class labels for Testing
%                     Curves.
%     options       - Any Options required such as:
%                           - k, Value for kNN
%                           - M, B x B metric matrix
%
% Outputs:
%     predictLabel  - M x 1 vector that has the predicted class labels for 
%                     the Testing curves.
%     accuracy      - Compares the predictLabel and testing true labels and
%                     then computes: (# correct labels)/(total labels), the
%                     value returned is between 0 and 1, where 1 is a 100%
%                     accuracy.
%
% Usage:
%   1 - After computing the DF of the curves
%   2 - set options.kValue 
%   3 - call [predict_Label_DF, accuracy_DF] = ...
%             DF_KNN (trainConvDF, trainLabels, testConvDF,...
%                                 testLabels, options);
%
% Update(s):
%   DEC2015 -- Rana Haber
%              Initial Code
%   JAN2016 -- Rana Haber
%              Updated the distance equation to include M (the metric)
%-------------------------------------------------------------------------
function [predictLabel, accuracy] = DF_KNN (trainConvDF, ...
                                            trainLabels, ...
                                            testConvDF, ...
                                            testLabels, ...
                                            options)
kValue = options.kValue;

[~, ~, curvesNum_train] = size(trainConvDF);
[~, ~, curvesNum_test] = size(testConvDF);

if (curvesNum_test > 1)

    metric = options.M;
	predictLabel = zeros(curvesNum_test,1);

	parfor n = 1 : curvesNum_test
        tempDist = zeros(1,curvesNum_train);
        temp1 = testConvDF(:,:,n);
        for m = 1 : curvesNum_train
            temp = temp1-trainConvDF(:,:,m);
            tempDist(1,m) = sqrt(sum(diag(temp'*metric*temp)));
        end

        [~, index ] = sort(tempDist);
        distLabels = trainLabels(index(1:kValue));
        [~,predictLabel(n)] = max(histc(distLabels,unique(trainLabels)));
	end
else
    [curvesNum_test, ~] = size(testConvDF);
    
    predictLabel = zeros(curvesNum_test,1);
    dist = pdist2(testConvDF, trainConvDF, 'euclidean');
    
    [~, index] = sort(dist,2);
    parfor n = 1 :  curvesNum_test
        distLabels = trainLabels(index(n,1:kValue));
        [~,predictLabel(n)] = max(histc(distLabels,unique(trainLabels)));
    end
end

accuracy = sum(predictLabel == testLabels)/ curvesNum_test;
end
