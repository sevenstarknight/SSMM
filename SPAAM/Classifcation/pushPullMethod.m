function [M, pullTime, pushTime] = pushPullMethod(trainConvDF, trainLabels, gamma, lambda)

%% Split Matrixes of different classes. 
uniqueClasses = unique(trainLabels);
convDFCell = cell(length(uniqueClasses),1);
for c = 1 : length(uniqueClasses)
    convDFCell{c,1} = trainConvDF(:,:,strcmp(uniqueClasses(c),trainLabels));
end

%% Pull Terms
[B, ~, N] = size(trainConvDF);
totalPull = zeros(B);
tic
for c = 1 : length(uniqueClasses)
    N_c = size(convDFCell{c,1},3)/2;
    tempMatrix = convDFCell{c,1};
    for i = 1 : N_c
        tmpSubMatrix = tempMatrix(:,:,i);
        for j = 1 : N_c
            temp = tmpSubMatrix - tempMatrix(:,:,j);
            totalPull = totalPull + (temp*temp');
        end
    end
    totalPull = totalPull/(N_c - 1);
end
totalPull = -totalPull/gamma;
pullTime = toc;

%% Computing the Push terms 
tic
totalPush = zeros(B);
for c = 1 : length(uniqueClasses)
    total = zeros(B);
    N_c = size(convDFCell{c,1},3);
    tempMatrixC = convDFCell{c,1};
    for i = 1 : N_c
        tempMatrixCI = tempMatrixC(:,:,i);
        for k = 1 : length(uniqueClasses)
            if (c ~= k) % not the same class
                N_k = size(convDFCell{k,1},3);
                tempMatrixK = convDFCell{k,1};
                for j = 1 : N_k
                    temp = tempMatrixCI - tempMatrixK(:,:,j);
                    total = total + (temp*temp');
                end
            end
        end
    end
    totalPush = totalPush + (total/(N-N_c));
end
totalPush = totalPush/gamma;
pushTime = toc;

%% Process
M = (1-lambda)*totalPull + lambda*totalPush;

%% Project to PSD Space
[U, D] = eig(M);
D(D<0)=0;
M = U*D*U'; 
end
