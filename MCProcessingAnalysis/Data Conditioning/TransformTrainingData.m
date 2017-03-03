function [fltPatternArrayMerge, canonicalWeights, mx, vectorMean, vectorStd, structRandom, fltPatternPhotometric] = ...
    TransformTrainingData(fltPatternArray_TrainCross, grpSource_TrainCross, structMCReduct, idxTrainCross, grpSource)



%% =====================================
% Extended CVA Applied to Training and Crossval (Total, with 5-Fold Validation)      
ECVAm = ecva('model', fltPatternArray_TrainCross, grpSource_TrainCross, 20, 'none', 'syst123', 5, [], 0);
fltPatternArrayMC = ECVAm.CanonicalVariates;
   
mx = ECVAm.Detail.mx; % To be used for predictions
canonicalWeights = ECVAm.CanonicalWeights; % To be used for predictions

%% Merge with photometric data
counter = 1;
fltPatternPhotometric = [];
for i = idxTrainCross

    if(~strcmp(grpSource{i}, 'No Variablity'))
        fltPatternPhotometric(counter,:) = structMCReduct(i).parameters;
        counter = counter + 1;
    else
        counter = counter + 1;
    end
end

% Generate Random for Photometric
structRandom = [];
for i = 1:1:length(fltPatternPhotometric(1,:))
    
    x_Bar = mean(fltPatternPhotometric(:,i));
    x_std = std(fltPatternPhotometric(:,i));
    
    tmpArray = fltPatternPhotometric(fltPatternPhotometric(:,i) > (x_Bar - 2*x_std) & ...
        fltPatternPhotometric(:,i) < (x_Bar + 2*x_std),i);
    
    [f,x] = ecdf(tmpArray);
    
    [C,ia,ic] = unique(x);
    
    structRandom(i).f = f(ia);
    structRandom(i).x = C;
    
end

%
counter = 1;
for i = idxTrainCross
    if(strcmp(grpSource{i}, 'No Variablity'))
        
        diffPhoto = zeros(1,length(structMCReduct(1).parameters));
        for j = 1:1:length(fltPatternPhotometric(1,:))
            xq = rand;
            diffPhoto(j) = interp1(structRandom(j).f,structRandom(j).x,xq);
        end
        
        fltPatternPhotometric(counter,:) = diffPhoto;
    end
    
    counter = counter + 1;
end

fltPatternArrayMerge = horzcat(fltPatternArrayMC, fltPatternPhotometric);

% standardize the merged dataset
vectorMean = mean(fltPatternArrayMerge,1);
vectorStd = std(fltPatternArrayMerge, 1);
for i = 1:1:length(fltPatternArrayMerge)
    
    fltPatternArrayMerge(i,:) = (fltPatternArrayMerge(i,:) - vectorMean)./vectorStd;
    
end


end