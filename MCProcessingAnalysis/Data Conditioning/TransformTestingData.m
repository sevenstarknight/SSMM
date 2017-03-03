function [fltPatternArray_Testing, grpSource_Testing] = ...
    TransformTestingData(spaceTotal, grpSource, indexTesting, ...
    canonicalWeights, mx, vectorMean, vectorStd, structRandom, structMCReduct, fltPatternPhotometric)


%% ==================================================================
% Grab Testing Data
fltPatternArray_Testing = spaceTotal(indexTesting, :);
grpSource_Testing = grpSource(indexTesting);

% transform testing data
for i = 1:1:length(grpSource_Testing)
    ECVATesting(i, :) = (fltPatternArray_Testing(i, :) - mx);
end
fltPatternArray_Testing_Canonical = ECVATesting*canonicalWeights;

counter = 1;
fltPatternPhotometricTesting = [];
for i = indexTesting

     if(strcmp(grpSource{i}, 'No Variablity'))
        
        diffPhoto = zeros(1,length(structMCReduct(1).parameters));
        for j = 1:1:length(fltPatternPhotometric(1,:))
            xq = rand;
            diffPhoto(j) = interp1(structRandom(j).f,structRandom(j).x,xq);
        end
        
        fltPatternPhotometricTesting(counter,:) = diffPhoto;
     else
        fltPatternPhotometricTesting(counter,:) = structMCReduct(i).parameters;

     end
     counter = counter + 1;

end

fltPatternArray_Testing = horzcat(fltPatternArray_Testing_Canonical, fltPatternPhotometricTesting);

for i = 1:1:length(fltPatternArray_Testing)
    fltPatternArray_Testing(i,:) = (fltPatternArray_Testing(i,:) - vectorMean)./vectorStd;
end
grpSource_Testing = grpSource(indexTesting);


end