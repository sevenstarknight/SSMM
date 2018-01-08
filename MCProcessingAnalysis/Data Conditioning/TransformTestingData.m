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