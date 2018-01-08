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
function [structPatternArray_TrainingCV, indexTesting] =   ...
    Generate_5FoldCrossVal(grpSource)

intLength = length(grpSource);
uniqueTypes = unique(grpSource);
intNumTypes = length(uniqueTypes);

indexTrainingCrossVal = [];
indexTesting = [];

%% Pipe in each class into a separate 
for i = 1:1:intLength
    X = rand;
    
    for j = 1:1:intNumTypes
        if(strcmp(grpSource{i}, uniqueTypes{j}))
            %% Random sample into both training/cross-val and testing
            if(X > 0.5)
                if(isempty(indexTrainingCrossVal))
                    indexTrainingCrossVal(1) = i; 
                else
                    n = length(indexTrainingCrossVal);
                    indexTrainingCrossVal(n + 1) = i; 
                end
            else
                if(isempty(indexTesting))
                    indexTesting(1) = i;
                else
                    n = length(indexTesting);
                    indexTesting(n + 1) = i;
                end
            end
        end
    end
end

%% Randomly distribute into five sets
p = randperm(length(indexTrainingCrossVal));
sampleSize = round(length(p)/5) - 1;

structPatternArray_TrainingCV = [];
count = 1;
for j = 1:1:5
    indexSet = p(count:count + sampleSize - 1);
    count = count + sampleSize;
    structPatternArray_TrainingCV(j).indexSet = indexSet;
end


end