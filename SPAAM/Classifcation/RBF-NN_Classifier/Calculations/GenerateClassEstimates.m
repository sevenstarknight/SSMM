function [classEstimate] = GenerateClassEstimates(fltPatternArray, uniqueGrps, structLRC, fltRejectionValue)

intLengthArray = length(fltPatternArray(:,1));
classEstimate = cell(intLengthArray,1);

% Transform
for i = 1:1:intLengthArray
    fltPattern = fltPatternArray(i,:);
    [fltProbEst] = ComputeSoftMax(structLRC, fltPattern);
    
    if nargin == 3,
        index = max(fltProbEst) == fltProbEst;
        classEstimate{i} = uniqueGrps{index};
    else
        if(fltRejectionValue > max(fltProbEst) || length(max(fltProbEst)) > 1 || sum(max(fltProbEst) == fltProbEst) == 0)
            classEstimate{i} = 'Reject';
        else
            index = max(fltProbEst) == fltProbEst;
            classEstimate{i} = uniqueGrps{index};
        end
    end
    
end


end