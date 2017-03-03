function [fltProbEst] = ComputeSoftMax(structLRC, fltPatternArray)
% function [fltProbEst] = ComputeSoftMax(structLRC, fltPatternArray)
%
% Author: Kyle Johnston     14 Jan 2010
%
% Usage: This function computes the postierior probability of the target
% (fltPatternArray) being of type clutter or SOI, given the set of linear
% coefficient provided in the structure structLRC. Note, this code is set
% up to allow for multi-class hypothesis classification, currently it is
% only implimented in a detection-like set up (intUniqueGroups = 2). 
%
% Meaningful improvements to the linear coefficients cannot be made by
% hand, and can only be done inlab using the LRC classifier. LRC
% coefficents here have been computed using "optimal" conditions, i.e.
% minimization of error and false alarm = missed detection
% 
% Input: 
%       structLRC
%       fltPatternArray
% Output:
%       fltProbEst: array of postieror probabilities

%% Initialize
intUniqueGroups = structLRC.intUniqueGroups;
fltThetaParams = structLRC.fltThetaParams;

fltProbEst = zeros(1,length(intUniqueGroups));

%% Compute Soft-Max Probability
for i = 1:1:intUniqueGroups
    
    w_i = fltThetaParams(i,:);
    
    a = w_i*fltPatternArray';
    fltNumer = exp(a);
    
    fltDenom = 0;
    for j = 1:1:intUniqueGroups
        w_j = fltThetaParams(j,:);
        b = w_j*fltPatternArray';
        
        fltDenom = fltDenom + exp(b);
    end
    
    fltProbEst(i) = fltNumer/(fltDenom);
end

end