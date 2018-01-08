% =======================================================
% Copyright (c) 2005, Kyle Johnston
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% The views and conclusions contained in the software and documentation are those
% of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of the FreeBSD Project.
% =======================================================
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