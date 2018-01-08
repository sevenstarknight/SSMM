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
function [structQDAParam] = GenerateQDAParams(fltPatternArray_Training,...
    grpSource_Training, switch_expr)
% function [structQDAParam] = GenerateQDAParams(fltPatternArray_Training, grpSource_Training,
% switch_expr)
%
% Author: Kyle Johnston     14 Jan 2010
%
% Usage: This function generates class estimates based on the QDA
%   algorithm. Ties are represented by the "missValue" string in the
%   output labels. 
% 
% Input: 
%       fltPatternArray_Training: the training dataset
%       grpSource_Training: the labels for the training dataset
%       crit_alpha: int representing the kernel selected

% Output:
%       structQDAParam: distribution parameters 
%% Initialize Parameters
uniqueGrp = unique(grpSource_Training);
intNumberGroups = length(uniqueGrp);

structQDAParam= struct([]);

% Estimate and store parameters
for i = 1:1:intNumberGroups
    boolDecision =  strcmp(uniqueGrp{i},grpSource_Training);
    
    fltReducedSet = fltPatternArray_Training(boolDecision,:);
    
    structQDAParam(i).mean = mean(fltReducedSet, 1);
    
    % Posterior Probabilities
    structQDAParam(i).n = sum(boolDecision);
    
    % Select method for generation of covariance matrix
    if(switch_expr == 1)
        [covMatrix] = GenerateGeneralCaseQDA(fltReducedSet);
    elseif(switch_expr == 2)
        [covMatrix] = GenerateNaiveCaseQDA(fltReducedSet);
    elseif(switch_expr == 3)
        [covMatrix] = GenerateIsotropicQDA(fltReducedSet);
    else
        covMatrix = [];
    end
    
    structQDAParam(i).cov = covMatrix;
    structQDAParam(i).invCov = inv(structQDAParam(i).cov);
    structQDAParam(i).logDet = log(det(structQDAParam(i).cov));
    
    structQDAParam(i).constant = log(structQDAParam(i).n/length(grpSource_Training)) - 0.5*structQDAParam(i).logDet;
    structQDAParam(i).type = uniqueGrp(i);
    
end

end


%% Subfunction General QDA
function [covMatrix] = GenerateGeneralCaseQDA(fltReducedSet)

    covMatrix = cov(fltReducedSet);

end

%% Subfunction Naive QDA
function [covMatrix] = GenerateNaiveCaseQDA(fltReducedSet)

    dimen = length(fltReducedSet(1,:));
    covMatrix = zeros(dimen);
    
    for i = 1:1:dimen
        covMatrix(i,i) = var(fltReducedSet(:,i));
    end

end

%% Subfunction Isotropic QDA
function [covMatrix] = GenerateIsotropicQDA(fltReducedSet)

    dimen = length(fltReducedSet(1,:));
    varMatrix = zeros(1,dimen);
    
    for i = 1:1:dimen
        varMatrix(i) = var(fltReducedSet(:,i));
    end

    estVar = mean(varMatrix);
    covMatrix = eye(dimen)*estVar;
end