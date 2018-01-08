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
function [fltResponse, classEstimate] = Parzen_Window_Classifier(...
    fltPatternArray_Training, grpSource_Training, fltPatternArray_CrossVal, ...
    intInterpType, fltSpread, strMissedClass)
% function [fltResponse, classEstimate] = Parzen_Window_Classifier(...
%     fltPatternArray_Training, grpSource_Training, fltPatternArray_CrossVal, ...
%     intInterpType, fltSpread, strMissedClass)
%
% Author: Kyle Johnston     14 Jan 2010
%
% Usage: This function generates class estimates based on the PWC
%   algorithm. Ties are represented by the "missValue" string in the
%   output labels. This particular implimentation is a weighted PWC
%   algorithm based on the distance (given the PDF and Kernel shape from 
%   the input pattern to other training pattterns.
% 
% Input: 
%       fltPatternArray_TrainingSet: the training dataset
%       grpSource_TrainingSet: the labels for the training dataset
%       fltTestSet: the test dataset
%       intInterpType: int representing the kernel selected
%       fltSpread: for the kernel
%       missValue: strnig label for rejected points (ties)
% Output:
%       fltResponse: the set of posterior probabilities
%       classEstimate: estimated labels for fltTestSet
tic
%% Initialize Variables
grpUniqueTypes = unique(grpSource_Training);

fltPriorProb = zeros(1,length(grpUniqueTypes));
fltResponse = zeros(length(fltPatternArray_CrossVal), length(grpUniqueTypes));

classEstimate = cell(1,length(fltPatternArray_CrossVal));

intD = length(fltPatternArray_Training(1,:));

for i = 1:1:length(grpUniqueTypes)
    fltPriorProb(i) = sum(strcmp(grpUniqueTypes(i), grpSource_Training));  
    fltPriorProb(i) = fltPriorProb(i)*(fltSpread^intD);
    fltPriorProb(i) = 1/fltPriorProb(i);
end

%% Loop through testing data
for k = 1:1:length(fltPatternArray_CrossVal)
    
    %% Initialize PWC Summation
    kernEstArray = zeros(1, length(grpUniqueTypes));
  
    for i = 1:1:length(fltPatternArray_Training)
        
        x = (fltPatternArray_Training(i,:) - fltPatternArray_CrossVal(k,:))./fltSpread;
        [est] = KernelEstimate(intInterpType, x);
        index = find(strcmp(grpUniqueTypes,grpSource_Training(i)));
        kernEstArray(index) = kernEstArray(index) + prod(est); 
        
    end
    
    %% Select Class Type
    fltResponse(k,:) = kernEstArray.*fltPriorProb;
    fltSumProbs = sum(fltResponse(k,:));
    
    if(fltSumProbs ~= 0)
        
        fltResponse(k,:) = fltResponse(k,:)./fltSumProbs;
        fltLikelyType = max(fltResponse(k,:));
        classEstimate{k} = grpUniqueTypes{fltLikelyType == fltResponse(k,:)};
        
    else
        
        fltResponse(k,:) = 0;
        classEstimate{k} = strMissedClass;
        
    end
    
end

toc
end

function [est] = KernelEstimate(intInterpType, x)
switch intInterpType
    case 1
        [est] = uniformFunction(x);
    case 2
        [est] = triangleFunction(x);
    case 3
        [est] = gaussianFunction(x);
    case 4
        [est] = laplaceFunction(x);
    case 5
        [est] = cauchyFunction(x);
    case 6
        [est] = sinhFunction(x);
    otherwise
        disp('Method is unknown')
end

end

function [fltModelReplica] = uniformFunction(x)
fltModelReplica = zeros(1,length(x));

    for i = 1:1:length(x)
        if(x(i) < -0.5 || x(i) > 0.5)
            fltModelReplica(i) = 0;
        else
            fltModelReplica(i) = 1;
        end
    end
    
end

function [fltModelReplica] = triangleFunction(x)
fltModelReplica = zeros(1,length(x));

for indexTmp = 1:1:length(x)
    if(x < -1 || x > 1)
        fltModelReplica(indexTmp) = 0;
    else
        fltModelReplica(indexTmp) = 1 - abs(x(indexTmp));
    end
end
    
end

function [fltModelReplica] = gaussianFunction(x)
fltModelReplica = (1/sqrt(2*pi)).*exp(-(x.^2)/2);
end

function [fltModelReplica] = laplaceFunction(x)
fltModelReplica = (1/2).*exp(-abs(x));
end

function [fltModelReplica] = cauchyFunction(x)
fltModelReplica = (1/pi).*(1./(1+x.^2));
end

function [fltModelReplica] = sinhFunction(x)
fltModelReplica = zeros(1,length(x));

for i = 1:1:length(x)
    if(x(i) ~= 0)
        fltModelReplica(i) = (1/pi)*(sin(x(i))/(x(i)))^2;
    else
        fltModelReplica(i) = (1/pi);
    end
end

end