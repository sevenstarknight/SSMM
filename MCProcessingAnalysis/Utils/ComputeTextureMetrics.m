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
function [fltSet] = ComputeTextureMetrics(col)

% Should be a Square matrix coming in
[n_g] = length(col);
fltSet = zeros(1,10);

%% Compute Marginal Probability Vectors
p_x = sum(col, 2);
p_y = sum(col, 1);

structMarginal = struct(...
    'p_x', p_x, ...
    'p_y', p_y, ...
    'mu_x', mean(p_x), ...
    'mu_y', mean(p_y), ...
    'std_x', sqrt(var(p_x)), ...
    'std_y', sqrt(var(p_y)));


%% Compute Metrics
[asm] = ComputeAngularSecondMoment(col, n_g);

[contrast] = ComputeContrast(col, n_g);

% [corr] = ComputeCorrolation(col, n_g, structMarginal);

[ssq] = ComputeSSQ(col, n_g);

 [idm] = ComputeIDM(col, n_g);

 [sum_avg] = ComputeSumAvg(col, n_g);

[sum_ent] = ComputeSumEntropy(col, n_g);

[sum_var] = ComputeSumVariance(col, n_g, sum_ent);

[ent] = ComputeEntropy(col, n_g);

[diff_var] = ComputeDiffVar(col,n_g);

[diff_entropy] = ComputeDiffEnt(col, n_g);

% [f12, f13] = ComputeInformationMetric(col, n_g, structMarginal);

%% Initialize Output Array

% 13 metrics
% structMetrics = struct('asm', asm, ...
%     'contrast', contrast, ...
%     'corr', corr, ...
%     'ssq', ssq, ...
%     'idm', idm, ...
%     'sum_avg', sum_avg, ...
%     'sum_ent', sum_ent, ...
%     'sum_var', sum_var,...
%     'ent', ent, ...
%     'diff_var', diff_var, ...
%     'diff_entropy', diff_entropy, ...
%     'f12', f12, ...
%     'f13', f13);
   
    fltSet(1) =  asm;
    fltSet(2) =  contrast;
%     fltSet(3) =  corr;
    fltSet(3) =  ssq;
     fltSet(4) =  idm;
     fltSet(5) =  sum_avg;
    fltSet(6) =  diff_var;
    fltSet(7) =  sum_var;
    fltSet(8) =  ent;
    fltSet(9) =  diff_entropy;
%     fltSet(11) =  f12;
%     fltSet(12) =  f13;
    fltSet(10) = sum_ent;
    
% structMetrics = struct('asm', asm, ...
%     'contrast', contrast, ...
%     'corr', corr, ...
%     'ssq', ssq, ...
%     'idm', idm, ...
%     'sum_avg', sum_avg, ...
%     'diff_var', diff_var);


end

function [asm] = ComputeAngularSecondMoment(col, n_g)
%% Compute Angular Second Moment
asm = 0;
for i = 1:1:n_g
    for j = 1:1:n_g
        asm = asm + col(i,j)^2;
    end
end
end

function [contrast] = ComputeContrast(col, n_g)
%% Compute Contrast
contrast = 0;
for n = 0:1:(n_g - 1)
    contrast = contrast + (n^2)*P_X_Minus_Y(col, n_g, n); 
end
end

function [corr] = ComputeCorrolation(col, n_g, structMarginal)
corr = 0;
for i = 1:1:n_g
    for j = 1:1:n_g
        corr = corr + ((i*j)*col(i,j) - structMarginal.mu_x*structMarginal.mu_y);
    end
end
corr = corr/(structMarginal.std_x*structMarginal.std_y);
end


function [ssq] = ComputeSSQ(col, n_g)
%% Compute Sum Of Squares (Variance)
ssq = 0;
for i = 1:1:n_g
    for j = 1:1:n_g
        ssq = ssq + ((i-j)^2)*col(i,j);
    end
end
end


function [idm] = ComputeIDM(col, n_g)
%% Inverse Difference Moment
idm = 0;
for i = 1:1:n_g
    for j = 1:1:n_g
        idm = idm + (1/(1 + (i-j)^2))*col(i,j);
    end
end
end


function [sum_avg] = ComputeSumAvg(col, n_g)
%% Sum Average
sum_avg = 0;
for n = 2:1:2*n_g
    sum_avg = sum_avg + n*P_X_Plus_Y(col, n_g, n);
end
end

function [sum_ent] = ComputeSumEntropy(col, n_g)
%% Sum Entropy
sum_ent = 0;
for n = 2:1:2*n_g
    tmpP_X_Plus_Y = P_X_Plus_Y(col, n_g, n);
    sum_ent = sum_ent + tmpP_X_Plus_Y*log(tmpP_X_Plus_Y +eps(tmpP_X_Plus_Y));
end
sum_ent = -sum_ent;
end

function [sum_var] = ComputeSumVariance(col, n_g, sum_ent)
%% Sum Variance
sum_var = 0;
for n = 2:1:2*n_g
    x = P_X_Plus_Y(col, n_g, n);
    sum_var = sum_var + (n - sum_ent)*log(x + eps(x));
end
end


function [ent] = ComputeEntropy(col, n_g)
%% Entropy
ent = 0;
for i = 1:1:n_g
    for j = 1:1:n_g
        ent = ent + col(i,j)*log(col(i,j) + eps(col(i,j)));
    end
end
ent = -ent;
end


function [diff_var] = ComputeDiffVar(col,n_g)
%% Difference Variance
diff_var_array = zeros(1,n_g-1);

for n = 0:1:(n_g - 1)
    diff_var_array(n + 1) = P_X_Minus_Y(col, n_g, n); 
end

diff_var = var(diff_var_array);
end


function [diff_entropy] = ComputeDiffEnt(col, n_g)
%% Difference Entropy
diff_entropy = 0;
for n = 0:1:(n_g - 1)
    tmpP_X_Minus_Y = P_X_Minus_Y(col, n_g, n);
    diff_entropy = diff_entropy + tmpP_X_Minus_Y*log(tmpP_X_Minus_Y + eps(tmpP_X_Minus_Y)); 
end
diff_entropy = -diff_entropy;
end




function [f12, f13] = ComputeInformationMetric(col, n_g, structMarginal)
%% Information Measures of Correlation
HXY = 0;
HXY1 = 0;
HXY2 = 0;

for i = 1:1:n_g
    for j = 1:1:n_g
        
        x = structMarginal.p_x(i)*structMarginal.p_y(j);
        
        HXY = HXY + col(i,j)*log(col(i,j)+ eps(x));
        HXY1 = HXY1 + col(i,j)*log(x + eps(x));
        HXY2 = HXY2 + structMarginal.p_x(i)*structMarginal.p_y(j)*log(x + eps(x));
    end
end

HX = 0;
for i = 1:1:n_g
    x = structMarginal.p_x(i);
    HX = HX + structMarginal.p_x(i)*log(x + eps(x));
end

HY = 0;
for j = 1:1:n_g
    x = structMarginal.p_y(j);
    HY = HY + structMarginal.p_y(j)*log(x + eps(x));
end

HXY = -HXY;
HXY1 = -HXY1;
HXY2 = -HXY2;

f12 = (HXY-HXY1)/max(HX,HY);
f13 = sqrt(1 - exp(-2*(HXY2 - HXY)));

end



function [pxy] = P_X_Plus_Y(col, n_g, k)
pxy = 0;

for i = 1:1:n_g
    for j = 1:1:n_g
        if((i+j) == k )
            pxy = pxy + col(i,j);
        end
    end
end

end

function [pxy] = P_X_Minus_Y(col, n_g, k)
pxy = 0;

for i = 1:1:n_g
    for j = 1:1:n_g
        if(abs(i-j) == k )
            pxy = pxy + col(i,j);
        end
    end
end

end