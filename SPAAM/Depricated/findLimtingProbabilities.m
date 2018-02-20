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
function[est] = findLimtingProbabilities(A)

X = A;
[n,m] = size(X);

for(i = 1:1:n)
    
    X(i,i) = A(i,i) - 1;
    
end

X(n + 1,:) = ones(1,m);
Y = [zeros(1,n), 1];

est = (X'*X)^(-1)*X'*Y';

end