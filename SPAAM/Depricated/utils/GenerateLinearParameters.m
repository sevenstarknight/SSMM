function [A, B] = GenerateLinearParameters(x, y)

N = length(x);
ssqX = sum(x.^2);
sxy = sum(x.*y);
sx = sum(x);
sy = sum(y);

delta = N*ssqX - (sx^2);

% y = A + Bx
A = (ssqX*sy - sx*sxy)/delta;
B = (N*sxy - sx*sy)/delta;

end