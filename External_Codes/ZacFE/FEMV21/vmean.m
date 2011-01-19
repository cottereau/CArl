function y=vmean(x)
%VMEAN Mean value of vector elements.
%   Y = VMEAN(X) where length(Y) = length(X)-1
%   If X = [1 2 6 0] then Y = [1.5 4 3]
%
%   See also MEDIAN, STD, MIN, MAX, COV, DIFF, SUM.

% Copyright (c) 2002-03-24, B. Rasmus Anthin.

y=(x(1:end-1)+x(2:end))/2;