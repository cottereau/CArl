function y = polyvalN( P, x )
% POLYVALN to generalize polyval function to vectors of polynomials, and
% possibly at several sets of points
%
% syntax: y = polyvalN( P, x )
% 
% P is an [N*Np matrix] of N polynomials of order Np-1
% x is an N*d matrix of coordinates
Np = size(P,2)-1;
y = P(:,end);
for i2=1:Np
    y = y + (P(:,Np+1-i2) .* x.^i2);
end
