function d = X2X( X1, X2 )
% X2X to compute the distances from a set of points to another set of
% points
%
% syntax: d = X2X( X1, X2 )
%
%  X1,X2: lists of points [N1*2 matrix]Ê[N2*2 matrix]
%
%  d: distance matrix [N1*N2 matrix]
n1 = size(X1,1);
n2 = size(X2,1);
d=sqrt( ( repmat(X2(:,1)',n1,1) - repmat(X1(:,1),1,n2) ).^2 + ...
        ( repmat(X2(:,2)',n1,1) - repmat(X1(:,2),1,n2) ).^2);
