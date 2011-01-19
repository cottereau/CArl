function [ X, T, Xn ] = ReduceMesh( X, T )
% REDUCEMESH to get rid of unused nodes and renumber the connectivity
% matrix
%
% syntax: [ X, T, Xn ] = ReduceMesh( X, T )
%
%  T,X: connectivity matrices and nodal coordinates [Ne*3 matrix] and 
%       [Nn*2 matrix]
%  Xn: indices into the input connectivity matrix
Xn = unique( T(:) );
Nn = size(Xn,1);
for i1 = 1:Nn
    T( T(:)==Xn(i1) ) = i1;
end
X = X(Xn,:);
