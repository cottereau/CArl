function alpha = CondensateAlpha( n, c2m, alpha0 )
% CONDENSATEALPHA to take the alpha functions computed independently for
% each coupling, and construct the corresponding alpha function for model
% n, taking into account all the couplings with different models
%
% syntax: alpha = CondensateAlpha( n, model, c2m, alpha )
%
%  n: index of the model being considered
%  c2m: [Nc*2 matrix] indicating the two models at play for each coupling
%  alpha1,alpha2: [Nc*1 cell] describing the alpha function for the first 
%        and second models of a coupling, respectively

% R. Cottereau 01/2011

[i1,j1] = find( c2m==n );
alpha = alpha0{ i1(1), j1(1) };
for i2 = 2:length(i1)
    error('Several coupling not implemented yet')
end
