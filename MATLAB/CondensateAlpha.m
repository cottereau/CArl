function alpha = CondensateAlpha( n, Cpl )
% CONDENSATE to take the alpha functions computed independently for
% each coupling, and construct the corresponding alpha function for model
% n, taking into account all the couplings with different modesl
%
% syntax: alpha = CondensateAlpha( n, model, c2m, alpha )
%
%  n: index of the model being considered
%  model: structured array describing the model, including field 'code'
%  c2m: [Nc*2 matrix] indicating the two models at play for each coupling
%  alpha1,alpha2: [Nc*1 cell] describing the alpha function for the first 
%        and second models of a coupling, respectively

% R. Cottereau 01/2011

% "classical" FE case: multiply the alpha for each coupling
[i1,j1] = find( Cpl{1}.models==n );
if length(i1)==1
    alpha = Cpl{i1}.alpha{j1};
else
    error('Several coupling not implemented yet')
end
