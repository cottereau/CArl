function alpha = CondensateAlpha( n, model, c2m, alpha0 )
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

% switch on the different possible codes
switch model.code
    
    % "classical" FE case: multiply the alpha for each coupling
    case {'HomeFE', 'Comsol'}
        [i1,j1] = find( c2m==n );
        alpha = alpha0{ i1(1), j1(1) };
        for i2 = 2:length(i1)
            alpha = convvec( alpha, alpha0{ i1(i2), j1(i2) } );
        end
       
    otherwise
        error( 'not implemented yet' )
        
end

% convolution of matrices
function a = convvec( a, b )
Na = size(a);
Nb = size(b);
a = fft( [a zeros(Na(1),Nb(2)-1)], [], 2 );
b = fft( [b zeros(Nb(1),Na(2)-1)], [], 2 );
a = ifft( a .* b, [], 2 );
