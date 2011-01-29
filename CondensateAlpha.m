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
            error('change this: alpha is now defined at nodes!!!')
            alpha = convvec( alpha, alpha0{ i1(i2), j1(i2) } );
        end
       
    % "random" FE case: the alpha must be corrected 
    case {'MonteCarloHomeFE'}
        if size(c2m,1)>1
            error(['this coupling has not been looked at for more' ...
                   ' than one coupling' ]);
        end
        switch model.random.law
            
            % uniform first-order marginal law
            case 'uniform'
                a = model.random.min;
                b = model.random.max;
                alpha = alpha0{ 1, find(c2m==n) };
                for i1=1:size(alpha,1)
                    alpha(i1,1)=fzero(@(x) ...
                 (exp((b-a)*x)-1)*((1-alpha(i1,1))+a*x)-(b-a)*x,alpha(i1,1));
                    alpha(i1,2)=fzero(@(x) ...
                 (exp((b-a)*x)-1)*((1-alpha(i1,2))+a*x)-(b-a)*x,alpha(i1,2));
                end
                alpha(alpha>1)=1;
                alpha(alpha<0)=0;
                
            otherwise
                error('unknown random field law')
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

% polynomial value to coefficients
function a = valpolyN( y, x )
y1 = y(:,1);
y2 = y(:,2);
x1 = x(:,1);
x2 = x(:,2);
a = [ (y1-y2)'./(x1-x2)';(x1.*y2-x2.*y1)'./(x1-x2)']';