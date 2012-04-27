function alpha = CondensateAlpha( n, model, Cpl )
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

% constant
Nc = length(Cpl);

% "classical" FE case: multiply the alpha for each coupling
[i1,j1] = find( Cpl{1}.models==n );
if j1==1
    alpha = Cpl{i1}.alpha1;
else
    alpha = Cpl{i1}.alpha2;
end
for i2 = 2:Nc
    error('Several coupling not implemented yet')
end

% switch on the different possible codes
switch model.code
    
    % "random" FE case: the alpha must be corrected 
    case {'MonteCarloHomeFE'}
        switch model.random.law
            
            % uniform first-order marginal law
            case 'uniform'
                % bounds of the uniform distribution
                a = model.random.min;
                b = model.random.max;
                for i1=1:size(alpha,1)
                    alpha(i1,1)=fzero(@(x) ...
                 (exp((b-a)*x)-1)*((1-alpha(i1,1))+a*x)-(b-a)*x,alpha(i1,1));
                    alpha(i1,2)=fzero(@(x) ...
                 (exp((b-a)*x)-1)*((1-alpha(i1,2))+a*x)-(b-a)*x,alpha(i1,2));
                end
                alpha(alpha>1)=1;
                alpha(alpha<0)=0;
                
            case 'lognormal'
                % variance of the underlying gaussian distribution
                Nmc = 1e4;
                Nx = 30;
                sig2 = model.random.sig2;
                k = exp(randn(Nmc,1)*sqrt(sig2)+sig2/2);
                a = linspace( 0, 1, Nx );
                b = 1-a;
                for i1 = 2:Nx-1
                    a(i1) = fzero(@(a)mean(1./(b(i1)+a*k))-1,a(i1));
                end
                alpha(:,1) = interp1( 1-b, a, alpha(:,1) );
                alpha(:,2) = interp1( 1-b, a, alpha(:,2) );
                
            otherwise
                error('unknown random field law')
        end
        
end
