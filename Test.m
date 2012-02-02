function Test(type)
% launch a test or a series of tests
% The existing tests are the following
%  'zoom1D' for a zoom case in 1D (deterministic continuum-continuum)
%  'join1D' for a junction case in 1D (deterministic continuum-continuum)
%  'MC1D' for a joint case in 1D (stochastic continuum-continuum)
%  'join2D' for  a junction case in 2D (deterministic continuum-continuum)
%
% the possible series are
%  'short' for only the short tests
%  'long' for the longer tests

% R. Cottereau 05/2010
% TO DO
% - a test with three models and two couplings

% default
if nargin==0
    type = 'short';
end

% selection of type of test
switch type
    
    case {'zoom1D','join1D','join1D_fine'}
        load(['Tests/' type '.mat']);
        sol = CArl( model, coupling, solver );
        plottest1D( model, sol );
        
    case {'join2D','join2D_fine'}
        load(['Tests/' type '.mat']);
        sol = CArl( model, coupling, solver );
        plottest2D( model, sol );

    case {'MC1D','MC1D_fine'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plotteststochastic( model, sol, out );
        
    case 'short'
        Test('zoom1D');
        Test('join1D');
        Test('MC1D');
        Test('join2D');
        
    case 'full'
        Test('short');
        Test('join1D_fine');
        Test('MC1D_fine');

end

% FUNCTION PLOTTEST 1D
function plottest1D( model, sol )
figure; plot(model{1}.mesh.X,sol{1}, 'bx-', model{2}.mesh.X,sol{2},'ro--')

% FUNCTION PLOTTEST 2D
function plottest2D( model, sol )
figure;
trimesh( model{1}.mesh.T, model{1}.mesh.X(:,1), model{1}.mesh.X(:,2), sol{1});
hold on ;
trisurf( model{2}.mesh.T, model{2}.mesh.X(:,1), model{2}.mesh.X(:,2), sol{2});


% FUNCTION PLOTTESTSTOCHASTIC
function plotteststochastic( model, sol, out )
X1 = model{1}.mesh.X;
s1 = diff( sol{1} ) / mean(diff(X1)); s1 = [s1 s1]'; s1 = s1(:);
x1 = [X1(1:(end-1)) X1(2:end)]'; x1 = x1(:);
X2 = model{2}.mesh.X;
x2 = [X2(1:(end-1)) X2(2:end)]'; x2 = x2(:);
s2 = zeros( length(x2), size(out.u.MC.u,2));
s2(1:2:end,:) = diff( out.u.MC.u(1:length(X2),:), 1, 1 ) / mean(diff(X2));
s2(2:2:end,:) = s2(1:2:end,:);
ms2 = mean(s2,2);
ss2 = std(s2,[],2);
pms2 = [ ms2+ss2/sqrt(1-0.9); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-0.9) ];
figure; fill( [x2;x2(end:-1:1)], pms2, 'y' )
hold on; plot( x1, s1, 'bx-', x2, ms2, 'ro--' );

