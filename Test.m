function Test(type)
% launch a test or a series of tests
% The existing tests are the following
%  'zoom1D' for a zoom case in 1D (deterministic continuum-continuum)
%  'join1D' for a junction case in 1D (deterministic continuum-continuum)
%  'MC1D' for a joint case in 1D (stochastic continuum-continuum)
%  'join2D' for  a junction case in 2D (deterministic continuum-continuum)
%
% the possible series are
%  'short' for only the short tests: zoom1D, join1D, MC1D, join2D
%  'long' for the longer tests

% R. Cottereau 05/2010
% TO DO
% - a test with three models and two couplings

% default
if nargin==0
    type = 'short';
end

% selection of type of test
switch lower(type)
    
    case {'zoom1d', 'join1d', 'join1d_fine'}
        load(['Tests/' type '.mat']);
        sol = CArl( model, coupling, solver );
        plottest1D( model, sol );
        
    case {'join2d', 'zoom2d', 'join2d_fine', 'comsol2d', 'nonembedded2d_1'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plottest2D( out.model, sol );

    case {'mc1d', 'mc1d_fine'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plotteststochastic( model, sol, out );
        
    case 'short'
        Test('zoom1D');
        Test('join1D');
        Test('MC1D');
        Test('join2D');
        Test('zoom2D')
        
    case 'full'
        Test('short');
        Test('join1D_fine');
        Test('MC1D_fine');
        
    case 'comsol'
        Test('comsol2d');
        
    otherwise
        error('unknown test case')

end

% FUNCTION PLOTTEST 1D
function plottest1D( model, sol )
x1 = model{1}.mesh.X( model{1}.mesh.T )';
x2 = model{2}.mesh.X( model{2}.mesh.T )';
u1 = sol{1}';
u2 = sol{2}';
figure; plot( x1, u1, 'bx-', x2, u2, 'ro--' )
s1 = diff(u1) ./ diff(x1); s1 = [s1; s1];
s2 = diff(u2) ./ diff(x2); s2 = [s2; s2];
figure; plot( x1, s1, 'bx-', x2, s2, 'ro--' );

% FUNCTION PLOTTEST 2D
function plottest2D( model, sol )
T1 = model{1}.mesh.Triangulation;
T2 = model{2}.mesh.Triangulation;
x1 = model{1}.mesh.X(:,1); x1 = x1(T1);
y1 = model{1}.mesh.X(:,2); y1 = y1(T1);
x2 = model{2}.mesh.X(:,1); x2 = x2(T2);
y2 = model{2}.mesh.X(:,2); y2 = y2(T2);
figure; patch(x1',y1',sol{1}',sol{1}')
hold on; patch(x2',y2',sol{2}',sol{2}')
view(3)

% FUNCTION PLOTTESTSTOCHASTIC
function plotteststochastic( model, sol, out )
pc = 0.9;
x1 = model{1}.mesh.X( model{1}.mesh.T )';
x2 = model{2}.mesh.X( model{2}.mesh.T )';
dx2 = diff(x2);
u1 = sol{1}';
s1 = diff(u1) ./ diff(x1); s1 = [s1; s1];
s2 = squeeze(diff(out.model{2}.auMC,[],2)) ...
           ./ repmat(dx2',[1 size(out.model{2}.auMC,3)]);
ms2 = mean(s2,2); ms2 = [ms2 ms2]'; ms2 = ms2(:);
ss2 = std(s2,[],2); ss2 = [ss2 ss2]'; ss2 = ss2(:);
pms2 = [ ms2+ss2/sqrt(1-pc); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-pc) ];
xx2 = x2(:); xx2 = [ xx2; xx2(end:-1:1) ];
figure; fill( xx2, pms2, 'y' )
hold on; 
plot( x1, s1, 'bx-', x2(:), ms2, 'ro--' );

