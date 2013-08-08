function Test( type )
% launch a test or a series of tests
% 
%   syntax: Test(type)
%
% ONE-DIMENSIONAL TESTS
%   ____________________ _____________ _____________ ______________
%  |                    |             |             |              |
%  |       type         |  Model 1    |  Model 2    |    remark    |
%  |____________________|_____________|_____________|______________|
%  |                    |             |             |              |
%  | 'zoom1D'           | acoustic    | acoustic    | zoom case    |
%  | 'join1D'           | acoustic    | acoustic    | join case    |
%  | 'force1D'          | acoustic    | acoustic    | bulk load    |
%  | 'stoDet1D'         | acoustic    | acoustic(S) |              |
%  | 'stoSto1D'         | acoustic(S) | acoustic(S) |              |
%  |____________________|_____________|_____________|______________|
%
% TWO-DIMENSIONAL TESTS
%   ____________________ _____________ _____________ ______________
%  |                    |             |             |              |
%  |       type         |  Model 1    |  Model 2    |    remark    |
%  |____________________|_____________|_____________|______________|
%  |                    |             |             |              |
%  | 'zoom2D'           | acoustic    | acoustic    | zoom mesh    |
%  | 'join2D'           | acoustic    | acoustic    | join mesh    |
%  | 'force2D'          | acoustic    | acoustic    | bulk load    |
%  | 'comsol2D'         | comsol      | comsol      | acoustic     |
%  | 'stoDet2D    '     | acoustic(S) | acoustic    |              |
%  | 'stoSto2D    '     | acoustic(S) | acoustic(S) |              |
%  | 'FE2D'             | elastic     | elastic     |              |
%  | 'beam2D'           | beam        | beam        |              |
%  | 'beamFE2D'         | beam        | elastic     |              |
%  |____________________|_____________|_____________|______________|
%
% (S) = stochastic
%
% TECHNICAL TESTS
%   ____________________ ______________________________________________
%  |                    |                                              |
%  |       type         |                  remark                      |
%  |____________________|______________________________________________|
%  |                    |                                              |
%  | 'nonembedded2D_1'  | Intersection of non-embedded meshes          |
%  | 'indent2D'         | Intersection of a mesh with acute angles     |
%  |____________________|______________________________________________|
%
% SERIES
%  type='short' launches in series the series '1D', '2D' and 'mesh'
%  type='1D' launches in series the short 1D tests: 'zoom1D', 'join1D',
%               'force1D', 'MC1D', 'zoom1Dstosto'
%  type='2D' launches in series the short 2D tests: 'join2D', 'zoom2D',
%               'force2D', 'stoDet2D', 'stoSto2D'
%  type='mesh' launches in series the meshing tests: 'indent2D', 
%               'NonEmbedded2D_1'

% default
if nargin==0
    type = 'short';
end

% selection of type of test
switch lower(type)
    
    case {'zoom1d', 'force1d', 'join1d'}
        load(['Tests/' type '.mat']);
        [sol,out] = CArl( model, coupling, solver );
        plottest1D( out.model, sol );
        title(type)
        
    case {'join2d', 'force2d', 'zoom2d', 'stodet2d', 'stosto2d', ...
          'comsol2d'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plottest2D( out.model, sol );
        title(type)

    case {'fe2d'}
        load(['Tests/' type '.mat']);
        model{2}.FE2D.load = model{2}.FE2D.load*10;
        [ sol, out ] = CArl( model, coupling, solver );
        plottestvec2D( out.model, sol );

    case {'stodet1d','stosto1d'}
        load(['Tests/' type '.mat']);
        sol = CArl( model, coupling, solver );
        plotteststochastic( model, sol );
        title(type)
        
    case {'beam2d'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plottest1D( out.model, sol );
        title(type)

    case {'nonembedded2d_1', 'indent2d'}
        load(['Tests/' type '.mat']);
        [ ~, out ] = CArl( model, coupling, solver );
        plotmesh2D( out );
        title(type)

    case 'short'
        Test('1D');
        Test('2D');
        Test('mesh');
        
    case '1d'
        Test('zoom1D');
        Test('join1D');
        Test('force1D');
        Test('stoDet1D');
        Test('stoSto1D');

    case '2d'
        Test('join2D');
        Test('zoom2D');
        Test('force2D');
        Test('stoDet2D');
        Test('stoSto2D');
        Test('FE2D');

    case 'mesh'
        Test('indent2D');
        Test('NonEmbedded2D_1');
        
    otherwise
        error('unknown test case')

end

% FUNCTION PLOTTEST 1D
function plottest1D( model, sol )
x1 = model{1}.mesh.Points( model{1}.mesh.ConnectivityList )';
x2 = model{2}.mesh.Points( model{2}.mesh.ConnectivityList )';
u1 = sol{1}( model{1}.mesh.ConnectivityList )';
u2 = sol{2}( model{2}.mesh.ConnectivityList )';
figure; plot( x1(:), u1(:), 'bx-', x2(:), u2(:), 'ro--' )
s1 = diff(u1) ./ diff(x1); s1 = [s1; s1];
s2 = diff(u2) ./ diff(x2); s2 = [s2; s2];
figure; plot( x1, s1, 'bx-', x2(1,:), s2(1,:), 'r' );

% FUNCTION PLOTTEST 2D
function plottest2D( model, sol )
T1 = model{1}.mesh.ConnectivityList;
T2 = model{2}.mesh.ConnectivityList;
X1 = model{1}.mesh.Points;
X2 = model{2}.mesh.Points;
N1 = size(X1,1);
N2 = size(X2,1);
figure; trimesh( T1, X1(:,1), X1(:,2), sol{1}(1:N1) ); 
hold on; trisurf( T2, X2(:,1), X2(:,2), sol{2}(1:N2) ); 
colorbar; view(3)

% FUNCTION PLOTTESTSTOCHASTIC
function plotteststochastic( model, sol )
pc = 0.9;
x1 = model{1}.HomeFE.mesh.X;
N1 = length(x1);
mu1 = mean(sol{1}(1:N1,:),2);
stdu1 = std(sol{1}(1:N1,:),0,2);
probu1=[mu1+stdu1/sqrt(1-pc);(mu1(end:-1:1)-stdu1(end:-1:1)/sqrt(1-pc))];
xxx1 = [ x1;x1(end:-1:1)];
figure
fill( xxx1, probu1, 'y' )
x2 = model{2}.HomeFE.mesh.X;
N2 = length(x2);
mu2 = mean(sol{2}(1:N2,:),2);
stdu2 = std(sol{2}(1:N2,:),0,2);
probu2=[mu2+stdu2/sqrt(1-pc);mu2(end:-1:1)-stdu2(end:-1:1)/sqrt(1-pc)];
xxx2 = [ x2;x2(end:-1:1)];
hold on
fill( xxx2, probu2, 'y' )
plot( x1, mu1, 'bx-', x2, mu2, 'r--' );
dx2 = diff(x2);
dx1 = diff(x1);
s1 = diff(sol{1}(1:N1,:),[],1) ./ repmat(dx1,[1 size(sol{1},2)]);
ms1 = mean(s1,2); ms1 = [ms1 ms1]'; ms1 = ms1(:);
ss1 = std(s1,[],2); ss1 = [ss1 ss1]'; ss1 = ss1(:);
pms1 = [ ms1+ss1/sqrt(1-pc); ms1(end:-1:1)-ss1(end:-1:1)/sqrt(1-pc) ];
x1 = model{1}.HomeFE.mesh.X(model{1}.HomeFE.mesh.T)';
x2 = model{2}.HomeFE.mesh.X(model{2}.HomeFE.mesh.T)';
xx1 = x1(:); xx1 = [ xx1; xx1(end:-1:1) ];
figure; fill( xx1, pms1, 'y' )
s2 = diff(sol{2}(1:N2,:),[],1) ./ repmat(dx2,[1 size(sol{2},2)]);
ms2 = mean(s2,2); ms2 = [ms2 ms2]'; ms2 = ms2(:);
ss2 = std(s2,[],2); ss2 = [ss2 ss2]'; ss2 = ss2(:);
pms2 = [ ms2+ss2/sqrt(1-pc); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-pc) ];
xx2 = x2(:); xx2 = [ xx2; xx2(end:-1:1) ];
hold on;
fill( xx2, pms2, 'y' )
plot( x1(:), ms1, 'bx-', x2(:), ms2, 'ro--' );

% FUNCTION PLOTTESTVEC2D
function plottestvec2D( model, sol )
T1 = model{1}.FE2D.T;
X1 = model{1}.FE2D.X;
N1 = numel(X1)/2;
u1 = reshape( sol{1}(1:2*N1), 2, N1 )';
dX1 = X1+u1;
T2 = model{2}.FE2D.T;
X2 = model{2}.FE2D.X;
N2 = numel(X2)/2;
u2 = reshape( sol{2}(1:2*N2), 2, N2 )';
dX2 = X2+u2;
figure; trimesh( T1, X1(:,1), X1(:,2), 'color', 'k', 'LineStyle', ':' );
hold on; trimesh( T1, dX1(:,1), dX1(:,2), 'color', 'b' );
hold on; trimesh( T2, X2(:,1), X2(:,2), 'color', 'k', 'LineStyle', ':' );
hold on; trimesh( T2, dX2(:,1), dX2(:,2), 'color', 'r' );
