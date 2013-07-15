function Test( type )
% launch a test or a series of tests
% 
%   syntax: Test(type)
%
% ONE-DIMENSIONAL TESTS
%   ____________________ ____________ ____________ _____ ______________
%  |                    |            |            |     |              |
%  |       type         |  Model 1   |  Model 2   | S/D |    remark    |
%  |____________________|____________|____________|_____|______________|
%  |                    |            |            |     |              |
%  | 'zoom1D'           | acoustic   | acoustic   |  D  | zoom case    |
%  | 'join1D'           | acoustic   | acoustic   |  D  | join case    |
%  | 'force1D'          | acoustic   | acoustic   |  D  | bulk load    |
%  | 'MC1D'             | acoustic   | acoustic   | D/S |              |
%  | 'zoom1Dstosto'     | acoustic   | acoustic   |  S  | zoom case    |
%  |____________________|____________|____________|_____|______________|
%
% TWO-DIMENSIONAL TESTS
%   ____________________ ____________ ____________ _____ ______________
%  |                    |            |            |     |              |
%  |       type         |  Model 1   |  Model 2   | S/D |    remark    |
%  |____________________|____________|____________|_____|______________|
%  |                    |            |            |     |              |
%  | 'zoom2D'           | acoustic   | acoustic   |  D  | zoom case    |
%  | 'join2D'           | acoustic   | acoustic   |  D  | join case    |
%  | 'force2D'          | acoustic   | acoustic   |  D  | bulk load    |
%  | 'comsol2D'         | comsol     | comsol     |  D  | acoustic     |
%  | 'zoom2Dindent'     | acoustic   | acoustic   |  S  | ??           |
%  | 'beam2D'           | beam       | elastic    | D/S |              |
%  |____________________|____________|____________|_____|______________|
%
% D = deterministic
% S = stochastic
%
% TECHNICAL TESTS
%   ____________________ ______________________________________________
%  |                    |                                              |
%  |       type         |                  remark                      |
%  |____________________|______________________________________________|
%  |                    |                                              |
%  | 'nonembedded2D_1'  | Intersection of non-embedded meshes          |
%  |____________________|______________________________________________|
%
% SERIES
%  type='short' launches in series the series '1D' and the tests :
%               'join2D', 'zoom2D', 'zoom2Dindent', 
%               'NonEmbedded2D_1', 'force2D', 'zoom2dstosto'
%  type='1D' launches in series the short 1D tests: 'zoom1D', 'join1D',
%               'force1D', 'MC1D', 'zoom1Dstosto'

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
        
    case {'join2d', 'force2d', 'zoom2d', 'zoom2dstosto', 'comsol2d', ...
          'nonembedded2d_1', 'zoom2dindent'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plottest2D( out.model, sol );

    case {'mc1d','zoom1dstosto'}
        load(['Tests/' type '.mat']);
        sol = CArl( model, coupling, solver );
        plotteststochastic( model, sol );
        
    case {'beam2d'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plottest2D( out.model, sol );

    case 'short'
        Test('zoom1D');
        Test('join1D');
        Test('force1D');
        Test('MC1D');
        Test('zoom1Dstosto');
        Test('join2D');
        Test('zoom2D');
        Test('zoom2Dindent');
        Test('NonEmbedded2D_1');
        Test('force2D');
        Test('zoom2dstosto');
        
    case '1d'
        Test('zoom1D');
        Test('join1D');
        Test('force1D');
        Test('MC1D');
        Test('zoom1Dstosto');

    otherwise
        error('unknown test case')

end

% FUNCTION PLOTTEST 1D
function plottest1D( model, sol, ref )
if nargin<3
    ref = [];
end
x1 = model{1}.mesh.X( model{1}.mesh.T )';
x2 = model{2}.mesh.X( model{2}.mesh.T )';
u1 = sol{1}( model{1}.mesh.T )';
u2 = sol{2}( model{2}.mesh.T )';
figure; plot( x1(:), u1(:), 'bx-', x2(:), u2(:), 'ro--' )
if ~isempty(ref)
    hold on; plot(ref.X,ref.u,'k--')
end
s1 = diff(u1) ./ diff(x1); s1 = [s1; s1];
s2 = diff(u2) ./ diff(x2); s2 = [s2; s2];
figure; plot( x1, s1, 'bx-', x2(1,:), s2(1,:), 'r' );
if ~isempty(ref)
    xref = zeros( 2*(length(ref.X)-1), 1);
    xref(1:2:end,:) = ref.X(1:end-1) ;
    xref(2:2:end,:) = ref.X(2:end) ;
    sref = zeros( 2*(length(ref.X)-1), 1);
    sref(1:2:end,:) = diff( ref.u, 1 ) / mean(diff(ref.X));
    sref(2:2:end,:) = sref(1:2:end,:);
    hold on; plot(xref,sref,'k--')
end

% FUNCTION PLOTTEST 2D
function plottest2D( model, sol )
T1 = model{1}.mesh.T3;
T2 = model{2}.mesh.T3;
X1 = model{1}.mesh.X3;
X2 = model{2}.mesh.X3;
figure; trimesh( T1, X1(:,1), X1(:,2), sol{1} ); 
hold on; trisurf( T2, X2(:,1), X2(:,2), sol{2} ); 
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

