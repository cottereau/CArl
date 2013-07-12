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
%  | 'MC1Dstosto'       | acoustic   | acoustic   |  S  |              |
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
%  type='short' launches in series the following tests: 'zoom1D', 'join1D',
%               'force1D', 'MC1D', 'MC1D_BC_u10', 'MC1D_stosto', 
%               'zoom1Dstosto', 'join2D', 'zoom2D', 'zoom2Dindent', 
%               'NonEmbedded2D_1', 'force2D', 'zoom2dstosto'

% creation: R. Cottereau 05/2010
% TO DO
% - a test with three models and two couplings

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
        [ sol, out ] = CArl( model, coupling, solver );
        plotteststochastic2( model, sol, out, ref );
        
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
function plotteststochastic2( model, ~, out, ref )
pc = 0.9;
x1 = model{1}.mesh.X( model{1}.mesh.T )';
x2 = model{2}.mesh.X( model{2}.mesh.T )';
dx2 = diff(x2);
dx1 = diff(x1);
u1=mean(out.model{1}.auMC,3)';uu1=u1(:);
stdu1=std(out.model{1}.auMC,0,3)';stdu1=stdu1(:);
probu1=[uu1+stdu1/sqrt(1-pc);(uu1(end:-1:1)-stdu1(end:-1:1)/sqrt(1-pc))];
xxx1 = [ (x1(:));(x1(end:-1:1))'];
figure;
if (strcmp(out.coupling{1}.mediator.type,'mesomicro'))
    fill( xxx1, probu1, 'y' )
end
hold on
%u1b = (sol{1})';
%u2b=(sol{2})';
u2=mean(out.model{2}.auMC,3)';uu2=u2(:);%mean u2
stdu2=std(out.model{2}.auMC,0,3)';stdu2=stdu2(:);%std u2
probu2=[uu2+stdu2/sqrt(1-pc);uu2(end:-1:1)-stdu2(end:-1:1)/sqrt(1-pc)];
xxx2 = [ x2(:);x2(end:-1:1)'];
fill( xxx2(:), probu2(:), 'y' )
plot(x1,u1,'bx-',x2,u2,'r--');


%s1b = diff(u1) ./ diff(x1);% s1b = [s1b; s1b];
%s2b = diff(u2) ./ diff(x2);% s2b = [s2b; s2b];

s1 = squeeze(diff(out.model{1}.auMC,[],2)) ...
           ./ repmat(dx1',[1 size(out.model{1}.auMC,3)]);
ms1 = mean(s1,2); ms1 = [ms1 ms1]'; ms1 = ms1(:);
ss1 = std(s1,[],2); ss1 = [ss1 ss1]'; ss1 = ss1(:);
pms1 = [ ms1+ss1/sqrt(1-pc); ms1(end:-1:1)-ss1(end:-1:1)/sqrt(1-pc) ];
xx1 = x1(:); xx1 = [ xx1; xx1(end:-1:1) ];
figure;
if (strcmp(out.coupling{1}.mediator.type,'mesomicro'))
fill( xx1, pms1, 'y' )
end
hold on;
s2 = squeeze(diff(out.model{2}.auMC,[],2)) ...
           ./ repmat(dx2',[1 size(out.model{2}.auMC,3)]);
ms2 = mean(s2,2); ms2 = [ms2 ms2]'; ms2 = ms2(:);
ss2 = std(s2,[],2); ss2 = [ss2 ss2]'; ss2 = ss2(:);
pms2 = [ ms2+ss2/sqrt(1-pc); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-pc) ];
xx2 = x2(:); xx2 = [ xx2; xx2(end:-1:1) ];
fill( xx2, pms2, 'y' )
plot( x1(:), ms1, 'bx-', x2(:), ms2, 'ro--' );
if ~isempty(ref)
    xref = [ref.X(1:(end-1)) ref.X(2:end)]'; xref = xref(:);
    sref = diff(ref.u) ./ diff(ref.X); sref = [sref; sref];
    hold on; plot(xref,sref,'k--')
end

