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
%  | 'beam2D'           | beam        | beam        |              |
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
%  type='1D' launches in series the short 1D tests: 'zoom1D', 'join1D',
%               'force1D', 'MC1D', 'zoom1Dstosto'
%  type='2D' launches in series the short 2D tests: 'join2D', 'zoom2D',
%               'force2D', 'stoDet2D', 'stoSto2D'
%  type='mesh' launches in series the meshing tests: 'indent2D', 
%               'NonEmbedded2D_1'
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% default
if nargin==0
    type = 'short';
end
type = lower(type);

% load case parameters
if ~( strcmp(type,'1d') || strcmpi(type,'2d') || strcmpi(type,'mesh') )
    
    % get model and solve
    load( ['Tests/' type '.mat'] );
    sol = CArl( model, coupling, solver );
    
    % plot results
    figure; plotTestCArl( model{1}, sol{1}, 'k', '--x' )
    hold on; plotTestCArl( model{2}, sol{2}, 'r', '--o' );
    subplot(2,1,1); title( [type ' - displacement'] ); box on; grid on
    subplot(2,1,2); title( [type ' - gradient'] ); box on; grid on
    
    return
end

% selection of type of test
switch type
    
    case '1d'
        Test('zoom1D');
        Test('join1D');
        Test('force1D');
        Test('stoDet1D');
        Test('stoSto1D');
        Test('Beam2D');
        
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

% FUNCTION PLOTTESTCARL
function plotTestCArl( model, u, coul, style )

switch model.code
    
    % ACOUSTIC
    case 'HomeFE'
        X = model.HomeFE.mesh.X;
        T = model.HomeFE.mesh.T;
        d = size(X,2);
        % 1D
        if d==1
            du = diff(u) ./ diff(X); du = [du du]'; 
            dX = [X(1:end-1,1) X(2:end,1)]';
            subplot(2,1,1); hold on; plot( X, u, [coul style] );
            subplot(2,1,2); hold on; plot( dX(:), du(:), [coul style] );
        elseif d==2
            subplot(2,1,1); hold on; trimesh( T, X(:,1), X(:,2), u );
            view(3)
        end

    % ACOUSTIC STOCHASTIC
    case 'MonteCarloHomeFE'
        X = model.HomeFE.mesh.X;
        T = model.HomeFE.mesh.T;
        d = size(X,2);
        % 1D
        if d==1
            Nmc = size(u,2);
            mu = mean(u,2);
            su = std(u,[],2)/sqrt(1+0.9);
            du = diff(u,1,1) ./ repmat(diff(X),[1 Nmc]); 
            mdu = mean(du,2); mdu = [mdu mdu]'; mdu = mdu(:);
            sdu = std(du,[],2); sdu = [sdu sdu]'; sdu = sdu(:);
            dX = [X(1:end-1,1) X(2:end,1)]'; dX2 = [dX(:);dX(end:-1:1)'];
            subplot(2,1,1); hold on; 
            fill( [X;X(end:-1:1)], [mu+su;mu(end:-1:1)-su(end:-1:1)], 'y' );
            subplot(2,1,2); hold on; 
            fill( dX2, [mdu+sdu;mdu(end:-1:1,1)-sdu(end:-1:1,1)], 'y' );
            subplot(2,1,1); hold on; plot( X, mu, [coul style] );
            subplot(2,1,2); hold on; plot( dX(:), mdu(:), [coul style] );
        elseif d==2
            mu = mean(u,2);
            subplot(2,1,1); hold on; trimesh( T, X(:,1), X(:,2), mu );
            view(3)            
        end
        
    % ELASTIC - FE2D
    case 'FE2D'
        X = model.FE2D.X;
        T = model.FE2D.T;
        dx = X(:,1)+u(:,1); dy = X(:,2)+u(:,2);
        subplot(2,1,1); 
        trimesh( T, X(:,1), X(:,2), 'color', coul, 'LineStyle', ':' );
        hold on; trimesh( T, dx, dy, 'color', coul, 'LineStyle', '-' );

    % BEAM
    case 'Beam'
        X = model.Beam.X;
        h = 0.02*(max(X)-min(X));
        theta = u(:,2);
        u = u(:,1);
        du = diff(u) ./ diff(X); du = [du du]';
        dX = [X(1:end-1,1) X(2:end,1)]';
        xt = [ X-h*sin(theta) X+h*sin(theta) ];
        yt = [ u+h*cos(theta) u-h*cos(theta) ];
        subplot(2,1,1); hold on; plot( X, u, [coul style] );
        hold on; plot( xt', yt', coul );
        subplot(2,1,2); hold on; plot( dX(:), du(:), [coul style] );

    % NOT IMPLEMENTED YET
    otherwise
        error('code not implemented yet in plotTestCArl');
end
