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
%  | 'beamcomp2D'       | beam        | beam        | compression  |
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
%  | 'beamFEComp2D'     | beam        | elastic     | compression  |
%  | 'beamFE2D'         | beam        | elastic     |              |
%  | 'StoDetFE2D'       | elastic     | elastic (S) |              |
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
    type = '1d';
end
type = lower(type);

% load case parameters
if ~( strcmp(type,'1d') || strcmpi(type,'2d') || strcmpi(type,'mesh') )
    
    % get model and solve
    load( ['Tests/' type '.mat'] );
    sol = CArl( model, coupling, solver );
    
    % plot results
    figure; plotTestCArl( model{1}, sol{1}, 'k', '--x' )
    hold on; plotTestCArl( model{2}, sol{2}, 'b', '--o' );
    subplot(2,1,1); title( [type ' - primal'] ); box on; grid on
    subplot(2,1,2); title( [type ' - dual'] ); box on; grid on
    if exist( 'ref', 'var' )
        hold on; plotTestCArl( ref, ref.sol, 'r', ':' );
    end
    
    return
end

% selection of type of test
switch type
    
    case '1d'
        Test('zoom1D');
        Test('join1D');
        Test('force1D');
%        Test('stoDet1D');
%        Test('stoSto1D');
        Test('BeamComp2D');
        Test('BeamFEComp2D');
        Test('BeamFEComp2Dcoarse');
        Test('BeamFEComp2Dcoarse2');
        
    case '2d'
        Test('join2D');
        Test('zoom2D');
        Test('force2D');
%        Test('stoDet2D');
%        Test('stoSto2D');
        Test('FE2D');
        Test('Beam2D');
        Test('BeamFE2D');
%        Test('stoDetFE2D');
        
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
            u = mean(u,2);
            Y = X(:,2);
            X = X(:,1);
            subplot(2,1,1); hold on; trimesh( T, X, Y, u ); view(3)
            Ne = size(T,1);
            strain = zeros(Ne,2);
            for i1=1:Ne
                grad = [X(T(i1,:)) Y(T(i1,:)) ones(3,1)]\u(T(i1,:),:);
                strain(i1,:) = grad(1:2,:)';
            end
            Xc = mean(X(T),2); Yc = mean(Y(T),2);
            subplot(2,1,2); hold on; trimesh( T, X, Y, 'color', coul );
            hold on; quiver( Xc, Yc, strain(:,1), strain(:,2),'color',coul );
        end
        
    % ELASTIC - FE2D
    case 'FE2D'
        u = mean(u,3);
        X = model.FE2D.X(:,1);
        Y = model.FE2D.X(:,2);
        T = model.FE2D.T;
        dx = X+u(:,1); dy = Y+u(:,2);
        subplot(2,1,1); 
        trimesh( T, X, Y, 'color', coul, 'LineStyle', ':' );
        hold on; trimesh( T, dx, dy, 'color', coul, 'LineStyle', '-' );
        Ne = size(T,1);
        strain = zeros(Ne,4);
        for i1=1:Ne
            grad = [X(T(i1,:)) Y(T(i1,:)) ones(3,1)]\u(T(i1,:),:);
            strain(i1,:) = reshape( grad(1:2,:), 1, 4 );
        end
        % trace of strain operator
        strain = strain(:,1)+strain(:,4);
        Xc = mean(X(T),2); Yc = mean(Y(T),2);
        subplot(2,1,2); 
        trimesh( T, X, Y, 'color', coul, 'LineStyle', ':' );
        hold on; scatter( Xc, Yc, 50, strain, 'full' );
        colorbar

    % BEAM
    case 'Beam'
        B = model.Beam;
        X = B.X;
        h = 0.02*(max(X)-min(X));
        theta = u(:,3);
        v = u(:,2);
        u = u(:,1);
        dv = diff(v) ./ diff(X); dv = [dv dv]';
        dv = dv + repmat((theta(1:end-1)+theta(2:end))/2,[1 2])';
        q = B.kappa*B.h*1*B.young/2/(1+B.poisson) * dv;
        dX = [X(1:end-1,1) X(2:end,1)]';
        xt0 = [ X X ];
        yt0 = ones(size(X))*[ 1 -1 ]*h;
        xt = [ X+u+h*sin(theta) X+u-h*sin(theta) ];
        yt = [ v+h*cos(theta) v-h*cos(theta) ];
        subplot(2,1,1); hold on; plot( X, zeros(size(X)), [coul ':'] );
        hold on; plot( xt0', yt0', [coul ':'] );
        hold on; plot( X+u, v, [coul style] );
        hold on; plot( xt', yt', coul );
        subplot(2,1,2); hold on; plot( dX(:), q(:), [coul style] );

    % NOT IMPLEMENTED YET
    otherwise
        error('code not implemented yet in plotTestCArl');
end
