function [ K, F, MC ] = StiffnessMatrix( model )
% STIFFNESSMATRIX to construct the basic stiffness matrix and force vector
% by calling an external code
%
% syntax: [ x, y, K, z, F ] = StiffnessMatrix( model )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       - 'code': code to be used to construct the stiffness matrices. For
%                 now, only 'Comsol' and 'HomeFE' available
%       - 'mesh': array dependent on the type of model.
%       - 'BC': structured array dependant on the type of models and
%               describing the boundary conditions
%       - 'load': structured array dependant on the type of models and
%                 describing the loads
%
%  output: the format is that of sparse matrices. The matrix of stiffness
%          and the vector of force are such that, schematically:
%               Stiffness( x, y ) = K
%               Force( z ) = F
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% Initialization
MC = [];

% switch on the external code
switch model.code
    
    % HOMEFE
    case 'HomeFE'
        stiff = model.HomeFE;
        alpha = interp(model.alpha,incenters(model.mesh));
        alpha = repmat(alpha',size(stiff.mesh.T,2),1)';
        stiff.property = stiff.property .* alpha;
        stiff.load = stiff.load .* alpha;   
        keyboard
        [ x, y, K, z, F ] = StiffnessMatrixHomeFE( stiff );
        
    % MONTE CARLO VERSION OF HOME FE
    case 'MonteCarloHomeFE'
        
        % modify the properties in each element according to alpha
        stiff = model.HomeFE;
        property = stiff.property;
        Nmc = size( property, 3 );
        alpha = interp(model.alpha,incenters(model.mesh));
        alpha = repmat(alpha',size(stiff.mesh.T,2),1)';
        stiff.load = stiff.load .* alpha;
        MC = cell( Nmc, 1 );
        for i1 = 1:Nmc
            stiff.property = property(:,:,i1) .* alpha;
            [ x, y, tmpK, z, F ] = StiffnessMatrixHomeFE( stiff );
            MC{i1} = sparse( x, y, tmpK );
        end
        K = zeros(size(tmpK));
        
    % TIMOSCHENKO BEAM MODEL
    case 'Beam'
        keyboard
        inc = mean(model.mesh.X(model.mesh.T),2);
        alpha = interp(model.alpha,[inc;zeros(1,length(inc))]', ...
                                   [inc;zeros(1,length(inc))]');
        model.load = model.HomeFE.load .* repmat(alpha,[1,2,3]);
        model.property = model.property .* repmat(alpha,1,2);
        [ x, y, K, z, F ] = Timostiff( model );
        
        % compute the modified value of the model property values
        Nmc = size( model.HomeFE.property, 3 );
        MC = cell( Nmc, 1 );
        for i1 = 1:Nmc
            model.property = property(:,:,i1) .* repmat(alpha,1,2);
            [ x, y, tmpK, z, F ] = Timostiff( model );
            MC{i1} = sparse( x, y, tmpK );
        end
        K = zeros(size(tmpK));
        
    % 2D ELASTIC MODEL
    case 'FE2D'
        % modify the properties in each element according to alpha
        model.mesh = model.mesh.tri3;
        x = model.mesh.X(:,1); y = model.mesh.X(:,2);
        Ne = size(model.mesh.Triangulation,1);
        X = [reshape(x(model.mesh.Triangulation'),3*Ne,1) reshape(y(model.mesh.Triangulation'),3*Ne,1)];
        inc = model.mesh.incenters;
        xe = repmat( inc(:,1)', [3 1]);
        ye = repmat( inc(:,2)', [3 1]);
        Xe = [xe(:) ye(:)];
        alpha = interp(model.alpha,X,Xe);
        model.load = model.HomeFE.load .* repmat(reshape(alpha,3,Ne)',[1,1,2]);
        property = model.HomeFE.property;

        % compute the modified value of the model property values
        Nmc = size( model.HomeFE.property, 3 );
        MC = cell( Nmc, 1 );
        for i1 = 1:Nmc
            model.property = property(:,:,i1) .* reshape(alpha,3,Ne)';
            [ x, y, tmpK, z, F ] = K_2D_elasticity( model );
            MC{i1} = sparse( x, y, tmpK );
        end
        K = zeros(size(tmpK));
        
    % COMSOL
    case 'Comsol'
        wd = pwd;
        eval( ['cd ' model.meshpath ';' ]);
        eval( ['out = ' model.meshfile ';' ]);
        eval( ['out = ' model.solvefile '(out);' ] );
        eval( ['cd ' wd ';' ]);
        K = mphmatrix( out, 'sol1', 'Out', {'L','K'} );
        info = mphxmeshinfo( out,'solname', 'sol1' );
        ind = double( info.dofs.nodes+1 );
        [ z, ~, F ] = find( K.L ); z = ind(z);
        [ x, y, K ] = find( K.K ); x = ind(x); y = ind(y);
        
    % unknown case
    otherwise
        error('this external code is not supported')
 
end

% output
K = sparse( x, y, K );
F = sparse( z, 1, F );
