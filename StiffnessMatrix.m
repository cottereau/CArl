function [ K, F ] = StiffnessMatrix( model )
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

% switch on the external code
switch model.code
    
    % HOMEFE
    case 'HomeFE'
        % modify the properties in each element according to alpha
        model.mesh = model.mesh.tri3;
        model.BC = model.HomeFE.BC;
        alpha = repmat(interp( model.alpha, model.mesh.incenters),[1 3]);
        model.property = model.HomeFE.property .* alpha;
        model.load = model.HomeFE.load .* alpha;

        % compute the modified value of the model property values
        [ x, y, K, z, F ] = StiffnessMatrixHomeFE( model );
                
    % MONTE CARLO VERSION OF HOME FE
    case 'MonteCarloHomeFE'

        % modify the properties in each element according to alpha
        model.mesh = model.mesh.tri3;
        alpha = repmat(interp( model.alpha, model.mesh.incenters),[1 3]);
        model.load = model.HomeFE.load .* alpha;
        property = model.HomeFE.property;
        model.BC = model.HomeFE.BC;
        model.property = [];

        % compute the modified value of the model property values
        Nmc = size( model.property, 3 );
        Ktot = cell( Nmc, 1 );
        for i1 = 1:Nmc
            model.property = property(:,:,i1) .* alpha;
            [ x, y, Ktot{i1}, z, F ] = StiffnessMatrixHomeFE( model );
        end
        K = Ktot{1};
                
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
K = struct( 'x', x, 'y', y, 'val', K );
F = struct( 'x', z, 'y', ones(size(z)), 'val', F );
if exist( 'Ktot', 'var' )
    K.MC = Ktot;
end
