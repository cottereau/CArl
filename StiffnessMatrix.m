function [ K, F, MC ] = StiffnessMatrix( model )
% STIFFNESSMATRIX to construct the basic stiffness matrix and force vector
% by calling an external code
%
% syntax: [ K, F, MC ] = StiffnessMatrix( model )
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
%  K, F: sparse matrices of stiffness and load
%  MC: realizations of the stiffness matrix when stochastic models are
%      considered (empty otherwise). When non-empty, K is then a zero
%      matrix
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
        [ K, F ] = StiffnessMatrixHomeFE( model.HomeFE, model.alpha );

    % MONTE CARLO VERSION OF HOME FE
    case 'MonteCarloHomeFE'
        stiff = model.HomeFE;
        property = stiff.property;
        Nmc = size( property, 3 );
        MC = cell( Nmc, 1 );
        for i1 = 1:Nmc
            stiff.property = property(:,:,i1);
            [ K, F ] = StiffnessMatrixHomeFE( stiff, model.alpha );
            MC{i1} = K;
        end
        K = sparse( size(K,1), size(K,2) );
        
    % 2D ELASTIC MODEL
    case 'FE2D'
        [ K, F ] = StiffnessMatrixFE2D( model.FE2D, model.alpha );

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
        K = sparse( x, y, K );
        F = sparse( z, 1, F );

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
        K = sparse( x, y, zeros(size(tmpK)) );
        F = sparse( z, 1, F );
        
        
    % unknown case
    otherwise
        error('this external code is not supported')
 
end

