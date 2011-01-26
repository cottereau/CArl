function [ K, F ] = StiffnessMatrix( model, alpha )
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
        alpha1 = polyvalN( alpha, model.mesh.X(model.mesh.T(:,1)) );
        alpha2 = polyvalN( alpha, model.mesh.X(model.mesh.T(:,2)) );
        alpha = [alpha1 alpha2];
        model.property = model.property .* alpha;
        model.load = model.load .* alpha;
        % compute the modified value of the model property values
        [ x, y, K, z, F, k ] = StiffnessMatrixHomeFE( model );
        
    % ZacFE
    % based on the Rasmus Anthin Finite Element toolbox
    case 'ZacFE'
        % compute the modified value of the model property values
        [ x, y, K, z, F, k ] = StiffnessMatrixZacFE( model );
        
    % MONTE CARLO VERSION OF HOME FE
    case 'MonteCarloHomeFE'
        alpha1 = polyvalN( alpha, model.mesh.X(model.mesh.T(:,1)) );
        alpha2 = polyvalN( alpha, model.mesh.X(model.mesh.T(:,2)) );
        alpha = [alpha1 alpha2];
        model.load = model.load .* alpha;
        Nmc = size( model.property, 3 );
        Ktot = cell( Nmc, 1 );
        property = model.property;
        for i1 = 1:Nmc
            model.property = property(:,:,i1) .* alpha;
            [ x, y, Ktot{i1}, z, F, k ] = StiffnessMatrixHomeFE( model );
        end
        K = Ktot{1};
        
    % MONTE CARLO VERSION OF ZAC FE
    case 'MonteCarloZacFE'
        Nmc = size( model.property, 2 );
        Ktot = cell( Nmc, 1 );
        property = model.property;
        for i1 = 1:Nmc
            model.property = property(:,i1);
            [ x, y, Ktot{i1}, z, F, k ] = StiffnessMatrixZacFE( model );
        end
        K = Ktot{1};
        
    % COMSOL
    case 'Comsol'
        error( 'COMSOL is not supported yet' );
        
    % unknown case
    otherwise
        error('this external code is not supported')
        
end

% output
K = struct( 'x', x, 'y', y, 'val', K );
F = struct( 'x', z, 'y', k, 'val', F );
if exist( 'Ktot', 'var' )
    K.MC = Ktot;
end
