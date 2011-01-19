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

% compute the modified value of the model property values
if numel(model.property)>1
    model.property = diag(model.alpha) * model.property;
else
    model.property = model.alpha * model.property;
end    
model.load = model.alpha .* model.load;

% switch on the external code
switch model.code
    
    % HOMEFE
    case 'HomeFE'
        % compute the modified value of the model property values
        [ x, y, K, z, F, k ] = StiffnessMatrixHomeFE( model );
        
    % ZacFE
    % based on the Rasmus Anthin Finite Element toolbox
    case 'ZacFE'
        % compute the modified value of the model property values
        [ x, y, K, z, F, k ] = StiffnessMatrixZacFE( model );
        
    % MONTE CARLO VERSION OF HOME FE
    case 'MonteCarloHomeFE'
        Nmc = size( model.property, 2 );
        Ktot = cell( Nmc, 1 );
        property = model.property;
        for i1 = 1:Nmc
            model.property = property(:,i1);
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
