function [ sol, out ] = CArl( model, coupling, solver )
% CARL to solve a system of models coupled in the Arlequin way
%
%  syntax: sol = CArl( model, coupling, solver )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       - 'code': code to be used to construct the stiffness matrices.
%                 Implemented: {'Comsol' 'HomeFE' 'ZacFE'}
%       - 'mesh': array with the field 'type' = 'FE' or 'discrete'
%                 In the case 'FE' it also contains fields 'X' and 'T',
%                   containing respectively the nodal coordinates and
%                   connectivity matrix.
%                 In the case 'discrete', it also contains fields 'X' and 
%                   'T', containing respectively the nodes of the grid and
%                   the closest nodes.
%       - other fields may be present depending on the type of code used
%
%  coupling: cell of structured array describing the coupling options, 
%       - 'models': 1*2 vector naming the models to be coupled (the number
%                   k indicates model{k}
%       - 'weight': structured array indicating the coupling domain and the
%                   way the weights should be computed. The arrays are
%             --'interior' and 'exterior', each with the fields
%                   * 'levelSet': set of points creating a closed surface
%                   (the last point is linked to the first one : it 
%                   has to be the same point). The dimension of the closed
%                   space has to be big enough
%                   * 'value' is the value of the weight for the FIRST
%                   model in coupling.models inside or outside 'levelSet',
%                   as indicated below. The 'value' corresponding to the
%                   SECOND model in that area is 1-'value'
%                'interior' means that the 'value' should be set for the
%                first model inside the 'levelSet' indicated, while
%                'exterior' sets the 'value' outside the model indicated.
%             --'type': 'linear' or 'constant' indicates how the weights
%               evolve between the values set inside the interior levelSet
%               and outside the exterior levelSet.
%             --'value': value at which the weight is set when
%               'type'='constant'
%       - 'mediator': structured array for the definition of the mediator 
%                     space. The fields are
%             --'type': 'deterministic' or 'stochastic' ('determinstic by
%                     default)
%             --'support': 1 or 2. indicates which model should be used as
%                     the support for the physical basis functions. Note
%                     that usually the coarser support should be used.
%       - 'operator': type of coupling operator used 'H1' or 'L2'.
%
%  solver: 'direct', 'FETI', 'MonteCarlo'
%
%  sol: structured array containing the solution, including
%       - 'models': cell of structured array containing the solution
%                   respective to each model independently
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% constants
Nm = length( model );    % number of models
Nc = length( coupling ); % number of coupling operations

% initialization
K = cell(Nm,1);
F = cell(Nm,1);
C1 = cell(Nc,1);
C2 = cell(Nc,1);
alpha = cell(Nc,2);
levelSet = cell(Nc,2);
c2m = zeros(Nc,2);

% reading the mesh in CArl format from exterior model
if  strcmp (model{1}.code,'Comsol')
        Tint1 = model{1}.femcomsol.mesh.t' ;
        Tint2 = model{2}.femcomsol.mesh.t' ;
        model{1}.mesh.X = model{1}.femcomsol.mesh.p' ;
        model{1}.mesh.T = Tint1(:,1:end-1) ;
        model{2}.mesh.X = model{2}.femcomsol.mesh.p' ;
        model{2}.mesh.T = Tint2(:,1:end-1) ;   
end

% construction of coupling operators
disp('Creating coupling matrices ...')
for i1 = 1:Nc

    % definition of models for this coupling
    couple = coupling{i1};
    c2m( i1, : ) = couple.models;
    model1 = model{ couple.models(1) };
    model2 = model{ couple.models(2) };
    
    % compute weights for each model
    [alpha{i1,1}, levelSet{i1,1}.ext, levelSet{i1,1}.int] = ...
                           ArlequinWeight( model1.mesh, couple.weight, 1 );
    [alpha{i1,2}, levelSet{i1,2}.ext, levelSet{i1,2}.int] = ...
                           ArlequinWeight( model2.mesh, couple.weight, 2 );
    
    % create intersection of meshes (for both representation and 
    % integration purposes)
    [ Int, Rep ] = MeshIntersect( model1, model2, levelSet{i1,1}, ...
                                                          levelSet{i1,2} );
    
    % definition of the mediator space
    Int.M = MediatorSpace( couple.mediator, Int, Rep ); 
    
    % construction of coupling operators
    C1{i1} = CouplingOperator( couple, model1.code, Int, Rep{1} );
    C2{i1} = CouplingOperator( couple, model2.code, Int, Rep{2} );
end

% construct stiffness and force matrices
disp('Creating stiffness matrices ...')
for i1 = 1:Nm
    model{i1}.alpha = prod( cat( 2, alpha{ find( c2m(:) == i1 ) } ), 2 );
    [ K{i1}, F{i1} ] = StiffnessMatrix( model{i1} );
end

% assemble sparse matrix system
disp('Assembling system ...')
[ K, F, opt ] = AssembleArlequin( K, F, C1, C2, model, coupling );
clear C1 C2;

% solve the Arlequin system
disp('Inverting system ...')
[ sol, out ] = SolveArlequin( K, F, solver, opt );

% preparing output
out = struct( 'models', {model}, ...
              'coupling', coupling, ...
              'K', K, ...
              'F', F, ...
              'u', out );
