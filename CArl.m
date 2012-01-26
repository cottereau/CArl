function [ sol, out ] = CArl( model, coupling, solver )
% CARL to solve a system of models coupled in the Arlequin way
%
%  syntax: sol = CArl( model, coupling, solver )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       -'code': code to be used to construct the stiffness matrices.
%                Implemented: {'Comsol' 'HomeFE' 'ZacFE'}
%       - other fields should be in format appropriate for the code used
%
%  coupling: cell of structured array describing the coupling options, 
%       -'models': 1*2 vector naming the models to be coupled (the number
%                  k indicates model{k}
%       -'weight1': structured array defining the relevant information for
%                 model{coupling.models(1)} and containing the fields
%            --'ext': external geometrical limit of the coupling zone given 
%                 as a mesh in (X,T) format, with one dimension less than 
%                 the volume mesh model{coupling.models(1)}.mesh
%            --'int': internal geometrical limit of the coupling zone given
%                 in the same format as ext. The curve int should be 
%                 located inside the curve ext.
%              NB: in 1D, the level sets should include +/-Inf to give
%              sense to 'interior' and 'exterior' definitions
%            --'value': value of the weight function inside the coupling
%                       zone, given as a vector of coefficients of a
%                       polynomial (a scalar for a constant weight, a 2*1
%                       vector for a linear weight, etc ...)
%            --'extvalue': value of the weight outside the exterior curve
%            --'intvalue': value of the weight inside the interior curve
%       -'weight2': same as 'weight1' for the other model. When not
%                   defined, weight2=weight1, except for the values that
%                   are taken such that (value2+value1=1), and int and ext
%                   that are inversed.
%       -'mediator': structured array for the definition of the mediator 
%                    space. The fields are
%             --'type': 'deterministic' or 'stochastic' ('deterministic by
%                     default)
%             --'support': 1 or 2. indicates which model should be used as
%                     the support for the physical basis functions. Note
%                     that usually the coarser support should be used.
%       -'operator': type of coupling operator used ('H1' or 'L2').
%
%  solver: 'direct', 'FETI', 'MonteCarlo'
%
%  sol: structured array containing the solution, including
%       - 'models': cell of structured array containing the solution
%                   respective to each model independently
%
% developed at 
% Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579, 
% grande voie des vignes
% F-92295 Chatenay-Malabry
% FRANCE
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
alpha1 = cell(Nc,1);
alpha2 = cell(Nc,1);
c2m = zeros(Nc,2);

% reading the code-dependant mesh into CArl format
for i1 = 1:Nm
    model{i1}.mesh = ReadCodeMesh( model{i1} );
end

% construction of coupling operators
disp('Creating coupling matrices ...')
for i1 = 1:Nc
    % definition of models for this coupling
    couple = coupling{i1};
    c2m( i1, : ) = couple.models;
    model1 = model{ couple.models(1) };
    model2 = model{ couple.models(2) };
    
    % define level sets
    LSet1 = DefineLevelSet( model1.mesh.X, couple.weight1 );
    LSet2 = DefineLevelSet( model2.mesh.X, couple.weight2 );

    % compute weights for each model
    alpha1{i1} = ArlequinWeight( model1.mesh, couple.weight1, LSet1 );
    alpha2{i1} = ArlequinWeight( model2.mesh, couple.weight2, LSet2 );

    % create intersection of meshes (for both representation and 
    % integration purposes)
    [ Int, Rep ] = MeshIntersect( model1.mesh, model2.mesh, LSet1, LSet2 );

    % definition of the mediator space
    Int.M = MediatorSpace( couple.mediator, Int, Rep ); 
    
    % construction of coupling operators
    C1{i1} = CouplingOperator( couple, model1.code, Int, Rep{1} );
    C2{i1} = CouplingOperator( couple, model2.code, Int, Rep{2} );
end

% construct stiffness and force matrices
disp('Creating stiffness matrices ...')
for i1 = 1:Nm
    % condensate alpha functions for each model
    alpha = CondensateAlpha( i1, model{i1}, c2m, [alpha1 alpha2] );
    % compute stiffness and force matrices
    [ K{i1}, F{i1} ] = StiffnessMatrix( model{i1}, alpha );
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
              'u', out, ...
              'opt', opt );
