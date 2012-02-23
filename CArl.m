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
%  coupling : cell of structured array containing the description of each
%        of the couplings, with the fields
%       -'type' : either 'zoom' or 'join'
%       -'models' : 1*2 vector naming the models to be coupled (the number
%                  k indicates model{k}). When 'zoom' is used, the first
%                  model of this vector should be the coarse one.
%       -'mediator' : structured array for the definition of the mediator 
%                  space. The fields are
%             --'type': 'deterministic' or 'stochastic' ('deterministic by
%                     default)
%             --'support': 1 or 2. indicates which model should be used as
%                     the support for the physical basis functions. Note
%                     that usually the coarser support should be used.
%       -'operator' : type of coupling operator used ('H1' or 'L2').
%       -'LevelSet1'; 'LevelSet2' : definition of level-sets bounding the
%             coupling zone, given as a mesh in (X,T) format, with the
%             appropriate dimension
%       -'epsilon' : residual value of the non-proeminent model
%  solver: 'direct' or 'MonteCarlo'
%
%  NB: more complex coupling can be implemented 
%     (see function DefineClassicalCoupling.m)
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

if ~isfield( coupling{1}, 'c2m' )
    % construction of coupling operators
    disp('Creating coupling matrices ...')
    for i1 = 1:Nc
        
        % definition of models for this coupling
        couple = coupling{i1};
        c2m( i1, : ) = couple.models;
        model1 = model{ couple.models(1) };
        model2 = model{ couple.models(2) };
        
        % define level sets
        couple = DefineClassicalCoupling( couple, model1 );
        LSet1 = DefineLevelSet( model1.mesh.X, couple.weight1 );
        LSet2 = DefineLevelSet( model2.mesh.X, couple.weight2 );

        % compute weights for each model
        alpha1{i1} = ArlequinWeight( model1.mesh, couple.weight1, LSet1 );
        alpha2{i1} = ArlequinWeight( model2.mesh, couple.weight2, LSet2 );
        
        % create intersection of meshes (for both representation and
        % integration purposes)
        [Int,Rep] = MeshIntersect( model1.mesh, model2.mesh, LSet1, LSet2 );
        
        % definition of the mediator space
        Int.M = MediatorSpace( couple.mediator, Int, Rep );

        % construction of coupling operators
        C1{i1} = CouplingOperator( couple, model1.code, Int, Rep{1} );
        C2{i1} = CouplingOperator( couple, model2.code, Int, Rep{2} );
        
        % storing information for further use
        coupling{i1}.C1 = C1{i1};
        coupling{i1}.C2 = C2{i1};
        coupling{i1}.alpha1 = alpha1{i1};
        coupling{i1}.alpha2 = alpha2{i1};
        coupling{i1}.c2m = c2m(i1,:);
    end
else
    % construction of coupling operators
    disp('Loading coupling matrices ...')
    for i1=1:Nc
        C1{i1} = coupling{i1}.C1;
        C2{i1} = coupling{i1}.C2;
        alpha1{i1} = coupling{i1}.alpha1;
        alpha2{i1} = coupling{i1}.alpha2;
        c2m(i1,:) = coupling{i1}.c2m;
    end
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
              'coupling', {coupling}, ...
              'K', K, ...
              'F', F, ...
              'u', out, ...
              'opt', opt );
