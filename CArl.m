function [ sol, out ] = CArl( Mdl, Cpl, solver )
% CARL to solve a system of models coupled in the Arlequin way
%
%  syntax: sol = CArl( model, coupling, solver )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       -'code': code to be used to construct the stiffness matrices.
%               Implemented: {'Comsol' 'HomeFE' 'MonteCarloHomeFE}
%       - other fields should be in a format appropriate for the code used
%               to a return stiffness matrix
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
%  sol: cell containing (a part of) the solution for each model
%  out: structured array containing information on the computation
%
% developed at 
% Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579, 
% grande voie des vignes
% F-92295 Chatenay-Malabry
% FRANCE
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% initialization
[ Nm, Nc, opt ] = Initiatevalues( Mdl, Cpl ) ;

% reading the code-dependant mesh into CArl format
for i1 = 1:Nm
    Mdl{i1}.mesh = ReadCodeMesh( Mdl{i1} );
end

if ~isfield( Cpl{1}, 'C1' )
    % construction of coupling operators
    disp('Creating coupling matrices ...')
    for i1 = 1:Nc
        
        % definition of models for this coupling
        m1 = Mdl{ Cpl{i1}.models(1) };
        m2 = Mdl{ Cpl{i1}.models(2) };
        
        % define level sets
        Cpl{i1} = DefineClassicalCoupling( Cpl{i1}, m1 );
        LSet1 = DefineLevelSet( m1.mesh.X, Cpl{i1}.weight1 );
        LSet2 = DefineLevelSet( m2.mesh.X, Cpl{i1}.weight2 );

        % compute weights for each model
        Cpl{i1}.alpha1 = ArlequinWeight( m1.mesh, Cpl{i1}.weight1, LSet1, opt );
        Cpl{i1}.alpha2 = ArlequinWeight( m2.mesh, Cpl{i1}.weight2, LSet2, opt );
        
        % create intersection of meshes (for both representation and
        % integration purposes)
        [Int,Rep] = MeshIntersect( m1.mesh, m2.mesh, LSet1, LSet2, opt );
        
        % definition of the mediator space
        Int.M = MediatorSpace( Cpl{i1}.mediator, Int, Rep );

        % construction of coupling operators
        Cpl{i1}.C1 = CouplingOperator( Cpl{i1}, Int, Rep{1}, opt );
        Cpl{i1}.C2 = CouplingOperator( Cpl{i1}, Int, Rep{2}, opt );
                
    end
end

% construct stiffness and force matrices
disp('Creating stiffness matrices ...')
for i1 = 1:Nm
    % condensate alpha functions for each model
    alpha = CondensateAlpha( i1, Mdl{i1}, Cpl );
    % compute stiffness and force matrices
    [ Mdl{i1}.K, Mdl{i1}.F ] = StiffnessMatrix( Mdl{i1}, alpha );
end

% assemble sparse matrix system
disp('Assembling system ...')
[ K, F, opt ] = AssembleArlequin( Mdl, Cpl );
clear C1 C2;

% solve the Arlequin system
disp('Inverting system ...')
[ sol, out ] = SolveArlequin( K, F, solver, opt );

% preparing output
out = struct( 'models', {Mdl}, ...
              'coupling', {Cpl}, ...
              'K', K, ...
              'F', F, ...
              'u', out, ...
              'opt', opt );
