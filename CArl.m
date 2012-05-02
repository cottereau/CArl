function [ sol, out ] = CArl( Mdl, Cpl, solver, opt )
% CARL to solve a system of models coupled in the Arlequin way
%
%  syntax: sol = CArl( model, coupling, solver, opt )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       -'code': code to be used to construct the stiffness matrices.
%               Implemented: {'Comsol' 'HomeFE' 'MonteCarloHomeFE}
%       - other fields should be in a format appropriate for the code used
%               to return a stiffness matrix in sparse format
%
%  coupling : cell of structured arrays containing the description of each
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
%
%  solver: 'direct' or 'MonteCarlo'
%
%  opt: structured array with options. Possible options include
%       -'recomputeK' logical array of size Nm*1 where Nm is the number of
%             models. Indicates, for each model, whether the stiffness
%             matrix should be recomputed. If .false. the fields 'K' and 
%             'F' should be defined in corresponding cell of input 
%             structured array 'model'. Default is .true.
%       -'recomputeC' logical array of size Nc*1 where Nc is the number of
%             couplings. Indicates, for each coupling, whether the
%             corresponding elements should be recomputed. If .false. the
%             fields 'C1', 'C2', 'alpha1', and 'alpha2' should be defined
%             in the corresponding cell of input structured array 
%             'coupling'. Default is .true.
%       -'gerr' to set a threshold distance for points to be considered
%             geometrically equal for some of the meshing operations.
%             Default is 1e-9
%       -'kappa' to set a value or the ratio between H1 and L2 elements of
%             the H1 norm for the coupling operator. Default is 1e-3
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
if nargin<4; opt = []; end
[ Nm, Nc, opt ] = Initiatevalues( Mdl, Cpl, opt ) ;

% reading the code-dependant mesh into CArl format
for i1 = 1:Nm
    Mdl{i1}.mesh = ReadCodeMesh( Mdl{i1} );
end

% construction of coupling operators
disp('Creating coupling matrices ...')
for i1 = 1:Nc
    
    % check for recomputing
    if opt.recomputeC(i1)
    
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
    
    % check for recomputing
    if opt.recomputeK(i1)
        
        % condensate alpha functions for each model
        Mdl{i1}.alpha = CondensateAlpha( i1, Mdl{i1}, Cpl );
        
        % construct stiffness and force matrices
        [ Mdl{i1}.K, Mdl{i1}.F ] = StiffnessMatrix( Mdl{i1} );
        
    end
end

% assemble sparse matrix system
disp('Assembling and inverting system ...')
[ K, F, opt ] = AssembleArlequin( Mdl, Cpl );

% solving system
[ Mdl, Cpl ] = SolveArlequin( K, F, Mdl, Cpl, solver, opt );

% preparing output
[ sol, Mdl ] = ArlequinOutput( Mdl );
out = struct( 'model', {Mdl}, ...
              'coupling', {Cpl}, ...
              'K', K, ...
              'F', F, ...
              'opt', opt );
