function [ sol, out ] = CArl( Mdl, Cpl, solver, opt )
% CARL to solve a system of models coupled in the Arlequin manner
%
%  syntax: sol = CArl( model, coupling, solver, opt )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       -'code': code to be used to construct the stiffness matrices.
%               Implemented: 'Comsol', 'HomeFE', 'MonteCarloHomeFE','FE2D'
%       - other fields should be in a format appropriate for the code used
%               to return a stiffness matrix in sparse format
%       NB: fields 'mesh', 'alpha', 'K', 'F', 'Kmc', 'u' should not be used
%           because they are internally used by the program
%
%  coupling : cell of structured arrays containing the description of each
%        of the couplings, with the fields
%       -'models' : 1*2 vector naming the models to be coupled (the number
%                  k indicates model{k}).
%       -'code' : indicates the code that will compute the coupling
%                 matrices 'HomeFE', 'FE2D' or 'Beam'
%       -'mediator': 1 or 2. indicates which model should be used as the
%                    support for the physical basis functions. Note that
%                    usually the coarser support should be used.
%       -'operator' : type of coupling operator used ('H1' or 'L2').
%       -'levelSet' : definition of the coupling zone, given as a
%                     levelSet1D or levelSet object.
%       -'cplval' : value of the weight function in the coupling zone (cell
%                   of values for each model. When the value is a scalar, a
%                   constant function is used in the coupling area. When
%                   two values are given, the function is taken as linear
%                   between the two values given.
%       -'freeval' : cell of the (constant) value of the weight function
%                    for each model outside the coupling zone.
%
%  solver: 'direct' (default) or 'FETI' (not implemented yet)
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
%  sol: cell containing the solution for each model
%  out: structured array containing information on the computation
%
% developed at 
% Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579, 
% grande voie des vignes
% F-92295 Chatenay-Malabry
% FRANCE
% contact: regis.cottereau@ecp.fr

% initialization
if nargin<3 || isempty(solver); solver = 'direct'; end
if nargin<4; opt = []; end
[ Nm, Nc, opt ] = Initiatevalues( Mdl, Cpl, opt ) ;

% reading the code-dependant mesh into CArl format (TRI6)
for i1 = 1:Nm
    Mdl{i1} = ReadCodeMesh( Mdl{i1} );
end

% construction of coupling operators
disp('Creating coupling matrices ...')
for i1 = 1:Nc
    
    % check for recomputing
    if opt.recomputeC(i1)
    
        % definition of models for this coupling
        m1 = Mdl{ Cpl{i1}.models(1) };
        m2 = Mdl{ Cpl{i1}.models(2) };
        
        % define coupling and free volumes
        Cpl{i1} = DefineDomains( Cpl{i1}, m1, m2 );

        % compute weights for each model
        Cpl{i1}.alpha{1} = ArlequinWeight( Cpl{i1}, 1, m1 );
        Cpl{i1}.alpha{2} = ArlequinWeight( Cpl{i1}, 2, m2 );

        % create intersection of meshes (for both representation and
        % integration purposes)
        [Int,Rep] = MeshIntersect( m1, m2, Cpl{i1}.domain );

        % definition of the mediator space
        Int.M = MediatorSpace( Cpl{i1}.mediator, Rep );
        
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
        Mdl{i1}.alpha = CondensateAlpha( i1, Cpl );

        % construct stiffness and force matrices
        [ Mdl{i1}.K, Mdl{i1}.F, Mdl{i1}.Kmc ] = StiffnessMatrix( Mdl{i1} );
        
    end
end

% assemble sparse matrix system
disp('Assembling and inverting system ...')
[ K, F, opt, Kmc ] = AssembleArlequin( Mdl, Cpl );

% solving system
[ Mdl, Cpl ] = SolveArlequin( K, F, Mdl, Cpl, solver, opt, Kmc );

% preparing output
[ sol, Mdl ] = ArlequinOutput( Mdl );
out = struct( 'model', {Mdl}, ...
              'coupling', {Cpl}, ...
              'K', K, ...
              'F', F, ...
              'opt', opt );
