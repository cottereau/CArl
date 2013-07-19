function [ x, y, C ] = CouplingOperatorHomeFE( operator, mesh, opt )
% COUPLINGOPERATORHOMEFE to construct the Arlequin coupling matrix 
% by calling a home-made FE code
%
% syntax: [ x, y, C ] = CouplingOperatorHomeFE( model coupling, n )
%
%  model: structured array containing the information relative to the
%         model, in particular:
%       - 'type': describes the type of representation used for the
%                 geometry. This is used for all interpolation purposes,
%                 in particular for the definition of weight functions
%                 Implemented: {'FE' 'discrete'}
%       - 'code': code to be used to construct the stiffness matrices.
%                 Implemented: {'Comsol' 'HomeFE'}
%       - 'mesh': array dependent on 'type'.
%  coupling: cell of structured array describing the coupling options, in
%            particular
%       - 'mesh': indicates the subdomain of the mesh that corresponds to
%                 the coupling domain
%       - 'mediator': definition of the mediator space: 'default' or an
%                     integer indicating the optional choices (not always
%                     available)
%       - 'operator': coupling operator 'H1' 'L2'. The definition of these
%                     spaces may be dependent on the type of coupling
%  n: indicates what model is being considered [1 or 2]
%
%  output: the format is that of sparse matrices. The coupling matrix
%          is such that, schematically:
%               CouplingMatrix( x, y ) = C
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% constants
s = size(mesh);
Nne = s(1);
ne = s(2);

% computation of MassMatrix
m = struct( 'mesh', struct( 'X', mesh.Points, ...
                            'T', mesh.ConnectivityList), ...
            'property', ones(Nne,ne), ...
            'load', zeros(Nne,ne), ...
            'BC', [] );
[ x, y, C ] = MassMatrixHomeFE( m );

% choice of the coupling operator
switch operator
    
    % L2 coupling
    case 'L2'            
        
    % H1 coupling
    case 'H1'
        [ z, k, K ] = StiffnessMatrixHomeFE( m );
        x = [ z; x ];
        y = [ k; y ];
        C = [ opt.kappa*K; C ];
        
    % unknown coupling operator
    otherwise
        error('unknown coupling operator')
end

