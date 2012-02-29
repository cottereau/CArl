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
Nni = size(mesh.X,1);
[Nne,ne] = size(mesh.Triangulation);

% loop on the columns of the matrix and creation of a set of forces
f = zeros(Nne,ne,Nni);
for i1=1:Nni
    for i2 = 1:ne
        f( mesh.Triangulation(:,i2) == i1,i2,i1 ) = 1;
    end
end

m = struct( 'mesh', mesh, ...
            'property', ones(Nne,ne), ...
            'load', f, ...
            'BC', [] );

[ x, y, K, z, F, k ] = StiffnessMatrixHomeFE( m );

% choice of the coupling operator
switch operator
    
    % L2 coupling
    case 'L2'            
        x = z;
        y = k;
        C = F;
        
    % H1 coupling
    case 'H1'
        x = [ x; z ];
        y = [ y; k ];
        C = [ opt.kappa*K; F ];
        
    % unknown coupling operator
    otherwise
        error('unknown coupling operator')
end

