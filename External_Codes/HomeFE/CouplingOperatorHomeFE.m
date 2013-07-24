function C = CouplingOperatorHomeFE( operator, mesh, opt )
% COUPLINGOPERATORHOMEFE to construct the Arlequin coupling matrix 
% by calling a home-made FE code
%
% syntax: C = CouplingOperatorHomeFE( operator, mesh, opt )
%
%    operator: 'H1' or 'L2' [string]
%    mesh    : mesh structure [INT3 or TRI6 object]
%    opt     : structured array containing field 'kappa' (only used with
%              'H1' operator
%
%    C: the output matrix is in sparse format
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% computation of MassMatrix
m = struct( 'mesh', struct( 'X', mesh.Points, ...
                            'T', mesh.ConnectivityList), ...
            'property', ones(size(mesh)), ...
            'load', zeros(size(mesh)), ...
            'BC', [] );
[ x, y, C ] = MassMatrixHomeFE( m );
C = sparse( x, y, C );

% choice of the coupling operator
switch operator
    
    % L2 coupling
    case 'L2'            
        
    % H1 coupling
    case 'H1'
        [ z, k, K ] = StiffnessMatrixHomeFE( m );
        K = sparse( z, k, K );
        C = opt.kappa*K + C;
        
    % unknown coupling operator
    otherwise
        error('unknown coupling operator')
end

