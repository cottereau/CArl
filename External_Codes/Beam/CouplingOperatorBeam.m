function C = CouplingOperatorBeam( operator, mesh, opt )
% COUPLINGOPERATORBEAM to construct the Arlequin coupling matrix 
% by calling a beam code
%
% syntax: C = CouplingOperatorBeam( operator, mesh, opt )
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
%
% uses routines by Manuel Diaz and A. Ferreira downloadable at 
% https://github.com/wme7/aero-matlab/tree/master/FEM/Timoshenko_beam

% constants
X = mesh.Points;
T = mesh.ConnectivityList;    

% construction of stiffness and mass Matrix
[ K, ~, C ] = formStiffnessMassTimoshenkoBeam( 2*mesh.Nn, mesh.Ne, ...
                                       T, mesh.Nn, X, eye(2), 1, 1, 1, 1 );

% adding coupling for the longitudinal displacement
model  = struct( 'mesh', struct('X',X,'T',T), 'property', ones(size(T)));
Cn = MassMatrixHomeFE( model );
Kn = StiffnessMatrixHomeFE( model );
O = zeros( size(Kn,1), size(K,2) );
K = [Kn O; O' K];
C = [Cn O; O' C];

% choice of the coupling operator
switch operator
    
    % L2 coupling
    case 'L2'
        
    % H1 coupling
    case 'H1'
        C = opt.kappa*K + C;

    % unknown coupling operator
    otherwise
        error('unknown coupling operator')
end

% sparse matrix
C = sparse(C);
