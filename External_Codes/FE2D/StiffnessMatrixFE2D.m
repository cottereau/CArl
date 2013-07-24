function C = StiffnessMatrixFE2D( operator, mesh, opt )
% COUPLINGOPERATORFE2D to construct the Arlequin coupling matrix 
% by calling a the elastic code FE2D
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
%
% uses routines by Talal Rahman and Jan Valdman downloadable at 
% http://www.mathworks.in/matlabcentral/fileexchange/
% 27826-fast-assembly-of-stiffness-and-matrices-in-finite-element-method 

% construction of Stiffness Matrix
K = stifness_matrixP1_2D_elasticity( mesh.ConnectivityList, ...
                                                     mesh.Points, 0, 1/2 );

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