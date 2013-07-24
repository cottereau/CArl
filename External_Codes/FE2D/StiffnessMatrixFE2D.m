function [K,F] = StiffnessMatrixFE2D( m, alpha )
% STIFFNESSMATRIXFE2D to construct the stiffness matrix corresponding to an
% elastic 2D problem by calling the external code FE2D
%
% syntax: [K,F] = StiffnessMatrixFE2D( model, alpha )
%
%    model: structured array containing the fields 'X', 'T', 'lambda', 'mu'
%           'load', 
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

% construction of the variable parameter field
warning('variable alpha not implemented yet')

% construction of Stiffness Matrix
m.mu = 100;
[ K, areas ] = stifness_matrixP1_2D_elasticity( m.T, m.X, m.lambda, m.mu );

% construction of Mass Matrix and load vector
M = mass_matrixP1_2D_elasticity( m.T, areas );
F = sum(M,2) .* reshape( m.load', numel(m.load), 1 );

% enforcing boundary conditions
N = 2*length(m.dirichlet);
I = sparse(1:N, [ (m.dirichlet-1)*2+1; 2*m.dirichlet ], 1, N, size(K,1) );
K = [ K I'; I sparse(N,N) ];
F = [ F; sparse(N,1) ];