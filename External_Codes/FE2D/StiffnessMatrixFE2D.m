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
x = m.X(:,1);
y = m.X(:,2);
alpha = interp( alpha, [ mean(x(m.T),2) mean(y(m.T),2)] );
alpha = repmat(alpha,[1 3]);

% construction of Stiffness Matrix
% [ K, areas ] = stifness_matrixP1_2D_elasticity( m.T, m.X, m.lambda, m.mu );
K = stifness_matrixP1_2D_elasfluc( m.T, m.X, m.lambda, m.mu, alpha );

% construction of load vector
F = sparse( size(K,1), 1 );
if any(m.load~=0)
    error('bulk load not enforced yet in FE2D')
end
if ~isfield(m,'neumann')
    m.neumann = [];
end
if ~isempty(m.neumann)
    edges = m.neumann(:,1:2);
    load = m.neumann(:,3:4);
    d = m.X(edges(:,2),:) - m.X(edges(:,1),:);
    n = [-d(:,2) d(:,1)];
    load = load.*[n(:,1) d(:,1)] + load.*[n(:,2) d(:,2)];
    F(2*edges(:,1)-1) = F(2*edges(:,1)-1) + load(:,1)/2;
    F(2*edges(:,2)-1) = F(2*edges(:,2)-1) + load(:,1)/2;
    F(2*edges(:,1)) = F(2*edges(:,1)) + load(:,2)/2;
    F(2*edges(:,2)) = F(2*edges(:,2)) + load(:,2)/2;
end

% enforcing boundary conditions
N = 2*length(m.dirichlet);
I = sparse(1:N, [ (m.dirichlet-1)*2+1; 2*m.dirichlet ], 1, N, size(K,1) );
K = [ K I'; I sparse(N,N) ];
F = [ F; sparse(N,1) ];
