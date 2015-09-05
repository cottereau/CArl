function [K,F] = StiffnessMatrixFE2D( m, alpha )
% STIFFNESSMATRIXFE2D to construct the stiffness matrix corresponding to an
% elastic 2D problem by calling the external code FE2D
%
% syntax: [K,F] = StiffnessMatrixFE2D( model, alpha )
%
%    model: structured array containing the fields 'X', 'T', 'lambda', 'mu'
%           'load'
%    alpha: weight function [discontinuous] 
%
%    K,F: stiffness and force matrices in sparse format
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr
%
% uses routines by Talal Rahman and Jan Valdman downloadable at 
% http://www.mathworks.in/matlabcentral/fileexchange/
% 27826-fast-assembly-of-stiffness-and-matrices-in-finite-element-method 

% constants
N = 2*size(m.X,1);

% construction of the variable parameter field
x = m.X(:,1);
y = m.X(:,2);
if nargin>1
    alpha = interp( alpha, [ mean(x(m.T'),1); mean(y(m.T'),1)]' );
    alpha = repmat(alpha,[1 3]);
else
    alpha = ones(size(m.T));
end

% construction of Stiffness Matrix
if isfield(m,'cubic') && m.cubic
    K = stifness_matrixP1_2D_cubic( m.T, m.X, m.c1, m.c4, m.c12, alpha, m.theta );
else
    K = stifness_matrixP1_2D_elasfluc( m.T, m.X, m.lambda, m.mu, alpha );
end

% construction of load vector
F = sparse( N, 1 );
if isfield(m,'load') && any(m.load(:)~=0)
    error('bulk load not enforced yet in FE2D')
end
%    F = mass_matrixP1_2D_elasticity( m.T , S ) * reshape( m.load', N, 1 );

% add boundary conditions
if isfield(m,'BC')&&~isempty( m.BC )

    % constants
    type = m.BC.type;
    nodes = m.BC.nodes;
    val = m.BC.value;
    Nn = length(x);

    % Dirichlet Boundary Conditions in longitudinal displacement
    ind = type == 'U';
    Nbc = nnz(ind);
    BC = sparse( 2*nodes(ind)-1, 1:Nbc, 1, 2*Nn, Nbc );
    FBC = sparse( 1:Nbc, 1, val(ind), Nbc, 1 );

    % Dirichlet Boundary Conditions in transverse displacement
    ind = type == 'V';
    Nbc = nnz(ind);
    BC = [BC sparse( 2*nodes(ind), 1:Nbc, 1, 2*Nn, Nbc )];
    FBC = [FBC; sparse( 1:Nbc, 1, val(ind), Nbc, 1 )];

    % Neumann Boundary Conditions in longitudinal Force
    ind = type == 'N';
    Nbc = nnz(ind);
    if Nbc>0
        nod = 2*node(ind)-1;
        F(nod) = F(nod) + sparse( 1:Nbc, 1, val(ind), Nbc, 1 );
    end

    % Neumann Boundary Conditions in transverse Force
    ind = type == 'Q';
    Nbc = nnz(ind);
    if Nbc>0
        nod = 2*node(ind);
        F(nod) = F(nod) + sparse( 1:Nbc, 1, val(ind), Nbc, 1 );
    end
    
    % construct new matrices
    K = [ K BC; BC' sparse(size(BC,2),size(BC,2)) ];
    F = [ F; FBC ]; 
    
end
