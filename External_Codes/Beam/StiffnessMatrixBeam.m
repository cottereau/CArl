function [ K, F ] = StiffnessMatrixBeam( model, alpha )
% STIFFNESSMATRIXBEAM to construct the stiffness matrix corresponding to a
% Timochenko beam problem by calling an external code
%
% syntax: [K,F] = StiffnessMatrixBeam( model, alpha )
%
%    model: structured array containing the fields 'X', 'T', 'young',
%           'poisson', 'h', 'P' (optional), 'BC' (optional)
%    mesh    : mesh structure [INT3 or TRI6 object]
%
%    K, F: the output matrices are in sparse format
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr
%
% uses routines by Manuel Diaz and A. Ferreira downloadable at 
% https://github.com/wme7/aero-matlab/tree/master/FEM/Timoshenko_beam

% geometrical parameters
b = 1;
h = model.h;
S = b*h;
I = b*h^3/12;

% mechanical parameters
E = model.young;
G = model.young/2/(1+model.poisson);
C = [ E*I; 0; 0; S*G];
P = model.P;

% mesh
X = model.X;
T = model.T;
Nn = size(X,1);
Ne = size(model.T,1);

% multiply properties and load by alpha (when present)
if nargin==2 && ~isempty(alpha)
    Xc = mean( X(T), 2 );
    alphac = interp( alpha, Xc );
    C = reshape( C * alphac', 2, 2, Ne );
    P = P * alphac;
else
    C = repmat( reshape(C,2,2), [1 1 Ne] );
    P = P * ones(Ne,1);
    alpha = [];
end

% compute stiffness matrix
[ K, F ] = StiffnessMassTimoshenkoBeamVar(2*Nn,Ne,T,Nn,X,C,P,1,I,h);
K = sparse(K);
F = sparse(F);

% add longitudinal stiffness
modelN = struct( 'mesh', struct('T',T,'X',X), ...
                'property', E*ones(size(T)) );
[ Kn, Fn ] = StiffnessMatrixHomeFE( modelN, alpha );

% condensate
O = sparse(size(Kn,1),size(K,2));
K = [Kn O;O' K];
F = [Fn; F];

% add boundary conditions
if isfield(model,'BC')&&~isempty( model.BC )

    % Dirichlet Boundary Conditions in longitudinal displacement
    ind = find(model.BC.type == 'U');
    Nbc = length(ind);
    BC = sparse( model.BC.nodes(ind), 1:Nbc, 1, 3*Nn, Nbc );
    FBC = sparse( 1:Nbc, 1, model.BC.value(ind), Nbc, 1 );

    % Dirichlet Boundary Conditions in transverse displacement
    ind = find(model.BC.type == 'V');
    Nbc = length(ind);
    BC = [BC sparse( model.BC.nodes(ind)+Nn, 1:Nbc, 1, 3*Nn, Nbc )];
    FBC = [FBC; sparse( 1:Nbc, 1, model.BC.value(ind), Nbc, 1 )];

    % Dirichlet Boundary Conditions in rotation
    ind = find(model.BC.type == 'R');
    Nbc = length(ind);
    BC = [BC sparse( model.BC.nodes(ind)+2*Nn, 1:Nbc, 1, 3*Nn, Nbc )];
    FBC = [FBC; sparse( 1:Nbc, 1, model.BC.value(ind), Nbc, 1 )];

    % Neumann Boundary Conditions in longitudinal Force
    ind = find(model.BC.type == 'N');
    Nbc = length(ind);
    if Nbc>0
        val = model.BC.value(ind);
        ind = model.BC.node(ind);
        F(ind) = F(ind) + sparse( 1:Nbc, 1, val, Nbc, 1 );
    end

    % Neumann Boundary Conditions in transverse Force
    ind = find(model.BC.type == 'Q');
    Nbc = length(ind);
    if Nbc>0
        val = model.BC.value(ind);
        ind = model.BC.node(ind);
        F(ind) = F(ind) + sparse( 1:Nbc, 1, val, Nbc, 1 );
    end
    
    % Neumann Boundary Conditions in moment
    ind = find(model.BC.type == 'M');
    Nbc = length(ind);
    if Nbc>0
        val = model.BC.value(ind);
        ind = model.BC.node(ind) + Nn;
        F(ind) = F(ind) + sparse( 1:Nbc, 1, val, Nbc, 1 );
    end

    % construct new matrices
    K = [ K BC; BC' sparse(size(BC,2),size(BC,2)) ];
    F = [ F; FBC ]; 
    
end

