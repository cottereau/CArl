function [ K, F ] = StiffnessMatrixBeam( model, alpha )
% STIFFNESSMATRIXBEAM to construct the stiffness matrix corresponding to a
% Timochenko beam problem by calling an external code
%
% syntax: [K,F] = StiffnessMatrixBeam( model, alpha )
%
%    model: structured array containing the fields 'X', 'T', 'lambda', 'mu'
%           'load', 
%    mesh    : mesh structure [INT3 or TRI6 object]
%
%    C: the output matrix is in sparse format
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr
%
% uses routines by Talal Rahman and Jan Valdman downloadable at 
% http://www.mathworks.in/matlabcentral/fileexchange/
% 27826-fast-assembly-of-stiffness-and-matrices-in-finite-element-method 

        warning('varying alpha not implemented yet for beams')

% geometrical parameters
h = model.L;
I = h^3/12;
E = model.young;
AG = h * 1 * model.young/2/(1+model.poisson);
C = [ E*I 0; 0 AG];
P = 0; % loading (to be modified in the future)

% mesh
X = model.X;
T = model.T;
Nn = size(X,1);
Ne = size(model.T,1);
Nd = 2*Nn;

% compute stiffness matrix
[ K, F ] = formStiffnessMassTimoshenkoBeam(Nd,Ne,T,Nn,X,C,P,1,I,h);
K = sparse(K);
F = sparse(F);
