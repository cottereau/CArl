function [ x, y, K, z, F, k ] = StiffnessMatrixZacFE( model )
% STIFFNESSMATRIXZacFE to construct the basic stiffness matrix and force 
% vector by calling a home-made FE code using the finite element toolboox
% of Rasmus Anthin
%
% syntax: [ x, y, K, z, F ] = StiffnessMatrixHomeFE( model )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       - 'mesh': .X = matrix of coordinates of nodes [Nn*d matrix]
%                 .T = connectivity matrix [Ne*e matrix]
%       - 'property': vector of mechanical property [Ne*1 vector]
%       - 'BC': in 1D: .type = list of 'U' and 'G' for displacement and
%                              force imposed
%                      .nodes = list of nodes on which the BC is imposed
%                      .value = value of the BC imposed
%       - 'load': vector of load per element [Ne*1 vector]
%
%  output: the format is that of sparse matrices. The matrix of stiffness
%          and the vector of force are such that, schematically:
%               Stiffness( x, y ) = K
%               Force( z ) = F
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% C. Zaccardi 04/2010

% constants
X = model.mesh.X ;
T = model.mesh.T ;
[ Nn d ] = size( X ) ;
E = model.property ;
Fb = model.load ;

% initialization
C = zeros(Nn,1) ;

% computation of base matrices
if d==1
    % property has to be constant per element
    Eg = E ;
    Eg(1) = (1/2)*(E(1)+E(2)) ;
    Eg(2:end) = (1/2)*(E(1:end-1)+E(2:end)) ;
    [A,B,b]=genmat1(X',Eg',C',Fb') ;
elseif d==2
    [A,B,b]=genmat2(X,T,E,C,Fb) ;
end
Kint = A+B ;
[ x, y, K ] = find(Kint) ;
[z, k, F ] = find(b) ;

% add boundary conditions
if ~isempty( model.BC )
    Nbc = length( model.BC.type ) ;
    disp('warning: only dirichlet boundary conditions implemented') ;
    Nx = max(x) ;
    Nf = size(model.load,2) ;
    x = [ x ; Nx+(1:Nbc)'; model.BC.nodes' ] ;
    y = [ y ; model.BC.nodes'; Nx+(1:Nbc)' ] ;
    K = [ K ; ones( 2*Nbc, 1 ) ] ;
    z = [ z; reshape(repmat(Nx+(1:Nbc)',[1 Nf]),Nbc*Nf,1) ] ;
    k = [ k; reshape(repmat(1:Nf,[Nbc 1]),Nbc*Nf,1) ] ;
    F = [ F; reshape(repmat(model.BC.value',[1 Nf]),Nbc*Nf,1) ] ;
end
