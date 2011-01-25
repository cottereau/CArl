function [ x, y, K, z, F, k ] = StiffnessMatrixHomeFE( model )
% STIFFNESSMATRIXHOMEFE to construct the basic stiffness matrix and force 
% vector by calling a home-made FE code in 1D
%
% syntax: [ x, y, K, z, F ] = StiffnessMatrixHomeFE( model )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       - 'mesh': .X = matrix of coordinates of nodes [Nn*d matrix]
%                 .T = connectivity matrix [Ne*e matrix]
%       - 'property': vector of mechanical property [Ne*e matrix] of the
%                     value of the property at each node of each element
%       - 'BC': in 1D: .type = list of 'U' and 'G' for displacement and
%                              force imposed
%                      .nodes = list of nodes on which the BC is imposed
%                      .value = value of the BC imposed
%       - 'load': vector of load per element [Ne*e matrix] of the value of
%                 the load at each node of each element
%
%  output: the format is that of sparse matrices. The matrix of stiffness
%          and the vector of force are such that, schematically:
%               Stiffness( x, y ) = K
%               Force( z ) = F
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% description of the model
X = model.mesh.X;
T = model.mesh.T;
E = model.property;

% constants
[ Ne nnode ] = size( T );
Nn2 = nnode^2;
Nt2 = Ne*Nn2;
d = size( X, 2 );
Nt = Ne*nnode;

% initializations
x = reshape( repmat( T, [1 nnode] )', Nt2, 1 );
y = reshape( repmat( reshape(T',Nt,1)', [nnode 1] ), Nt2, 1 );

% initialization of force matrix
K = zeros( Nt2, 1 );
Nf = size(model.load,3);
F = zeros( Nt, Nf );
z = reshape( T', Nt, 1 );
[z,k] = ndgrid( z, (1:Nf)');
z = z(:);
k = k(:);

% gauss weights and shape functions
[ gaussX, gaussW ] = simplexquad( d+2, d );
[ N, Nxi, Neta ] = shapeFunction( d-1, nnode, gaussX );

% loop on elements
for i1 = 1:Ne
    Te = T( i1, : );
    Xe = X( Te, : );
    Ee = N * E( i1, : )';
    Fe = N * reshape(model.load( i1, :, : ),nnode,Nf);
    [ Ke, fe ] = elementStiffnessMatrixHomeFE( Xe, Ee, nnode, ...
                                        gaussX, gaussW, N, Nxi, Neta, Fe );
    K( (i1-1)*Nn2 + (1:Nn2) ) = Ke( : );
    F( (i1-1)*nnode + (1:nnode),: ) = fe;
end
F = F(:);

% add boundary conditions
if ~isempty( model.BC )
    Nbc = length( model.BC.type );
    disp('warning: only dirichlet boundary conditions implemented');
    Nx = max(x);
    Nf = size(model.load,2);
    x = [ x ; Nx+(1:Nbc)'; model.BC.nodes' ];
    y = [ y ; model.BC.nodes'; Nx+(1:Nbc)' ];
    K = [ K ; ones( 2*Nbc, 1 ) ];
    z = [ z; reshape(repmat(Nx+(1:Nbc)',[1 Nf]),Nbc*Nf,1) ];
    k = [ k; reshape(repmat(1:Nf,[Nbc 1]),Nbc*Nf,1) ];
    F = [ F; reshape(repmat(model.BC.value',[1 Nf]),Nbc*Nf,1) ];
end
