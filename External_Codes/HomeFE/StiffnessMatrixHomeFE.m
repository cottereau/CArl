function [ x, y, K, z, F, k ] = StiffnessMatrixHomeFE( model )
% STIFFNESSMATRIXHOMEFE to construct the basic stiffness matrix and force 
% vector by calling a home-made FE code
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
Xf = cell( 1, Ne );
Nf = zeros( 1, Ne );
z = cell( 1, Ne );
for i1 = 1:Ne
    Xf{i1} = find( any( model.load(T(i1,:),:)~=0, 1 ) );
    Nf(i1) = length(Xf{i1})*nnode;
    z{i1} = repmat( T(i1,:), [ 1 length(Xf{i1})] );
end
cNf = [0 cumsum(Nf) ];
Nft = cNf(end);
k = reshape( repmat( cat( 2, Xf{:} ), [nnode 1] ), Nft, 1 );
z = cat( 2, z{:} )';
F = zeros( Nft, 1 );   

% gauss weights and shape functions
[ gaussX, gaussW ] = simplexquad( d+2, d );
[ N, Nxi, Neta ] = shapeFunction( d-1, nnode, gaussX );

% loop on MonteCarlo trials
%Nmc = size(E,2);
%K = cell(Nmc,1);
%for i0 = 1:Nmc
K = zeros( Nt2, 1 );
%    K{i0} = zeros( Nt2, 1 );

    % loop on elements
    for i1 = 1:Ne
        Te = T( i1, : );
        Xe = X( Te, : );
        Ee = E( Te, : );
        Fe = N * model.load( Te, Xf{i1} );
        [ Ke, fe ] = elementStiffnessMatrixHomeFE( Xe, Ee, nnode, ...
            gaussX, gaussW, N, Nxi, Neta, Fe );
%        K{i0}( (i1-1)*Nn2 + (1:Nn2) ) = Ke( : );
        K( (i1-1)*Nn2 + (1:Nn2) ) = Ke( : );
        F( cNf(i1) + (1:Nf(i1)) ) = fe( : );
    end

%end

% add boundary conditions
if ~isempty( model.BC )
    Nbc = length( model.BC.type );
    disp('warning: only dirichlet boundary conditions implemented');
    Nx = max(x);
    Nf = size(model.load,2);
    x = [ x ; Nx+(1:Nbc)'; model.BC.nodes' ];
    y = [ y ; model.BC.nodes'; Nx+(1:Nbc)' ];
    K = [ K ; ones( 2*Nbc, 1 ) ];
%    for i0 = 1:Nmc
%        K{i0} = [ K{i0} ; ones( 2*Nbc, 1 ) ];
%    end
    z = [ z; reshape(repmat(Nx+(1:Nbc)',[1 Nf]),Nbc*Nf,1) ];
    k = [ k; reshape(repmat(1:Nf,[Nbc 1]),Nbc*Nf,1) ];
    F = [ F; reshape(repmat(model.BC.value',[1 Nf]),Nbc*Nf,1) ];
end

% simplify when there is no Monte Carlo trial
%if Nmc==1
%    K = K{1};
%end
