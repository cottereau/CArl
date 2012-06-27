function [ x, y, K, z, F ] = StiffnessMatrixHomeFE( model )
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
T = model.mesh.Triangulation;
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
F = zeros( Nt, 1 );
z = reshape( T', Nt, 1 );

% gauss weights and shape functions
[ gaussX, gaussW ] = simplexquad( d+2, d );
[ N, Nxi, Neta ] = shapeFunctions( d-1, nnode, gaussX );

% loop on elements - case without load
if all( model.load(:)==0 )
    for i1 = 1:Ne
        Te = T( i1, : );
        Xe = X( Te, : );
        Ee = N * E( i1, : )';
        Ke = elementStiffnessMatrixHomeFENoLoad( Xe, Ee, nnode, ...
            gaussX, gaussW, Nxi, Neta );
        K( (i1-1)*Nn2 + (1:Nn2) ) = Ke( : );
    end
% loop on elements - case with load
else
    for i1 = 1:Ne
        Te = T( i1, : );
        Xe = X( Te, : );
        Ee = N * E( i1, : )';
        Fe = N * model.load( i1, : )';
        [ Ke, fe ] = elementStiffnessMatrixHomeFE( Xe, Ee, nnode, ...
            gaussX, gaussW, N, Nxi, Neta, Fe );
        K( (i1-1)*Nn2 + (1:Nn2) ) = Ke( : );
        F( (i1-1)*nnode + (1:nnode) ) = fe;
    end
end    

% add boundary conditions
if ~isempty( model.BC )
    % Dirichlet Boundary Conditions
    ind = find( model.BC.type == 'U' );
    Nbc = length( ind );
    Nx = max(x);
    Nf = size(model.load,3);
    x = [ x ; Nx+(1:Nbc)'; model.BC.nodes(ind)' ];
    y = [ y ; model.BC.nodes(ind)'; Nx+(1:Nbc)' ];
    K = [ K ; ones( 2*Nbc, 1 ) ];
    z = [ z; reshape(repmat(Nx+(1:Nbc)',[1 Nf]),Nbc*Nf,1) ];
    F = [ F; reshape(repmat(model.BC.value(ind)',[1 Nf]),Nbc*Nf,1) ];
    % Neumann Boundary Conditions
    ind = find( model.BC.type == 'F' );
    Nbc = length( ind );
    z = [ z; reshape(repmat(model.BC.nodes(ind)',[1 Nf]),Nbc*Nf,1) ];
    F = [ F; reshape(repmat(model.BC.value(ind)',[1 Nf]),Nbc*Nf,1) ];
end
%==========================================================================
% ELEMENTSTIFFNESSMATRIXHOMEFENOLOAD
%==========================================================================
function Ke = elementStiffnessMatrixHomeFENoLoad( Xe, Epg, numberOfNodes, ...
                           pospg, pespg, Nxi, Neta ) 
% Ke = elementStiffnessMatrixHomeFENoLoad( Xe, numberOfNodes, pospg, pespg, ...
%   Nxi, Neta ) 
%
% creates an elemental matrix
%
% INPUT
%   Xe             nodal coords
%   numberOfNodes  number of element nodes
%   pospg,pespg   position and weigth of gauss points 
%   N,Nxi,Neta    shape functions and its derivetives
%
% OUTPUT
%   K
%

d = size(Xe,2);
numberOfGaussPoints = size( pospg, 1 );

Ke = zeros( numberOfNodes );

if d==1
    for igaus = 1:numberOfGaussPoints
        dN = Nxi(igaus,:);
        jacob = dN*Xe ;
        dvolu = pespg(igaus) * det( jacob );
        res = jacob \ dN;
        Nx = res(1,:);
        Ke = Ke + ( Nx'*Nx ) * dvolu * Epg(igaus);
    end

elseif d==2
    for igaus = 1:numberOfGaussPoints
        dN = [ Nxi(igaus,:) ; Neta(igaus,:) ];
        jacob = dN*Xe ;
        dvolu = pespg(igaus) * det( jacob );
        res = jacob \ dN;
        Nx = res(1,:);
        Ny = res(2,:);
        Ke = Ke + ( Nx'*Nx + Ny'*Ny ) * dvolu * Epg(igaus);
    end
end

%==========================================================================
% ELEMENTSTIFFNESSMATRIXHOMEFE
%==========================================================================
function [Ke,fe] = elementStiffnessMatrixHomeFE( Xe, Epg, numberOfNodes, ...
                           pospg, pespg, N, Nxi, Neta, load ) 
% [Ke,fe] = elementStiffnessMatrixHomeFE( Xe, numberOfNodes, pospg, pespg, ...
%   N, Nxi, Neta ) 
%
% creates an elemental matrix
%
% INPUT
%   Xe             nodal coords
%   numberOfNodes  number of element nodes
%   pospg,pespg   position and weigth of gauss points 
%   N,Nxi,Neta    shape functions and its derivetives
%
% OUTPUT
%   K
%   f
%

d = size(Xe,2);
numberOfGaussPoints = size( pospg, 1 );

Ke = zeros( numberOfNodes );
fe = zeros( numberOfNodes, 1 );

if d==1
    for igaus = 1:numberOfGaussPoints
        dN = Nxi(igaus,:);
        jacob = dN*Xe ;
        dvolu = pespg(igaus) * det( jacob );
        res = jacob \ dN;
        Nx = res(1,:);
        Ke = Ke + ( Nx'*Nx ) * dvolu * Epg(igaus);
        fe = fe + N(igaus,:)' * load(igaus) * dvolu;
    end

elseif d==2
    for igaus = 1:numberOfGaussPoints
        dN = [ Nxi(igaus,:) ; Neta(igaus,:) ];
        jacob = dN*Xe ;
        dvolu = pespg(igaus) * det( jacob );
        res = jacob \ dN;
        Nx = res(1,:);
        Ny = res(2,:);
        Ke = Ke + ( Nx'*Nx + Ny'*Ny ) * dvolu * Epg(igaus);
        fe = fe + N(igaus,:)' * load(igaus) * dvolu;
    end
end



