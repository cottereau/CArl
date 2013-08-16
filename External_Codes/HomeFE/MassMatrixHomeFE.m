function M = MassMatrixHomeFE( model )
% MASSMATRIXHOMEFE to construct the basic mass matrix by calling a 
% home-made FE code
%
% syntax: M = MassMatrixHomeFE( model )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       - 'mesh': .X = matrix of coordinates of nodes [Nn*d matrix]
%                 .T = connectivity matrix [Ne*e matrix]
%
%  output matrix is sparse
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2012

% description of the model
X = model.mesh.X;
T = model.mesh.T;

% constants
[ Ne, nnode ] = size( T );
Nn2 = nnode^2;
Nt2 = Ne*Nn2;
d = size( X, 2 );
Nt = Ne*nnode;

% initializations
x = reshape( repmat( T, [1 nnode] )', Nt2, 1 );
y = reshape( repmat( reshape(T',Nt,1)', [nnode 1] ), Nt2, 1 );
M = zeros( Nt2, 1 );

% gauss weights and shape functions
[ gaussX, gaussW ] = simplexquad( d+2, d );
[ N, Nxi, Neta ] = shapeFunctions( d-1, nnode, gaussX );

% loop on elements
for i1 = 1:Ne
    Te = T( i1, : );
    Xe = X( Te, : );
    Me = elementMassMatrixHomeFE( Xe, nnode, ...
                                        gaussX, gaussW, N, Nxi, Neta );
    M( (i1-1)*Nn2 + (1:Nn2) ) = Me( : );
end

% sparse format
M = sparse( x, y, M );

%==========================================================================
% ELEMENTMASSMATRIXHOMEFE
%==========================================================================
function Me = elementMassMatrixHomeFE( Xe, numberOfNodes, ...
                           pospg, pespg, N, Nxi, Neta ) 
% Me = elementMassMatrixHomeFE( Xe, numberOfNodes, pospg, pespg, ...
%   N, Nxi, Neta ) 
%
% creates an elemental mass matrix
%
% INPUT
%   Xe             nodal coords
%   numberOfNodes  number of element nodes
%   pospg,pespg   position and weigth of gauss points 
%   N,Nxi,Neta    shape functions and its derivetives
%
% OUTPUT
%   M
%

% constant
d = size(Xe,2);
numberOfGaussPoints = size( pospg, 1 );

% initialization
Me = zeros( numberOfNodes, numberOfNodes );

if d==1
    for igaus = 1:numberOfGaussPoints
        dN = Nxi(igaus,:);
        jacob = dN*Xe ;
        dvolu = pespg(igaus) * det( jacob );
        Me = Me + N(igaus,:)' * N(igaus,:) * dvolu;
    end

elseif d==2
    for igaus = 1:numberOfGaussPoints
        dN = [ Nxi(igaus,:) ; Neta(igaus,:) ];
        jacob = dN*Xe ;
        dvolu = pespg(igaus) * det( jacob );
        Me = Me + N(igaus,:)' * N(igaus,:) * dvolu;
    end
end


