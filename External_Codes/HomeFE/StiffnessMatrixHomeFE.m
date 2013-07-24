function [ K, F ] = StiffnessMatrixHomeFE( model )
% STIFFNESSMATRIXHOMEFE to construct the basic stiffness matrix and force
% vector by calling a home-made FE code in 1D
%
% syntax: [ K, F ] = StiffnessMatrixHomeFE( model )
%
%  model: cell of structured arrays containing the information relative to
%         the different models, in particular:
%       - 'mesh': .X = matrix of coordinates of nodes [Nn*d matrix]
%                 .Triangulation = connectivity matrix [Ne*e matrix]
%       - 'property': vector of mechanical property [Ne*e matrix] of the
%                     value of the property at each node of each element
%       - 'BC': in 1D: .type = list of 'U' and 'G' for displacement and
%                              force imposed
%                      .nodes = list of nodes on which the BC is imposed
%                      .value = value of the BC imposed
%       - 'load': vector of load per element [Ne*e matrix] of the value of
%                 the load at each node of each element
%
%  K and F are sparse matrices
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% description of the model
X = model.mesh.X;
d = size(X,2);
T = model.mesh.T;
[ Ne, nnode ] = size( T );
E = model.property;
if ~isfield(model,'load')
    load = [];
else 
    load = model.load;
end

% initializations
x = reshape( repmat( T, [nnode 1] ), Ne*nnode^2, 1 );
y = reshape( repmat( T, [1 nnode] ), Ne*nnode^2, 1 );
z = reshape( T', Ne*nnode, 1 );

% gauss weights and shape functions
[ xg, wg ] = simplexquad( 2, d );
[ N, Nxi, Neta ] = shapeFunctions( d-1, nnode, xg );

% constructing stiffness matrix without boundary conditions
[ K, F ] = BaseStiffness( X, T, E, wg, N, Nxi, Neta, load );

% add boundary conditions
if isfield(model,'BC')&&~isempty( model.BC )
    
    % Dirichlet Boundary Conditions
    ind = find( model.BC.type == 'U' );
    Nbc = length( ind );
    Nx = max(x);
    Nf = size(load,3);
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
K = sparse( x, y, K );
F = sparse( z, 1, F );

%==========================================================================
% BASESTIFFNESSMATRIXHOMEFE
%==========================================================================
function [K,F] = BaseStiffness( X, T, E, pespg, N, Nxi, Neta, load )

% constants
d = size(X,2);
[Ne,nnode] = size(T);

% initialization
K = zeros( Ne*nnode, nnode );
F = zeros( nnode, Ne );

% loop on virtual fields
for i1 = 1:nnode
    
    % initialization
    lload = (i1==1) && ~(isempty(load)||all(load(:)==0));
    Ke = zeros( nnode, Ne );
    
    % loop on gauss points
    for ipg = 1:length(pespg)
        
        if d==1
            dN = Nxi(ipg,:);
            jac = dN*reshape( X(T'), nnode, Ne );
            dvol = ((N(ipg,:)*E')./jac)*pespg(ipg);
            Ke = Ke + dN'*dvol*dN(i1);
            
        elseif d==2
            dN = [Nxi(ipg,:);Neta(ipg,:)];
            jacx = dN*reshape( X(T',1), nnode, Ne );
            jacy = dN*reshape( X(T',2), nnode, Ne );
            jac = jacx(1,:).*jacy(2,:)-jacx(2,:).*jacy(1,:);
            Nx = dN(1,:)'*jacy(2,:)-dN(2,:)'*jacy(1,:);
            Ny = dN(2,:)'*jacx(1,:)-dN(1,:)'*jacx(2,:);
            dvol = ((pespg(ipg)*N(ipg,:))*E')./jac;
            Ke = Ke + repmat(Nx(i1,:).*dvol,[nnode 1]).*Nx ...
                + repmat(Ny(i1,:).*dvol,[nnode 1]).*Ny;
        end
        
        if lload
            F = F +  N(ipg,:)'*(N(ipg,:)*load'.*jac*pespg(ipg));
        end
        
    end
    % storing value of stiffness
    K( :, i1 ) = reshape(Ke',nnode*Ne,1);

end

% reshaping
K = K(:);
F = reshape(F,nnode*Ne,1);
