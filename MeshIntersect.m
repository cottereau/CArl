function [ Int, Rep] = MeshIntersect( mesh1, mesh2, LSet1, LSet2 )
% MESHINTERSECT to create the mesh at the intersection between two meshes,
% both in terms of support (for the definition of the interpolation
% functions) and for integration purposes
%
%  syntax: [ci,c1,c2] = MeshIntersect( model1, model2 )
%
%  model1, model2 : structured arrays containing the information relative
%         to the models, in particular:
%       - 'type': describes the type of representation used for the
%                 geometry. This is used for all interpolation purposes,
%                 in particular for the definition of weight functions
%                 Implemented: {'FE' 'discrete'}
%       - 'mesh': array dependent on 'type'.
%
%  ci: mesh used for integration purposes. It is a structured array
%     containing the fields 'X' and 'T'
%  c1,c2: meshes used for representation purposes for each of the two
%         models, and passage between representation and integration
%         meshes. Structured arrays containing the fields
%       - 'Tr' describes the mesh for representation purposes by giving the
%              appropriate indices of model.mesh.T
%       - 'r2i' cell giving the list of elements of the integration mesh Ti
%               that are in each element of the representation mesh Tr
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% % constants
[N1,d] = size(mesh1.X);
N2 = size(mesh2.X,1);

% coupling zone for representation purposes
[meshr1,Xrg1] = DefineCouplingMesh( mesh1, LSet1.int.*LSet1.ext );
[meshr2,Xrg2] = DefineCouplingMesh( mesh2, LSet2.int.*LSet2.ext );

% intersect the two meshes to get a first draft of the integration mesh
meshr = MergeMeshes( meshr1, meshr2 );

% create a mesh for the coupling area
area = LevelSetMesh( LSet1 );

% delimitate meshi by the level sets area
meshi = ElagateMesh( area, meshr );

% compute the passage matrices in terms of nodes
% = get the values of the basis functions in the
% representation meshes at the nodes of the integration
% mesh
[ Mx1, My1, Mval1 ] = XR2XI( meshr1, meshi );
[ Mx2, My2, Mval2 ] = XR2XI( meshr2, meshi );

% output
Int.mesh = meshi;
Ni = size(meshi.X,1);
Rep{1} = struct( 'mesh', meshr1, ...
                 'M', sparse(Xrg1(Mx1),My1,Mval1,N1,Ni) );
Rep{2} = struct( 'mesh', meshr2, ...
                 'M', sparse(Xrg2(Mx2),My2,Mval2,N2,Ni) );

%==========================================================================
function mesh = MergeDoubleNodes( mesh )
% merge nodes that are too close to each other
gerr = 1e-8;
d = size(mesh.X,2);
x = mesh.X(:,1);
T = mesh.Triangulation;

if d==1
    d = abs( x(T(:,1)) - x(T(:,2)) );
    T = T( d > gerr, : );
    [ X, T ] = ReduceMesh( x, T );
    mesh = struct( 'X', X, 'Triangulation', T );
    
elseif d==2
    
    y = mesh.X(:,2);
    
    % compute elements of zero area
    S = polyarea( x(T'), y(T') );
    
    % select candidate nodes for repeated
    ind = unique(T(S<gerr,:));
    [ X1, X2 ] = ndgrid( x(ind), x(ind) );
    [ Y1, Y2 ] = ndgrid( y(ind), y(ind) );
    
    % compute distance between all the points
    d = sqrt((X1-X2).^2+(Y1-Y2).^2) + triu( ones(length(ind)), 0 );
    [ indx, indy ] = find( d < gerr );
    ind = sort( [ind(indx) ind(indy)], 2 );
    
    % merge repeated nodes
    for i1=1:size(ind,1)
        T( T==ind(i1,2) ) = ind(i1,1);
    end
    
    % erase flat elements
    T = T( S>gerr, : );
    
    % erase unnecessary nodes
    [ X, T ] = ReduceMesh( mesh.X, T );
    mesh = TriRep( T, X );
end

%==========================================================================
function mesh = ElagateMesh( area, meshr )
[N,d]=size(area.X);

if d==1
    T = area.Triangulation;
    ind = false( size(meshr.X) );
    for i1 = 1:size(T,1)
        ind = ind | ( meshr.X>=area.X(T(i1,1)) & meshr.X<=area.X(T(i1,2)) );
    end
    T = meshr.Triangulation;
    T = T( all(ismember(T,find(ind)),2), : );
    mesh = struct( 'X', meshr.X, 'Triangulation', T );
    mesh = MergeDoubleNodes( mesh );
    
elseif d==2
    X = [ area.X ; meshr.X ];
    C = [ freeBoundary(area); meshr.Constraints+N ];
    mesh = DelaunayTri( X, C );
    mesh = MergeDoubleNodes( mesh );
    mesh = BoundedMesh( mesh, area );
    mesh = BoundedMesh( mesh, meshr );
    
end

%==========================================================================
function area = LevelSetMesh( LSet )
[N,d] = size(LSet.meshint.X);

if d==1 || N==1
    X = sort( [ LSet.meshint.X(:); LSet.meshext.X(:)]);
    if ~any(isinf(X))
        T = [1 2;3 4];
    else
        X = X(~isinf(X));
        T = [1 2];
    end
    area = struct( 'X', X, 'Triangulation', T );
    
elseif d==2
    X = [ LSet.meshint.X; LSet.meshext.X ];
    int = DelaunayTri( X, LSet.meshint.T );
    ext = DelaunayTri( X, LSet.meshext.T+N );
    area = DelaunayTri( X, [LSet.meshint.T;LSet.meshext.T+N] );
    if all(inOutStatus(int))||all(inOutStatus(ext))
        ind = inOutStatus(area);
    else
        ind = ~inOutStatus(area);
    end
    area = TriRep( area.Triangulation(ind,:), area.X );
    
end

%==========================================================================
function mesh = MergeMeshes( mesh1, mesh2 )
[N,d] = size(mesh1.X);
X = [ mesh1.X ; mesh2.X ];
T = [ mesh1.Triangulation ; mesh2.Triangulation+N ];
if d==1
    T1 = mesh1.Triangulation;
    X = []; T = [];
    ind = [0; find(diff(T1(:,1))>1); size(T1,1) ]+1;
    for i1=1:length(ind)-1
        x1 = mesh1.X(T1(ind(i1),1):T1(ind(i1+1)-1,2));
        ind2 = mesh2.X>=x1(1) & mesh2.X<=x1(end);
        x1 = unique( [x1; mesh2.X(ind2)] );
        N1 = length(x1);
        T = [T; [(1:N1-1)' (2:N1)']+size(X,1)];
        X = [X; x1];
    end
    mesh = struct( 'X', X, 'Triangulation', T );
elseif d==2
    mesh = DelaunayTri( X, [T(:,1:2); T(:,2:3); T(:,[3 1])] );
end

%==========================================================================
function mesh = BoundedMesh( mesh, Bnd )
% to extract the submesh bounded by the free boundary of other meshes
indT = inpoly( mesh.incenters, Bnd.X, Bnd.freeBoundary );
[ X, T ] = ReduceMesh( mesh.X, mesh.Triangulation(indT,:) );
mesh = TriRep( T, X );

%==========================================================================
function [mesh,indX] = DefineCouplingMesh( mesh, LSp )
% to extract the submesh where the product of level set
% functions is positive (ie the coupling zone)
gerr = 1e-9;
d = size(mesh.X,2);
T = mesh.Triangulation;
%indX = (abs(LSp) <= gerr) | (LSp >= gerr) ;
ind0 = abs(LSp) <= gerr;
indT = all( (LSp(T)>=gerr) | ind0(T), 2 );
[ X, T, indX ] = ReduceMesh( mesh.X, T(indT,:) );
if d==1
    mesh = struct( 'X', X, 'Triangulation', T );
elseif d==2
    mesh = DelaunayTri( X, [T(:,1:2); T(:,2:3); T(:,[1 3])] );
end

%==========================================================================
function [ Mx, My, Mval ] = XR2XI( meshr, meshi )
% for each node in meshi, find the nodes that are inside the elements
% that touch it, and the value of the linear FE basis function centered on
% Xr (=local barycentric coordinate of that node in the element)
% return a matrix in sparse format
d = size(meshr.X,2);
Ni = size(meshi.X,1);

if d==1
    indx = zeros(Ni,1);
    X1 = meshr.X(meshr.Triangulation(:,1));
    for i1 = 1:Ni
        indx(i1) = find( X1<=meshi.X(i1), 1, 'last' );
    end
    indx = meshr.Triangulation(indx,:);
    My = repmat( (1:size(indx,1))', [2 1] );
    Mx = indx(:);
    Mval = (meshi.X-meshr.X(indx(:,1)))./diff(meshr.X(indx),1,2);
    Mval = [1-Mval ; Mval];
    
elseif d==2
    [ indx, Mval ] = pointLocation( meshr, meshi.X );
    indx = meshr.Triangulation(indx,:);
    Mx = indx(:);
    My = repmat( (1:size(Mval,1))', [3 1] );
    Mval = Mval(:);
end
%==========================================================================
function [ X, T, Xn ] = ReduceMesh( X, T )
% REDUCEMESH to get rid of unused nodes and renumber the connectivity
% matrix
%
% syntax: [ X, T, Xn ] = ReduceMesh( X, T )
%
%  T,X: connectivity matrices and nodal coordinates [Ne*3 matrix] and 
%       [Nn*2 matrix]
%  Xn: indices into the input connectivity matrix
Xn = unique( T(:) );
Nn = size(Xn,1);
for i1 = 1:Nn
    T( T(:)==Xn(i1) ) = i1;
end
X = X(Xn,:);

