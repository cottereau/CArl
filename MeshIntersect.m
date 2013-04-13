function [ Int, Rep] = MeshIntersect( mesh1, mesh2, LSet )
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
N1 = size(mesh1.tri3.X,1);
N2 = size(mesh2.tri3.X,1);

% coupling zone for representation purposes (including elements of the
% original meshes that are cut by the level-set)
[meshr1,Xrg1] = subSet( mesh1, elementsInBoundary(mesh1,LSet,false) );
Xrg1 = Xrg1( ismember(Xrg1,mesh1.ind3v6) );
[meshr2,Xrg2] = subSet( mesh2, elementsInBoundary(mesh2,LSet,false) );
Xrg2 = Xrg2( ismember(Xrg2,mesh2.ind3v6) );

% intersect the two meshes to get a first draft of the integration mesh
meshi = MergeMeshes( meshr1, meshr2, LSet );

% compute the passage matrices in terms of nodes
% = get the values of the basis functions in the
% representation meshes at the nodes of the integration
% mesh
[ Mx1, My1, Mval1 ] = XR2XI( meshr1.tri3, meshi.tri3 );
[ Mx2, My2, Mval2 ] = XR2XI( meshr2.tri3, meshi.tri3 );

% output
Int.mesh = meshi;
Ni = size(meshi.tri3.X,1);
Rep{1} = struct( 'mesh', meshr1, ...
                 'M', sparse(Xrg1(Mx1),My1,Mval1,N1,Ni) );
Rep{2} = struct( 'mesh', meshr2, ...
                 'M', sparse(Xrg2(Mx2),My2,Mval2,N2,Ni) );

%==========================================================================
function mesh = MergeMeshes( mesh1, mesh2, LSet )
[N,d] = size(mesh1.tri3.X);
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
    X = [ mesh1.tri3.X ; mesh2.tri3.X ];
    T = [ mesh1.tri3.Triangulation ; mesh2.tri3.Triangulation+N ];
    C = [T(:,1:2); T(:,2:3); T(:,[3 1])];
    for i1=1:LSet.N
        C = [C;LSet.T{i1}+size(X,1)];
        X = [X;LSet.X{i1}];
    end
    ldt = clean( levelSet( C, X ) );
    mesh = DelaunayTri( ldt.X{1}, ldt.T{1} );
    mesh = TRI6( mesh.Triangulation, mesh.X );
    bnd1 = freeBoundary(mesh1);
    bnd2 = freeBoundary(mesh2);
    ind = inside(LSet,mesh.X);
    ind = elementsInBoundary(mesh,levelSet(bnd1.T{1},bnd1.X{1})) & ...
          elementsInBoundary(mesh,levelSet(bnd2.T{1},bnd2.X{1})) & ...
          all(ind(mesh.T),2);
    mesh = subSet( mesh, ind );
end

% %==========================================================================
% function [mesh,indX] = ReduceCouplingArea( mesh, bndCoupling )
% % to extract the submesh where the product of level set
% % functions is positive or null for at least one node
% %
% % output mesh is a TRI6 (possibly non-convex) mesh
% 
% % % constants
% % d = size(mesh.X,2);
% % T = mesh.T;
% 
% % selection of elements inside the coupling area (positive product of
% % level-sets)
% indT = elementsInBoundary( mesh, bndCoupling, false );
% [mesh,indX] = subSet( mesh, indT );
% 
% % % creation of output structure
% % [ X, T, indX ] = ReduceMesh( mesh.X, T(indT,:) );
% % if d==1
% %     mesh = struct( 'X', X, 'Triangulation', T );
% % elseif d==2
% %     mesh = TRI6( T, X );
% % end

%==========================================================================
function [ Mx, My, Mval ] = XR2XI( meshr, meshi )
% for each node in meshi, find the nodes that are inside the elements
% that touch it, and the value of the linear FE basis function centered on
% Xr (=local barycentric coordinate of that node in the element)
% return a matrix in sparse format
d = size(meshr.X,2);
Ni = size(meshi.X,1);
indx = zeros(Ni,1);
gerr = 1e-9;

if d==1
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
    Nr = meshr.size;
    ey = repmat((1:Nr(1)),[Ni 1]);
    cc = reshape( cartToBary( meshr, ey(:), ...
             repmat(meshi.X,[Nr(1) 1]) ), [Ni Nr(1) 3]);
    for i1 = 1:Ni
        cci = squeeze(cc(i1,:,:));
        indx(i1) = find( all(cci>=-gerr,2) & all(cci<=1+gerr,2), 1, 'first' );
    end
    Mval = cartToBary( meshr, indx, meshi.X );
    My = repmat( (1:size(Mval,1))', [3 1] );
    indx = meshr.Triangulation(indx,:);
    %  alleluiah !!
    ind = abs(Mval(:))>gerr;
    %  alleluiah !!
    Mx = indx(ind);
    Mval = Mval(ind);
    My = My(ind);
end
% %==========================================================================
% function [ X, T, Xn ] = ReduceMesh( X, T )
% % REDUCEMESH to get rid of unused nodes and renumber the connectivity
% % matrix
% %
% % syntax: [ X, T, Xn ] = ReduceMesh( X, T )
% %
% %  T,X: connectivity matrices and nodal coordinates [Ne*3 matrix] and 
% %       [Nn*2 matrix]
% %  Xn: indices into the input connectivity matrix
% Xn = unique( T(:) );
% Nn = size(Xn,1);
% for i1 = 1:Nn
%     T( T(:)==Xn(i1) ) = i1;
% end
% X = X(Xn,:);
