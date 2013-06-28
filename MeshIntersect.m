function [ Int, Rep] = MeshIntersect( mesh1, mesh2, LSet,code1,code2 )
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

meshi = MergeMeshes( meshr1, meshr2, LSet );
% compute the passage matrices in terms of nodes
% = get the values of the basis functions in the
% representation meshes at the nodes of the integration
% mesh
[ Mx1, My1, Mval1 ] = XR2XI( meshr1.tri3, meshi.tri3 );
Ni = size(meshi.tri3.X,1);

[ Mx2, My2, Mval2 ] = XR2XI( meshr2.tri3, meshi.tri3 );

M11 = sparse(Xrg1(Mx1),My1,Mval1,N1,Ni);
M22 = sparse(Xrg2(Mx2),My2,Mval2,N2,Ni);

if (strcmp(code1,code2)==0)
    indy0m1 = find(meshi.tri3.X(:,2)==0);
    indx0m1 = find(meshr1.tri3.X(:,2)==0);
    Mtemp1 = sparse(Xrg1(Mx1),My1,Mval1,N1,Ni);
    M1 = Mtemp1(indx0m1,indy0m1);
    
    indy0m2 = find(meshi.tri3.X(:,2)==0);
    indx0m2 = find(meshr2.tri3.X(:,2)==0);
    Mtemp2 = sparse(Xrg2(Mx2),My2,Mval2,N2,Ni);
    M2 = Mtemp2(indx0m2,indy0m2);
    
    M11b = zeros(length(indx0m1),size(Mtemp1,2));%Mtemp(indx0m1,:);
    M22b = zeros(length(indx0m2),size(Mtemp2,2));%Mtemp(indx0m1,:);
    pti = meshi.tri3.X;
    pti2 = mesh1.tri3.X;
    pti3 = mesh2.tri3.X;
    for ijk=1:size(meshi.tri3.X,1)
        aa = pti2(indx0m1,1)-pti(ijk,1);
        [val1 near1] = min(abs(aa));
        val1 = (val1);
        aa(near1)=1000000000;
        [val2 near2] = min(abs(aa));
        val2 = abs(val2);
        col = ijk;
        M11b([(near1) (near2)],col) = 1-[val1;val2]./abs(diff(pti2(indx0m1([near1;near2]),1)));
        
        aa = pti3(indx0m2,1)-pti(ijk,1);
        [val1 near1] = min(abs(aa));
        val1 = (val1);
        aa(near1)=1000000000;
        [val2 near2] = min(abs(aa));
        val2 = abs(val2);
        col = ijk;
        M22b([(near1) (near2)],col) = 1-[val1;val2]./abs(diff(pti3(indx0m2([near1;near2]),1)));
    end
    M11b=sparse(M11b);
    M22b=sparse(M22b);
    
    
end

% output
Int.mesh = meshi;
Rep{1} = struct( 'mesh', meshr1, ...
    'M',  M11 );
Rep{2} = struct( 'mesh', meshr2, ...
    'M', M22 );
if (strcmp(code1,code2)==0)
    Rep{1}.Mbeam=M1;
    Rep{2}.Mbeam=M2;
    Rep{1}.Mbeam2D=M11b;
    Rep{2}.Mbeam2D=M22b;
end
%==========================================================================
function mesh = MergeMeshes( mesh1, mesh2, LSet )
N = size(mesh1.tri3.X,1);
X = [ mesh1.tri3.X ; mesh2.tri3.X ];
T = [ mesh1.tri3.Triangulation ; mesh2.tri3.Triangulation+N ];
C = [T(:,1:2); T(:,2:3); T(:,[3 1])];
for i1=1:LSet.N
    C = [C;LSet.T{i1}+size(X,1)];
    X = [X;LSet.X{i1}];
end
[~,ind] = unique( sort(C,2) ,'rows', 'first');
C = C(ind,:);
ldt = clean( levelSet( C, X ) );
mesh = DelaunayTri( ldt.X{1}, ldt.T{1} );
mesh = TRI6( mesh.Triangulation, mesh.X );
mesh = subSet( mesh, elementsInBoundary(mesh,LSet,true) );

%==========================================================================
function [ Mx, My, Mval ] = XR2XI( meshr, meshi )
% for each node in meshi, find the nodes that are inside the elements
% that touch it, and the value of the linear FE basis function centered on
% Xr (=local barycentric coordinate of that node in the element)
% return a matrix in sparse format
Ni = size(meshi.X,1);
indx = zeros(Ni,1);
gerr = 1e-9;

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
