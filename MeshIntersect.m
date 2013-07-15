function [ Int, Rep ] = MeshIntersect( mesh1, mesh2, LSet )
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

% coupling zone for representation purposes (including elements of the
% original meshes that are cut by the level-set)
[meshr1,ind1] = subSet( mesh1, elementsInBoundary(mesh1,LSet,false) );
[meshr2,ind2] = subSet( mesh2, elementsInBoundary(mesh2,LSet,false) );

% coupling mesh for integration puroses
meshi = MergeMeshes( bounded(meshr1,LSet), bounded(meshr2,LSet), LSet );

% compute the passage matrices in terms of nodes = get the values of the 
% basis functions in the representation meshes at the nodes of the 
% integration mesh
[ Mx1, My1, Mval1 ] = XR2XI( meshr1, meshi );
[ Mx2, My2, Mval2 ] = XR2XI( meshr2, meshi );
keyboard 
M11 = sparse( ind1(Mx1), My1, Mval1, mesh1.Nn, meshi.Nn);
M22 = sparse( ind2(Mx2), My2, Mval2, mesh2.Nn, meshi.Nn);

% output
Int.mesh = meshi;
Rep{1} = struct( 'mesh', meshr1, 'M', M11 );
Rep{2} = struct( 'mesh', meshr2, 'M', M22 );
