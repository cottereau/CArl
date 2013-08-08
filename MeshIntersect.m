function [ Int, Rep ] = MeshIntersect( mdl1, mdl2, LSet )
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

% coupling zone for representation purposes (including elements of the
% original meshes that are cut by the level-set)
[meshr1,ind1] = subSet( mdl1.mesh, elementsInBoundary(mdl1.mesh,LSet,false) );
[meshr2,ind2] = subSet( mdl2.mesh, elementsInBoundary(mdl2.mesh,LSet,false) );

% coupling mesh for integration purposes
Int.mesh = bounded( mergeMeshes( meshr1, meshr2 ), LSet );

% compute passage matrix from integration mesh to model meshes
Rep{1} = integration2Model( mdl1, meshr1, Int.mesh, ind1 );
Rep{2} = integration2Model( mdl2, meshr2, Int.mesh, ind2 );
 
% passage matrix from integration mesh to model meshes
function R = integration2Model( mdl, m, mi, ind )
% compute passage matrix in terms of nodes of the mesh
[ Mx, My, Mval ] = XR2XI( m, mi );
M = sparse( ind(Mx), My, Mval, mdl.mesh.Nn, mi.Nn);
% compute passage matrix in terms of DOFs
switch mdl.code
    
    % elastic case
    case {'FE2D'}
        N = 2*size(M);
        [ x, y, M ] = find(M);
        M = sparse( (x-1)*2+1, (y-1)*2+1, M, N(1), N(2) ) + ...
            sparse( x*2,       y*2,       M, N(1), N(2) );

        
    % beam case
    case {'Beam'}
        N = size(M);
        [ x, y, M ] = find(M);
        M = sparse( [x;x+N(1)], [y;y+N(2)], [M;M], 2*N(1), 2*N(2) ); 
        
end
% output
R = struct( 'mesh', m, 'M', M );

