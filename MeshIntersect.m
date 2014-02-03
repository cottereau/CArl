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

% compute passage matrix in terms of DOFs
switch mdl.code
    
    % elastic case
    case {'FE2D'}
        % beam-solid coupling case
        if max(My)>mi.Nn
            M = sparse( ind(Mx), My, Mval );
            Mu = M(:,1:mi.Nn);
            Mr = M(:,mi.Nn+1:end);
            [ Mx, My, Mval ] = find(Mu);
            N = 2*[mdl.mesh.Nn  mi.Nn];
            Mu = sparse( (Mx-1)*2+1, My,       Mval, N(1), N(2) ) + ...
                 sparse( Mx*2,       My+mi.Nn, Mval, N(1), N(2) );
            [ Mx, My, Mval ] = find(Mr);
            Mr = sparse( (Mx-1)*2+1, My,       Mval, N(1), mi.Nn );
            M = [Mu Mr];
            
        % classical solid-solid coupling
        else
            N = 2*[mdl.mesh.Nn  mi.Nn];
            M = sparse( (ind(Mx)-1)*2+1, (My-1)*2+1, Mval, N(1), N(2) ) ...
              + sparse( ind(Mx)*2,       My*2,       Mval, N(1), N(2) );
        end
        
    % beam case
    case {'Beam'}
        M = sparse( ind(Mx), My, Mval, mdl.mesh.Nn, mi.Nn);
        N = size(M);
        [ x, y, M ] = find(M);
        M = sparse( [x;x+N(1);x+2*N(1)], ...
                    [y;y+N(2);y+2*N(2)], [M;M;M], 3*N(1), 3*N(2) );
        M = M * mdl.Beam.h;
      
    % other cases
    otherwise
        M = sparse( ind(Mx), My, Mval, mdl.mesh.Nn, mi.Nn);
        
end
% output
R = struct( 'mesh', m, 'M', M );

