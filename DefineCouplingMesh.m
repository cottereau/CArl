function [mesh,free12,free] = DefineCouplingMesh( LSet, m1, m2 )
% DefineCoupling Create the different subdomain meshes used in the
% Arlequin problem
%
% cplstruct = DefineCoupling( cplstruct, mesh1, mesh2 ) completes the
% cplstruct fields 'free{1}', 'free{2}','free12' and 'coupling'. All these 
% new fields are elements of TRI6 class and describe the different 
% subdomains considered in the computation. The input variables are:
%
%  cplstruct: structured array containing the field 'levelSet' describing
%             the interface of the coupling zone as an element of the
%             levelSet class.
%
%  mesh1,mesh2: TRI6 class variables describing the support domains of 
%               models 1 and 2.
%
% see also: TRI6, levelSet

% R. Cottereau 04/2013
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% default meshes
mesh1 = m1.mesh;
mesh2 = m2.mesh;

% choose meshes for incompatibles models
if strcmp( m1.code, 'Beam')&&strcmp( m2.code, 'FE2D')
    mesh1 = m1.virtualMesh2D;
elseif strcmp( m1.code, 'FE2D')&&strcmp( m2.code, 'Beam')
    mesh2 = m2.virtualMesh2D;
end

% definition of a level-set from boundary of the mesh
bnd1 = freeBoundary(mesh1);
bnd2 = freeBoundary(mesh2);

% definition of coupling area
mesh = intersection( intersection( bnd1, bnd2 ), LSet );

% determination of common free area
free12 = complement( intersection( bnd1, bnd2 ), mesh );

% definition of free areas
free{1} = complement( bnd1, bnd2 );
free{2} = complement( bnd2, bnd1 );
