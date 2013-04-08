function [mesh,free12,free] = DefineCouplingMesh( LSet, mesh1, mesh2 )
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

% definition of a level-set from boundary of the mesh
[Tbnd,Xbnd] = freeBoundary(mesh1);
bnd1 = levelSet( Tbnd, Xbnd );
[Tbnd,Xbnd] = freeBoundary(mesh2);
bnd2 = levelSet( Tbnd, Xbnd );

% definition of coupling area
mesh = LSintersect( LSintersect( bnd1, bnd2 ), LSet );

% determination of common free area
free12 = LScomplement( LSintersect( bnd1, bnd2 ), mesh );

% definition of free areas
free{1} = LScomplement( bnd1, bnd2 );
free{2} = LScomplement( bnd2, bnd1 );