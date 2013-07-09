function cpl = DefineDomains( cpl, m1, m2 )
% DefineDomains Create the different subdomain meshes used in the
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
    ibeam = 1;
elseif strcmp( m1.code, 'FE2D')&&strcmp( m2.code, 'Beam')
    mesh2 = m2.virtualMesh2D;
    ibeam = 2;
end

% definition of a level-set from boundary of the mesh
bnd1 = freeBoundary(mesh1);
bnd2 = freeBoundary(mesh2);

% definition of coupling area
cpl.domain = intersection( intersection( bnd1, bnd2 ), cpl.levelSet );

% determination of common free area
cpl.free12 = complement( intersection( bnd1, bnd2 ), cpl.domain );

% definition of free areas
cpl.free{1} = complement( bnd1, bnd2 );
cpl.free{2} = complement( bnd2, bnd1 );

% choose coupling mesh for incompatibles models
if strcmp( m1.code, 'Beam')&&strcmp( m2.code, 'FE2D') || ...
                         strcmp( m1.code, 'FE2D')&&strcmp( m2.code, 'Beam')
    cpl.virtual2D = cpl.mesh;
    cpl.mesh = struct( 'X', project1D( cpl.virtual2D )' );
    cpl.mesh.T = reshape( 1:(2*cpl.virtual2D.N), 2, cpl.virtual2D.N )';
    cpl.virtual2Dfree12 = cpl.free12;
    cpl.free12 = struct( 'X', project1D( cpl.virtual2Dfree12 )' );
    cpl.free12.T = reshape( 1:(2*cpl.virtual2Dfree12.N), 2, ...
                                                  cpl.virtual2Dfree12.N )';
    cpl.free{ibeam} = struct( 'X', project1D( cpl.free{ibeam} )' );
    N = size(cpl.free{ibeam}.X,1);
    cpl.free{ibeam}.T = reshape( 1:(2*N), 2, N )';
end

% projecting a 2D domain onto a curve
% we are assuming here that this curve is along x only
function Bnd = project1D( LSet )
    Bnd = zeros( LSet.N, 2 );
    for i1 = 1:LSet.N
        Bnd(i1,:) = [min(LSet.X{i1}(:,1)) max(LSet.X{i1}(:,1))];
    end
